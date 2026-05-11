"""Single-cell RNA-seq tools: barcode validation, QC metrics, UMI analysis."""

import math
import re
from collections import Counter


# ---------------------------------------------------------------------------
# Chemistry definitions
# ---------------------------------------------------------------------------

CHEMISTRY_BARCODE_LENGTH = {
    "10x_v2": 16,
    "10x_v3": 16,
    "10x_v3.1": 16,
    "dropseq": 12,
    "indrop": 8,
}

VALID_BASES = set("ACGT")


# ---------------------------------------------------------------------------
# Cell barcode validator
# ---------------------------------------------------------------------------

def validate_barcodes(
    barcodes: list,
    chemistry: str = "10x_v3",
    whitelist: str | None = None,
) -> dict:
    """
    Validate a list of cell barcodes for format and optionally against a whitelist.

    Checks expected barcode length for the given chemistry, flags barcodes with
    ambiguous bases (N), and — when a whitelist is supplied — reports exact matches
    and Hamming-distance-1 correctable barcodes (the same correction logic used
    by Cell Ranger).

    Args:
        barcodes: List of raw barcode strings (e.g. from R1 reads).
        chemistry: Sequencing chemistry — 10x_v2, 10x_v3, 10x_v3.1, dropseq, indrop.
        whitelist: Optional newline-separated string of valid barcodes. When provided,
                   enables whitelist matching and Hamming-1 correction.
    """
    chemistry = chemistry.lower()
    if chemistry not in CHEMISTRY_BARCODE_LENGTH:
        return {
            "error": f"Unknown chemistry '{chemistry}'. "
                     f"Choose from: {sorted(CHEMISTRY_BARCODE_LENGTH)}"
        }

    expected_len = CHEMISTRY_BARCODE_LENGTH[chemistry]
    wl_set = set(whitelist.split()) if whitelist else None

    wrong_length, has_n, exact_match, correctable, invalid = [], [], [], [], []

    for bc in barcodes:
        bc = bc.strip().upper()

        if len(bc) != expected_len:
            wrong_length.append(bc)
            continue

        if "N" in bc:
            has_n.append(bc)
            continue

        if not VALID_BASES.issuperset(bc):
            invalid.append(bc)
            continue

        if wl_set is None:
            exact_match.append(bc)
            continue

        if bc in wl_set:
            exact_match.append(bc)
        elif _hamming1_correction(bc, wl_set) is not None:
            correctable.append(bc)
        else:
            invalid.append(bc)

    total = len(barcodes)
    result = {
        "chemistry": chemistry,
        "expected_barcode_length": expected_len,
        "total_barcodes": total,
        "wrong_length": len(wrong_length),
        "contains_n": len(has_n),
        "valid_format": len(exact_match) + len(correctable),
        "exact_whitelist_match": len(exact_match) if wl_set else None,
        "hamming1_correctable": len(correctable) if wl_set else None,
        "invalid": len(invalid),
        "pass_rate_pct": round((len(exact_match) + len(correctable)) / total * 100, 2)
        if total > 0 else 0.0,
    }

    if wl_set is None:
        result.pop("exact_whitelist_match")
        result.pop("hamming1_correctable")
        result["note"] = "No whitelist provided — reporting format validity only."

    return result


def _hamming1_correction(barcode: str, whitelist: set) -> str | None:
    """Return a whitelist match within Hamming distance 1, or None."""
    for i, base in enumerate(barcode):
        for alt in VALID_BASES - {base}:
            candidate = barcode[:i] + alt + barcode[i + 1:]
            if candidate in whitelist:
                return candidate
    return None


# ---------------------------------------------------------------------------
# scRNA-seq QC metrics
# ---------------------------------------------------------------------------

def compute_scrna_qc(count_matrix_csv: str, mt_prefix: str = "MT-") -> dict:
    """
    Compute per-cell QC metrics from a genes-x-cells or cells-x-genes CSV count matrix.

    Mirrors the core metrics computed by scanpy.pp.calculate_qc_metrics:
    - total_counts: total UMI counts per cell
    - n_genes_by_counts: number of genes with at least 1 count
    - pct_counts_mt: percentage of counts from mitochondrial genes

    High pct_counts_mt (>20%) typically indicates dying or low-quality cells.
    Low n_genes_by_counts may indicate empty droplets or doublets at extremes.

    Args:
        count_matrix_csv: CSV string with gene names as one axis and cell barcodes
                          as the other. First row and first column are treated as labels.
        mt_prefix: Prefix used to identify mitochondrial genes (default: MT-).
    """
    lines = [l.strip() for l in count_matrix_csv.strip().splitlines() if l.strip()]
    if len(lines) < 2:
        return {"error": "Count matrix must have at least a header row and one data row."}

    header = lines[0].split(",")
    cell_barcodes = header[1:]

    if not cell_barcodes:
        return {"error": "No cell barcodes found in header row."}

    n_cells = len(cell_barcodes)
    total_counts = [0] * n_cells
    n_genes = [0] * n_cells
    mt_counts = [0] * n_cells
    genes_seen = 0

    for line in lines[1:]:
        parts = line.split(",")
        gene = parts[0]
        try:
            counts = [int(x) for x in parts[1:n_cells + 1]]
        except ValueError:
            return {"error": f"Non-integer count value in row for gene '{gene}'."}

        genes_seen += 1
        is_mt = gene.upper().startswith(mt_prefix.upper())

        for i, c in enumerate(counts):
            total_counts[i] += c
            if c > 0:
                n_genes[i] += 1
            if is_mt:
                mt_counts[i] += c

    per_cell = []
    for i, bc in enumerate(cell_barcodes):
        tc = total_counts[i]
        pct_mt = round(mt_counts[i] / tc * 100, 2) if tc > 0 else 0.0
        per_cell.append({
            "barcode": bc,
            "total_counts": tc,
            "n_genes_by_counts": n_genes[i],
            "pct_counts_mt": pct_mt,
        })

    all_tc = total_counts
    all_ng = n_genes
    all_mt = [c["pct_counts_mt"] for c in per_cell]

    return {
        "n_cells": n_cells,
        "n_genes": genes_seen,
        "summary": {
            "median_total_counts": _median(all_tc),
            "median_n_genes": _median(all_ng),
            "median_pct_counts_mt": _median(all_mt),
            "cells_high_mt_pct": sum(1 for x in all_mt if x > 20),
        },
        "per_cell": per_cell,
    }


# ---------------------------------------------------------------------------
# UMI analyzer
# ---------------------------------------------------------------------------

def analyze_umis(umis: list, chemistry: str = "10x_v3") -> dict:
    """
    Analyze a list of UMI sequences for quality, composition, and collision risk.

    UMIs (Unique Molecular Identifiers) are short random sequences (10–12 bp in
    10x Genomics) ligated to each mRNA molecule before PCR. They allow PCR
    duplicates to be collapsed, enabling accurate transcript counting.

    Reports nucleotide bias (deviations from 25% per base indicate synthesis bias),
    N-base rate (sequencing errors), and estimates PCR collision probability using
    the birthday problem approximation — relevant when UMI diversity is low relative
    to molecule count.

    Args:
        umis: List of UMI sequences (e.g. extracted from R1 reads, bases 17-28 for 10x v3).
        chemistry: Used to determine expected UMI length —
                   10x_v3/10x_v3.1 = 12 bp, 10x_v2 = 10 bp, dropseq = 8 bp, indrop = 6 bp.
    """
    UMI_LENGTHS = {"10x_v2": 10, "10x_v3": 12, "10x_v3.1": 12, "dropseq": 8, "indrop": 6}
    chemistry = chemistry.lower()
    if chemistry not in UMI_LENGTHS:
        return {"error": f"Unknown chemistry '{chemistry}'. Choose from: {sorted(UMI_LENGTHS)}"}

    expected_len = UMI_LENGTHS[chemistry]

    if not umis:
        return {"error": "No UMIs provided."}

    total = len(umis)
    wrong_len = sum(1 for u in umis if len(u) != expected_len)
    valid_umis = [u.upper() for u in umis if len(u) == expected_len]

    n_content = sum(u.count("N") for u in valid_umis)
    n_rate = round(n_content / (len(valid_umis) * expected_len) * 100, 3) if valid_umis else 0.0

    # Per-position nucleotide frequency
    pos_freqs = []
    for pos in range(expected_len):
        bases_at_pos = [u[pos] for u in valid_umis if u[pos] in VALID_BASES]
        cnt = Counter(bases_at_pos)
        total_at_pos = sum(cnt.values()) or 1
        pos_freqs.append({b: round(cnt[b] / total_at_pos * 100, 1) for b in "ACGT"})

    # Overall base composition across all positions
    all_bases = "".join(valid_umis)
    total_bases = len(all_bases) or 1
    composition = {b: round(all_bases.count(b) / total_bases * 100, 2) for b in "ACGT"}

    # Unique UMI count and duplication
    unique_umis = len(set(valid_umis))
    duplication_rate = round((1 - unique_umis / len(valid_umis)) * 100, 2) if valid_umis else 0.0

    # Collision probability estimate (birthday problem)
    umi_space = 4 ** expected_len
    collision_prob = round((1 - math.exp(-len(valid_umis) ** 2 / (2 * umi_space))) * 100, 4)

    return {
        "chemistry": chemistry,
        "expected_umi_length": expected_len,
        "total_umis": total,
        "wrong_length": wrong_len,
        "valid_umis": len(valid_umis),
        "unique_umis": unique_umis,
        "duplication_rate_pct": duplication_rate,
        "n_base_rate_pct": n_rate,
        "base_composition": composition,
        "per_position_frequencies": pos_freqs,
        "umi_space": umi_space,
        "collision_probability_pct": collision_prob,
    }


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _median(values: list) -> float:
    if not values:
        return 0.0
    s = sorted(values)
    n = len(s)
    mid = n // 2
    return round((s[mid] + s[~mid]) / 2, 2)
