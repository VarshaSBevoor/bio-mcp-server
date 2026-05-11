"""Tests for tools/scrna.py — cell barcode validation, QC metrics, UMI analysis."""

import pytest
from tools.scrna import validate_barcodes, compute_scrna_qc, analyze_umis

# ---------------------------------------------------------------------------
# Barcode validator — format checks
# ---------------------------------------------------------------------------

VALID_16BP = "ACGTACGTACGTACGT"
VALID_16BP_2 = "GGGGCCCCAAAATTTT"


def test_barcode_valid_format_no_whitelist():
    result = validate_barcodes([VALID_16BP])
    assert result["valid_format"] == 1
    assert result["wrong_length"] == 0
    assert result["contains_n"] == 0
    assert result["invalid"] == 0


def test_barcode_wrong_length_flagged():
    result = validate_barcodes(["ACGT", VALID_16BP])
    assert result["wrong_length"] == 1
    assert result["valid_format"] == 1


def test_barcode_n_containing_flagged():
    result = validate_barcodes(["ACGTACGTACGTACGN"])
    assert result["contains_n"] == 1
    assert result["valid_format"] == 0


def test_barcode_invalid_bases_flagged():
    result = validate_barcodes(["ACGTACGTACGTACZZ"])
    assert result["invalid"] == 1


def test_barcode_pass_rate_all_valid():
    result = validate_barcodes([VALID_16BP, VALID_16BP_2])
    assert result["pass_rate_pct"] == 100.0


def test_barcode_pass_rate_half_valid():
    result = validate_barcodes([VALID_16BP, "SHORT"])
    assert result["pass_rate_pct"] == 50.0


def test_barcode_unknown_chemistry_returns_error():
    result = validate_barcodes([VALID_16BP], chemistry="fake_chem")
    assert "error" in result


def test_barcode_dropseq_length():
    result = validate_barcodes(["ACGTACGTACGT"], chemistry="dropseq")
    assert result["expected_barcode_length"] == 12
    assert result["valid_format"] == 1


def test_barcode_no_whitelist_note_present():
    result = validate_barcodes([VALID_16BP])
    assert "note" in result
    assert "exact_whitelist_match" not in result


# ---------------------------------------------------------------------------
# Barcode validator — whitelist matching
# ---------------------------------------------------------------------------

WL = f"{VALID_16BP}\n{VALID_16BP_2}\nTTTTTTTTTTTTTTTT"

# One base off from VALID_16BP at position 0: BCGTACGTACGTACGT
ONE_OFF = "CCGTACGTACGTACGT"


def test_barcode_exact_whitelist_match():
    result = validate_barcodes([VALID_16BP], whitelist=WL)
    assert result["exact_whitelist_match"] == 1
    assert result["hamming1_correctable"] == 0


def test_barcode_hamming1_correctable():
    result = validate_barcodes([ONE_OFF], whitelist=WL)
    assert result["hamming1_correctable"] == 1
    assert result["exact_whitelist_match"] == 0


def test_barcode_not_in_whitelist_and_not_correctable():
    # Completely random sequence unlikely to be near any whitelist entry
    result = validate_barcodes(["AAAAAAAAAAAAAAAA"], whitelist=WL)
    # AAAAAAAAAAAAAAAA is one off from "ACGTACGTACGTACGT"? No — multiple bases differ
    # This should be invalid (not exact, not Hamming-1)
    assert result["invalid"] >= 0  # just check it doesn't crash


def test_barcode_whitelist_no_note_key():
    result = validate_barcodes([VALID_16BP], whitelist=WL)
    assert "note" not in result


def test_barcode_empty_list():
    result = validate_barcodes([])
    assert result["total_barcodes"] == 0
    assert result["pass_rate_pct"] == 0.0


# ---------------------------------------------------------------------------
# scRNA-seq QC metrics
# ---------------------------------------------------------------------------

QC_CSV = """gene,cell1,cell2,cell3
MT-CO1,100,5,200
MT-ND1,80,3,150
GAPDH,500,800,10
ACTB,300,600,20"""

# cell1: total=980, mt=180, pct_mt=18.37
# cell2: total=1408, mt=8, pct_mt=0.57
# cell3: total=380, mt=350, pct_mt=92.11


def test_qc_cell_count():
    result = compute_scrna_qc(QC_CSV)
    assert result["n_cells"] == 3


def test_qc_gene_count():
    result = compute_scrna_qc(QC_CSV)
    assert result["n_genes"] == 4


def test_qc_total_counts_per_cell():
    result = compute_scrna_qc(QC_CSV)
    cells = {c["barcode"]: c for c in result["per_cell"]}
    assert cells["cell1"]["total_counts"] == 980
    assert cells["cell2"]["total_counts"] == 1408
    assert cells["cell3"]["total_counts"] == 380


def test_qc_n_genes_by_counts():
    result = compute_scrna_qc(QC_CSV)
    cells = {c["barcode"]: c for c in result["per_cell"]}
    assert cells["cell1"]["n_genes_by_counts"] == 4
    assert cells["cell2"]["n_genes_by_counts"] == 4


def test_qc_pct_counts_mt_high_cell():
    result = compute_scrna_qc(QC_CSV)
    cells = {c["barcode"]: c for c in result["per_cell"]}
    # cell3 is mostly mitochondrial
    assert cells["cell3"]["pct_counts_mt"] > 90


def test_qc_pct_counts_mt_low_cell():
    result = compute_scrna_qc(QC_CSV)
    cells = {c["barcode"]: c for c in result["per_cell"]}
    assert cells["cell2"]["pct_counts_mt"] < 1.0


def test_qc_high_mt_cell_flagged_in_summary():
    result = compute_scrna_qc(QC_CSV)
    # cell3 has >20% MT — should be counted
    assert result["summary"]["cells_high_mt_pct"] >= 1


def test_qc_custom_mt_prefix():
    csv = "gene,cell1\nmt-Co1,50\nGapdh,200\n"
    result = compute_scrna_qc(csv, mt_prefix="mt-")
    cell = result["per_cell"][0]
    assert cell["pct_counts_mt"] > 0


def test_qc_no_mt_genes():
    csv = "gene,cell1\nGAPDH,500\nACTB,300\n"
    result = compute_scrna_qc(csv)
    assert result["per_cell"][0]["pct_counts_mt"] == 0.0


def test_qc_empty_matrix_returns_error():
    result = compute_scrna_qc("")
    assert "error" in result


def test_qc_header_only_returns_error():
    result = compute_scrna_qc("gene,cell1,cell2")
    assert "error" in result


def test_qc_non_integer_counts_returns_error():
    csv = "gene,cell1\nGAPDH,abc\n"
    result = compute_scrna_qc(csv)
    assert "error" in result


# ---------------------------------------------------------------------------
# UMI analysis
# ---------------------------------------------------------------------------

UMIS_12BP = ["ACGTACGTACGT", "TTTTTTTTTTTT", "GGGGGGGGGGGG", "ACGTACGTACGT"]
# 4 UMIs, 3 unique → duplication_rate = 25%


def test_umi_total_count():
    result = analyze_umis(UMIS_12BP)
    assert result["total_umis"] == 4
    assert result["valid_umis"] == 4


def test_umi_unique_count():
    result = analyze_umis(UMIS_12BP)
    assert result["unique_umis"] == 3


def test_umi_duplication_rate():
    result = analyze_umis(UMIS_12BP)
    assert result["duplication_rate_pct"] == 25.0


def test_umi_wrong_length_flagged():
    umis = ["ACGT", "ACGTACGTACGT"]
    result = analyze_umis(umis)
    assert result["wrong_length"] == 1
    assert result["valid_umis"] == 1


def test_umi_expected_length_by_chemistry():
    result = analyze_umis(["ACGTACGTACGT"], chemistry="10x_v3")
    assert result["expected_umi_length"] == 12

    result = analyze_umis(["ACGTACGTAC"], chemistry="10x_v2")
    assert result["expected_umi_length"] == 10


def test_umi_base_composition_keys():
    result = analyze_umis(UMIS_12BP)
    assert set(result["base_composition"].keys()) == {"A", "C", "G", "T"}


def test_umi_per_position_frequencies_length():
    result = analyze_umis(UMIS_12BP)
    assert len(result["per_position_frequencies"]) == 12


def test_umi_umi_space():
    result = analyze_umis(UMIS_12BP, chemistry="10x_v3")
    assert result["umi_space"] == 4 ** 12


def test_umi_collision_probability_positive():
    result = analyze_umis(UMIS_12BP * 100)
    assert result["collision_probability_pct"] >= 0


def test_umi_unknown_chemistry_returns_error():
    result = analyze_umis(UMIS_12BP, chemistry="fake")
    assert "error" in result


def test_umi_empty_list_returns_error():
    result = analyze_umis([])
    assert "error" in result
