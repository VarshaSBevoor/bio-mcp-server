"""FASTA and FASTQ file/string parsers."""

import io
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


def parse_fasta(content: str) -> dict:
    """Parse a FASTA string and return per-record stats plus a summary."""
    try:
        records = list(SeqIO.parse(io.StringIO(content), "fasta"))
    except ValueError as exc:
        return {"error": str(exc).splitlines()[0]}
    if not records:
        return {"error": "No valid FASTA records found."}

    parsed = []
    lengths = []
    for rec in records:
        seq = str(rec.seq).upper()
        lengths.append(len(seq))
        parsed.append({
            "id": rec.id,
            "description": rec.description,
            "length": len(seq),
            "gc_content_pct": round(gc_fraction(rec.seq) * 100, 2),
        })

    return {
        "record_count": len(records),
        "total_bases": sum(lengths),
        "min_length": min(lengths),
        "max_length": max(lengths),
        "mean_length": round(sum(lengths) / len(lengths), 1),
        "records": parsed,
    }


def parse_fastq(content: str) -> dict:
    """Parse a FASTQ string and return per-record stats plus quality summary."""
    try:
        records = list(SeqIO.parse(io.StringIO(content), "fastq"))
    except ValueError as exc:
        return {"error": str(exc).splitlines()[0]}
    if not records:
        return {"error": "No valid FASTQ records found."}

    parsed = []
    lengths = []
    mean_quals = []

    for rec in records:
        quals = rec.letter_annotations["phred_quality"]
        mean_q = round(sum(quals) / len(quals), 2)
        lengths.append(len(rec.seq))
        mean_quals.append(mean_q)
        parsed.append({
            "id": rec.id,
            "length": len(rec.seq),
            "mean_quality": mean_q,
            "min_quality": min(quals),
            "max_quality": max(quals),
            "gc_content_pct": round(gc_fraction(rec.seq) * 100, 2),
        })

    return {
        "record_count": len(records),
        "total_bases": sum(lengths),
        "mean_read_length": round(sum(lengths) / len(lengths), 1),
        "mean_quality_across_all": round(sum(mean_quals) / len(mean_quals), 2),
        "records": parsed,
    }
