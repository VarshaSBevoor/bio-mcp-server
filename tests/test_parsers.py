"""Tests for tools/parsers.py — FASTA and FASTQ parsing."""

import pytest
from tools.parsers import parse_fasta, parse_fastq

# ---------------------------------------------------------------------------
# FASTA
# ---------------------------------------------------------------------------

SINGLE_FASTA = ">seq1 human gene\nATGCATGCATGC\n"

MULTI_FASTA = (
    ">seq1\nAAAA\n"
    ">seq2\nGGGGCCCC\n"
    ">seq3\nATATATATATAT\n"
)

def test_fasta_single_record():
    result = parse_fasta(SINGLE_FASTA)
    assert result["record_count"] == 1
    assert result["total_bases"] == 12
    assert result["records"][0]["id"] == "seq1"
    assert result["records"][0]["length"] == 12


def test_fasta_gc_content():
    result = parse_fasta(SINGLE_FASTA)
    assert result["records"][0]["gc_content_pct"] == 50.0


def test_fasta_multi_record_count():
    result = parse_fasta(MULTI_FASTA)
    assert result["record_count"] == 3


def test_fasta_multi_total_bases():
    # 4 + 8 + 12 = 24
    result = parse_fasta(MULTI_FASTA)
    assert result["total_bases"] == 24


def test_fasta_min_max_length():
    result = parse_fasta(MULTI_FASTA)
    assert result["min_length"] == 4
    assert result["max_length"] == 12


def test_fasta_mean_length():
    result = parse_fasta(MULTI_FASTA)
    assert result["mean_length"] == 8.0


def test_fasta_description_preserved():
    result = parse_fasta(SINGLE_FASTA)
    assert "human gene" in result["records"][0]["description"]


def test_fasta_all_gc():
    result = parse_fasta(">gc\nGGGCCC\n")
    assert result["records"][0]["gc_content_pct"] == 100.0


def test_fasta_all_at():
    result = parse_fasta(">at\nAAAATTTT\n")
    assert result["records"][0]["gc_content_pct"] == 0.0


def test_fasta_empty_input():
    result = parse_fasta("")
    assert "error" in result


def test_fasta_invalid_input():
    result = parse_fasta("not a fasta file at all")
    assert "error" in result


# ---------------------------------------------------------------------------
# FASTQ
# ---------------------------------------------------------------------------

# 'I' = Phred 40, 'A' = Phred 32
SINGLE_FASTQ = "@read1\nATGCATGC\n+\nIIIIIIII\n"

MULTI_FASTQ = (
    "@read1\nATGC\n+\nIIII\n"     # mean Q=40, GC=50%
    "@read2\nGGGG\n+\nAAAA\n"     # mean Q=32, GC=100%
)


def test_fastq_single_record():
    result = parse_fastq(SINGLE_FASTQ)
    assert result["record_count"] == 1
    assert result["total_bases"] == 8


def test_fastq_mean_quality():
    result = parse_fastq(SINGLE_FASTQ)
    assert result["records"][0]["mean_quality"] == 40.0


def test_fastq_gc_content():
    result = parse_fastq(SINGLE_FASTQ)
    assert result["records"][0]["gc_content_pct"] == 50.0


def test_fastq_multi_record_count():
    result = parse_fastq(MULTI_FASTQ)
    assert result["record_count"] == 2


def test_fastq_total_bases():
    result = parse_fastq(MULTI_FASTQ)
    assert result["total_bases"] == 8


def test_fastq_mean_read_length():
    result = parse_fastq(MULTI_FASTQ)
    assert result["mean_read_length"] == 4.0


def test_fastq_overall_mean_quality():
    # read1 mean=40, read2 mean=32 → average=36
    result = parse_fastq(MULTI_FASTQ)
    assert result["mean_quality_across_all"] == 36.0


def test_fastq_per_read_min_max_quality():
    # All 'I' (Phred 40) → min==max==mean==40
    result = parse_fastq(SINGLE_FASTQ)
    rec = result["records"][0]
    assert rec["min_quality"] == 40
    assert rec["max_quality"] == 40


def test_fastq_high_gc_read():
    result = parse_fastq(MULTI_FASTQ)
    read2 = result["records"][1]
    assert read2["gc_content_pct"] == 100.0


def test_fastq_empty_input():
    result = parse_fastq("")
    assert "error" in result


def test_fastq_invalid_input():
    result = parse_fastq("this is not fastq")
    assert "error" in result
