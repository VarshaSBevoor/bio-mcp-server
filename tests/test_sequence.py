"""Tests for tools/sequence.py — DNA/RNA validation and stats."""

import pytest
from tools.sequence import analyze_sequence


# ---------------------------------------------------------------------------
# Happy path — DNA
# ---------------------------------------------------------------------------

def test_dna_basic_stats():
    result = analyze_sequence("ATGCATGC")
    assert result["valid"] is True
    assert result["type"] == "DNA"
    assert result["length"] == 8
    assert result["base_counts"] == {"A": 2, "C": 2, "G": 2, "T": 2}
    assert result["gc_content_pct"] == 50.0


def test_dna_lowercase_accepted():
    result = analyze_sequence("atgcatgc")
    assert result["valid"] is True
    assert result["length"] == 8


def test_dna_reverse_complement():
    # ATGC → complement TACG → reverse GCAT
    result = analyze_sequence("ATGC")
    assert result["reverse_complement"] == "GCAT"
    assert result["complement"] == "TACG"


def test_dna_gc_content_all_gc():
    result = analyze_sequence("GCGCGC")
    assert result["gc_content_pct"] == 100.0


def test_dna_gc_content_all_at():
    result = analyze_sequence("ATATAT")
    assert result["gc_content_pct"] == 0.0


def test_dna_molecular_weight():
    # 10 bp DNA → 10 * 330 = 3300 Da
    result = analyze_sequence("ATGCATGCAT")
    assert result["molecular_weight_approx_da"] == 3300.0


def test_dna_iupac_ambiguous_accepted():
    # N, R, Y etc. are valid IUPAC DNA
    result = analyze_sequence("ATGCNRYW")
    assert result["valid"] is True
    assert "invalid_characters" not in result


# ---------------------------------------------------------------------------
# Happy path — RNA
# ---------------------------------------------------------------------------

def test_rna_detection():
    # Contains U, no T → RNA
    result = analyze_sequence("AUGCAUGC")
    assert result["type"] == "RNA"
    assert result["valid"] is True


def test_rna_base_counts():
    result = analyze_sequence("AAUGCC")
    counts = result["base_counts"]
    assert counts["A"] == 2
    assert counts["U"] == 1
    assert counts["G"] == 1
    assert counts["C"] == 2


def test_rna_molecular_weight():
    # 8 nt RNA → 8 * 340 = 2720
    result = analyze_sequence("AUGCAUGC")
    assert result["molecular_weight_approx_da"] == 2720.0


def test_rna_reverse_complement_present():
    result = analyze_sequence("AUGC")
    assert "reverse_complement" in result
    assert "complement" not in result  # RNA path omits DNA complement key


# ---------------------------------------------------------------------------
# Error / edge cases
# ---------------------------------------------------------------------------

def test_empty_sequence():
    result = analyze_sequence("")
    assert "error" in result


def test_whitespace_only_sequence():
    result = analyze_sequence("   ")
    assert "error" in result


def test_invalid_characters_reported():
    result = analyze_sequence("ATGCXYZ")
    assert result["valid"] is False
    assert "invalid_characters" in result
    assert "X" in result["invalid_characters"]
    assert "Y" not in result["invalid_characters"]  # Y is valid IUPAC
    assert "Z" in result["invalid_characters"]


def test_invalid_chars_no_stats_returned():
    result = analyze_sequence("ATGCXXX")
    assert "gc_content_pct" not in result
    assert "base_counts" not in result


def test_single_base():
    result = analyze_sequence("A")
    assert result["valid"] is True
    assert result["length"] == 1
