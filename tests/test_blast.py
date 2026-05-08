"""Tests for tools/blast.py — BLAST validation and mocked network calls."""

import io
import pytest
from unittest.mock import patch, MagicMock
from tools.blast import run_blast

VALID_SEQ = "ATGCATGCATGCATGCATGC"  # 20 nt, passes length check

# ---------------------------------------------------------------------------
# Input validation (no network calls)
# ---------------------------------------------------------------------------

def test_invalid_program():
    result = run_blast(VALID_SEQ, program="badprog")
    assert "error" in result
    assert "badprog" in result["error"]


def test_invalid_database():
    result = run_blast(VALID_SEQ, database="unicorn_db")
    assert "error" in result
    assert "unicorn_db" in result["error"]


def test_sequence_too_short():
    result = run_blast("ATGC")  # < 10 chars
    assert "error" in result
    assert "too short" in result["error"].lower()


def test_sequence_exactly_at_limit():
    # 10 chars should pass validation (error will come from BLAST, not validation)
    # We mock the network so this just confirms no early-exit error
    mock_record = _make_mock_blast_record([])
    with patch("tools.blast.NCBIWWW.qblast") as mock_qblast, \
         patch("tools.blast.NCBIXML.read", return_value=mock_record):
        mock_qblast.return_value = io.StringIO("")
        result = run_blast("ATGCATGCAT")  # exactly 10
    assert "error" not in result


def test_max_hits_clamped_low():
    mock_record = _make_mock_blast_record([])
    with patch("tools.blast.NCBIWWW.qblast") as mock_qblast, \
         patch("tools.blast.NCBIXML.read", return_value=mock_record):
        mock_qblast.return_value = io.StringIO("")
        run_blast(VALID_SEQ, max_hits=0)
        # qblast should have been called with hitlist_size=1 (clamped from 0)
        _, kwargs = mock_qblast.call_args
        assert kwargs["hitlist_size"] == 1


def test_max_hits_clamped_high():
    mock_record = _make_mock_blast_record([])
    with patch("tools.blast.NCBIWWW.qblast") as mock_qblast, \
         patch("tools.blast.NCBIXML.read", return_value=mock_record):
        mock_qblast.return_value = io.StringIO("")
        run_blast(VALID_SEQ, max_hits=999)
        _, kwargs = mock_qblast.call_args
        assert kwargs["hitlist_size"] == 20


def test_program_case_insensitive():
    mock_record = _make_mock_blast_record([])
    with patch("tools.blast.NCBIWWW.qblast") as mock_qblast, \
         patch("tools.blast.NCBIXML.read", return_value=mock_record):
        mock_qblast.return_value = io.StringIO("")
        result = run_blast(VALID_SEQ, program="BLASTN")
    assert "error" not in result


# ---------------------------------------------------------------------------
# Successful BLAST response (mocked)
# ---------------------------------------------------------------------------

def test_successful_blast_returns_hits():
    hit = _make_hit(
        title="Homo sapiens chromosome 1",
        accession="NC_000001",
        length=248956422,
        score=95.0,
        expect=1e-20,
        identities=19,
        align_length=20,
        query_start=1,
        query_end=20,
    )
    mock_record = _make_mock_blast_record([hit])

    with patch("tools.blast.NCBIWWW.qblast") as mock_qblast, \
         patch("tools.blast.NCBIXML.read", return_value=mock_record):
        mock_qblast.return_value = io.StringIO("")
        result = run_blast(VALID_SEQ)

    assert result["hits_returned"] == 1
    assert result["hits"][0]["accession"] == "NC_000001"
    assert result["hits"][0]["e_value"] == 1e-20
    assert result["hits"][0]["identity_pct"] == 95.0
    assert result["hits"][0]["score"] == 95.0


def test_successful_blast_metadata():
    mock_record = _make_mock_blast_record([])
    with patch("tools.blast.NCBIWWW.qblast") as mock_qblast, \
         patch("tools.blast.NCBIXML.read", return_value=mock_record):
        mock_qblast.return_value = io.StringIO("")
        result = run_blast(VALID_SEQ, program="blastn", database="nt")

    assert result["program"] == "blastn"
    assert result["database"] == "nt"
    assert result["query_length"] == len(VALID_SEQ)


def test_no_hits_returns_empty_list():
    mock_record = _make_mock_blast_record([])
    with patch("tools.blast.NCBIWWW.qblast") as mock_qblast, \
         patch("tools.blast.NCBIXML.read", return_value=mock_record):
        mock_qblast.return_value = io.StringIO("")
        result = run_blast(VALID_SEQ)

    assert result["hits"] == []
    assert result["hits_returned"] == 0


def test_blast_network_error_returns_error_key():
    with patch("tools.blast.NCBIWWW.qblast", side_effect=Exception("timeout")):
        result = run_blast(VALID_SEQ)
    assert "error" in result
    assert "timeout" in result["error"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_hit(title, accession, length, score, expect,
              identities, align_length, query_start, query_end):
    hsp = MagicMock()
    hsp.score = score
    hsp.expect = expect
    hsp.identities = identities
    hsp.align_length = align_length
    hsp.query_start = query_start
    hsp.query_end = query_end

    alignment = MagicMock()
    alignment.title = title
    alignment.accession = accession
    alignment.length = length
    alignment.hsps = [hsp]
    return alignment


def _make_mock_blast_record(alignments: list):
    record = MagicMock()
    record.alignments = alignments
    return record
