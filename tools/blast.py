"""NCBI BLAST queries via Biopython's NCBIWWW interface."""

import time
from Bio.Blast import NCBIWWW, NCBIXML

VALID_PROGRAMS = {"blastn", "blastp", "blastx", "tblastn", "tblastx"}
VALID_DATABASES = {"nt", "nr", "refseq_rna", "refseq_protein", "swissprot"}


def run_blast(
    sequence: str,
    program: str = "blastn",
    database: str = "nt",
    max_hits: int = 5,
) -> dict:
    """
    Submit a sequence to NCBI BLAST and return the top hits.

    Args:
        sequence: Raw nucleotide or protein sequence string.
        program: BLAST program — blastn, blastp, blastx, tblastn, tblastx.
        database: Target database — nt, nr, refseq_rna, refseq_protein, swissprot.
        max_hits: Maximum number of hits to return (1–20).
    """
    program = program.lower()
    database = database.lower()

    if program not in VALID_PROGRAMS:
        return {"error": f"Invalid program '{program}'. Choose from: {sorted(VALID_PROGRAMS)}"}
    if database not in VALID_DATABASES:
        return {"error": f"Invalid database '{database}'. Choose from: {sorted(VALID_DATABASES)}"}

    max_hits = max(1, min(20, max_hits))
    sequence = sequence.strip()

    if len(sequence) < 10:
        return {"error": "Sequence too short for BLAST (minimum 10 characters)."}

    try:
        result_handle = NCBIWWW.qblast(
            program,
            database,
            sequence,
            hitlist_size=max_hits,
        )
    except Exception as exc:
        return {"error": f"BLAST submission failed: {exc}"}

    try:
        blast_record = NCBIXML.read(result_handle)
    except Exception as exc:
        return {"error": f"Failed to parse BLAST results: {exc}"}

    hits = []
    for alignment in blast_record.alignments[:max_hits]:
        hsp = alignment.hsps[0]
        hits.append({
            "title": alignment.title[:120],
            "accession": alignment.accession,
            "length": alignment.length,
            "score": hsp.score,
            "e_value": hsp.expect,
            "identity_pct": round(hsp.identities / hsp.align_length * 100, 1),
            "query_cover_pct": round(
                (hsp.query_end - hsp.query_start) / len(sequence) * 100, 1
            ),
        })

    return {
        "program": program,
        "database": database,
        "query_length": len(sequence),
        "hits_returned": len(hits),
        "hits": hits,
    }
