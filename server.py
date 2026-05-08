"""Bioinformatics MCP server — exposes DNA/RNA analysis, BLAST, and sequence file parsing."""

from mcp.server.fastmcp import FastMCP
from tools.sequence import analyze_sequence
from tools.parsers import parse_fasta, parse_fastq
from tools.blast import run_blast

mcp = FastMCP("bio-mcp-server")


@mcp.tool()
def analyze_dna_rna(sequence: str) -> dict:
    """
    Validate and compute statistics for a DNA or RNA sequence.

    Returns length, base counts, GC content, complement/reverse-complement,
    and an approximate molecular weight. Automatically detects DNA vs RNA.
    Invalid IUPAC characters are reported rather than raising an error.

    Args:
        sequence: Raw nucleotide sequence string (DNA or RNA, any case).
    """
    return analyze_sequence(sequence)


@mcp.tool()
def parse_fasta_sequences(fasta_content: str) -> dict:
    """
    Parse a multi-record FASTA string and return per-record stats and a summary.

    Stats include sequence ID, description, length, and GC content for each record,
    plus aggregate totals and length distribution.

    Args:
        fasta_content: Full text of a FASTA file (one or more records).
    """
    return parse_fasta(fasta_content)


@mcp.tool()
def parse_fastq_reads(fastq_content: str) -> dict:
    """
    Parse a FASTQ string and return per-read quality and GC stats plus a summary.

    Per-read output includes mean/min/max Phred quality score and GC content.
    Summary includes total reads, total bases, mean read length, and overall quality.

    Args:
        fastq_content: Full text of a FASTQ file (one or more reads).
    """
    return parse_fastq(fastq_content)


@mcp.tool()
def blast_sequence(
    sequence: str,
    program: str = "blastn",
    database: str = "nt",
    max_hits: int = 5,
) -> dict:
    """
    Submit a sequence to NCBI BLAST and return the top hits.

    Makes a live network request to NCBI — may take 30–60 seconds.
    Results include accession, title, E-value, percent identity, and query coverage
    for each hit.

    Args:
        sequence: Raw nucleotide or protein sequence.
        program: BLAST program to use — blastn, blastp, blastx, tblastn, tblastx.
        database: Database to search — nt, nr, refseq_rna, refseq_protein, swissprot.
        max_hits: Number of top hits to return (1–20, default 5).
    """
    return run_blast(sequence, program=program, database=database, max_hits=max_hits)


if __name__ == "__main__":
    mcp.run()
