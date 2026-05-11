"""Bioinformatics MCP server — exposes DNA/RNA analysis, BLAST, sequence parsing, and scRNA-seq tools."""

from mcp.server.fastmcp import FastMCP
from tools.sequence import analyze_sequence
from tools.parsers import parse_fasta, parse_fastq
from tools.blast import run_blast
from tools.scrna import validate_barcodes, compute_scrna_qc, analyze_umis

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


@mcp.tool()
def validate_cell_barcodes(
    barcodes: list,
    chemistry: str = "10x_v3",
    whitelist: str | None = None,
) -> dict:
    """
    Validate cell barcodes from a single-cell sequencing experiment.

    Checks barcode length against the expected value for the given chemistry,
    flags barcodes containing N bases, and — when a whitelist is provided —
    reports exact matches and barcodes correctable within Hamming distance 1
    (matching Cell Ranger's correction strategy).

    Args:
        barcodes: List of barcode strings extracted from R1 reads.
        chemistry: Library chemistry — 10x_v2, 10x_v3, 10x_v3.1, dropseq, indrop.
        whitelist: Optional newline-separated string of known valid barcodes.
    """
    return validate_barcodes(barcodes, chemistry=chemistry, whitelist=whitelist)


@mcp.tool()
def scrna_qc_metrics(count_matrix_csv: str, mt_prefix: str = "MT-") -> dict:
    """
    Compute per-cell QC metrics from a raw count matrix.

    Calculates the three core metrics used in scRNA-seq QC pipelines
    (e.g. Scanpy's calculate_qc_metrics):
    - total_counts: total UMI counts per cell
    - n_genes_by_counts: genes detected per cell (count > 0)
    - pct_counts_mt: fraction of counts from mitochondrial genes

    High pct_counts_mt (>20%) flags dying cells. Outlier total_counts or
    n_genes_by_counts flags empty droplets or doublets.

    Args:
        count_matrix_csv: Genes-x-cells or cells-x-genes CSV with row/column labels.
        mt_prefix: Gene name prefix for mitochondrial genes (default: MT-).
    """
    return compute_scrna_qc(count_matrix_csv, mt_prefix=mt_prefix)


@mcp.tool()
def umi_analysis(umis: list, chemistry: str = "10x_v3") -> dict:
    """
    Analyze UMI sequences for quality, nucleotide bias, and collision risk.

    UMIs allow PCR duplicates to be collapsed so each unique UMI+gene combination
    represents one original molecule. This tool reports UMI length consistency,
    per-position nucleotide frequencies (to detect synthesis bias), duplication rate,
    and estimates collision probability via the birthday problem approximation.

    Args:
        umis: List of UMI strings (for 10x v3: bases 17-28 of the R1 read).
        chemistry: Library chemistry — 10x_v2 (10 bp), 10x_v3/v3.1 (12 bp),
                   dropseq (8 bp), indrop (6 bp).
    """
    return analyze_umis(umis, chemistry=chemistry)


if __name__ == "__main__":
    mcp.run()
