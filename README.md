# bio-mcp-server

A [Model Context Protocol (MCP)](https://modelcontextprotocol.io) server that exposes bioinformatics tools as callable functions for AI assistants like Claude. Built with [FastMCP](https://github.com/modelcontextprotocol/python-sdk) and [Biopython](https://biopython.org).

---

## What is MCP?

MCP is an open protocol that lets AI assistants call external tools during a conversation — similar to how a browser extension adds capabilities to a web app. This server makes bioinformatics functions available to Claude (or any MCP-compatible client) so it can validate sequences, run BLAST searches, and compute scRNA-seq QC metrics on demand.

---

## Tools

### Sequence analysis

| Tool | Description |
|---|---|
| `analyze_dna_rna` | Validates IUPAC characters, computes GC content, base counts, complement/reverse-complement, and approximate molecular weight. Auto-detects DNA vs RNA. |
| `parse_fasta_sequences` | Parses multi-record FASTA text. Returns per-record length and GC content plus aggregate stats (total bases, min/max/mean length). |
| `parse_fastq_reads` | Parses FASTQ text. Returns per-read Phred quality (mean/min/max) and GC content, plus an overall quality summary. |
| `blast_sequence` | Submits a sequence to NCBI BLAST and returns top hits with E-value, percent identity, and query coverage. Supports all five BLAST programs and multiple databases. |

### Single-cell RNA-seq

| Tool | Description |
|---|---|
| `validate_cell_barcodes` | Validates cell barcodes against expected length for a given chemistry (10x v2/v3/v3.1, Drop-seq, inDrop). When a whitelist is provided, reports exact matches and barcodes correctable within Hamming distance 1 — the same correction strategy used by Cell Ranger. |
| `scrna_qc_metrics` | Computes the three core scRNA-seq QC metrics per cell from a raw count matrix CSV: `total_counts`, `n_genes_by_counts`, and `pct_counts_mt`. Mirrors `scanpy.pp.calculate_qc_metrics`. Flags cells with >20% mitochondrial counts. |
| `umi_analysis` | Analyzes UMI sequences for length consistency, per-position nucleotide bias, duplication rate, and PCR collision probability (birthday problem approximation). Supports 10x v2/v3, Drop-seq, and inDrop UMI lengths. |

---

## Setup

**Requirements:** Python 3.12+

```bash
git clone git@github.com:VarshaSBevoor/bio-mcp-server.git
cd bio-mcp-server
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

### Run the server

```bash
python server.py
```

### Connect to Claude Desktop

Add the following to your `claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "bio": {
      "command": "/path/to/bio-mcp-server/venv/bin/python",
      "args": ["/path/to/bio-mcp-server/server.py"]
    }
  }
}
```

On macOS the config file is at `~/Library/Application Support/Claude/claude_desktop_config.json`.

---

## Tests

```bash
source venv/bin/activate
pytest tests/ -v
```

86 tests covering happy paths, edge cases, and mocked BLAST network calls.

---

## Project structure

```
bio-mcp-server/
├── server.py          # FastMCP server — registers all tools
├── requirements.txt
└── tools/
    ├── sequence.py    # DNA/RNA validation and statistics
    ├── parsers.py     # FASTA and FASTQ parsing
    ├── blast.py       # NCBI BLAST via Biopython NCBIWWW
    └── scrna.py       # Single-cell barcode, QC, and UMI tools
```

---

## Design notes

**Why return dicts instead of raising exceptions?**  
MCP tools serialize their return values to JSON for the AI client. Returning a structured `{"error": "..."}` keeps the interface consistent and lets the model reason about what went wrong rather than receiving an opaque stack trace.

**BLAST tests use mocks**  
`tests/test_blast.py` patches `NCBIWWW.qblast` so the suite runs fast and deterministically without hitting NCBI. Validation logic (program names, database names, sequence length) is tested directly; network behavior is tested via `unittest.mock`.

**Barcode correction**  
The Hamming-1 correction in `validate_barcodes` iterates over all single-base substitutions for each input barcode and checks membership in a whitelist set (O(1) lookup). This is intentionally simple — production Cell Ranger uses a precomputed correction map, but the logic is equivalent for the barcode lengths in use.
