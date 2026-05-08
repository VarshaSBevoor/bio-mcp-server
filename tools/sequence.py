"""DNA/RNA sequence validation and statistics."""

from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

IUPAC_DNA = set("ACGTRYSWKMBDHVNacgtryswkmbdhvn")
IUPAC_RNA = set("ACGURYSWKMBDHVNacguryswkmbdhvn")

COMPLEMENT_TABLE = str.maketrans("ACGTacgt", "TGCAtgca")


def analyze_sequence(sequence: str) -> dict:
    """Return stats and validation for a DNA or RNA sequence string."""
    seq = sequence.strip().upper()

    if not seq:
        return {"error": "Empty sequence provided."}

    is_rna = "U" in seq and "T" not in seq
    alphabet = IUPAC_RNA if is_rna else IUPAC_DNA
    invalid = sorted({c for c in seq if c not in alphabet})

    result: dict = {
        "length": len(seq),
        "type": "RNA" if is_rna else "DNA",
        "valid": len(invalid) == 0,
    }

    if invalid:
        result["invalid_characters"] = invalid
        return result

    bio_seq = Seq(seq)
    counts = {base: seq.count(base) for base in ("A", "C", "G", "T" if not is_rna else "U")}
    result["base_counts"] = counts
    result["gc_content_pct"] = round(gc_fraction(bio_seq) * 100, 2)

    if not is_rna:
        result["reverse_complement"] = str(bio_seq.reverse_complement())
        result["complement"] = str(bio_seq.complement())
    else:
        result["reverse_complement"] = str(bio_seq.reverse_complement_rna())

    result["molecular_weight_approx_da"] = round(len(seq) * (330 if not is_rna else 340), 1)

    return result
