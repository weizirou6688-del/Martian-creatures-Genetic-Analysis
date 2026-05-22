"""
Stage 3 - Determine the Martian codon length.

The codon length is derived mathematically from the mRNA and protein
lengths. Two hypotheses are tested - a sequence with a stop codon and a
sequence without one - and whichever yields an integer codon length is
accepted as the correct biological rule.
"""

import re
from pathlib import Path

import sequence_analysis

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def get_codon_length(rna_seq: str, protein_seq: str) -> int | None:
    """Derive the codon length from the mRNA and protein lengths.

    Codon length = len(mRNA) / number_of_amino_acids. If that is not an
    integer, the calculation is retried assuming one extra (stop) codon.

    Args:
        rna_seq: The mRNA sequence.
        protein_seq: The protein sequence (element symbols, e.g. 'NaMg').

    Returns:
        The integer codon length, or ``None`` if none can be found.
    """
    rna_len = len(rna_seq)
    amino_acids = re.findall(r"[A-Z][a-z]?", protein_seq)
    protein_len = len(amino_acids)

    if protein_len == 0:
        print("Error: Protein sequence is empty")
        return None

    print(f"RNA length: {rna_len}, Protein length {protein_len}")

    calculated_codon_len_a = rna_len / protein_len
    calculated_codon_len_b = rna_len / (protein_len + 1)

    if calculated_codon_len_a % 1 == 0:
        codon_len = int(calculated_codon_len_a)
        print("Result: No termination codon")
        print(f"Codon length is {codon_len}")
        return codon_len
    elif calculated_codon_len_b % 1 == 0:
        codon_len = int(calculated_codon_len_b)
        print("Result: Have termination codon")
        print(f"Codon length is {codon_len}")
        return codon_len
    else:
        print("Error: Could not find integer codon length")
        return None


def main() -> None:
    """Run stage 3 standalone on the Gene A reference data."""
    print("Stage 3: determine codon length")
    rna_seq = sequence_analysis.read_file_sequence(DATA_DIR / "rna_a.fa")
    protein_seq = sequence_analysis.read_file_sequence(DATA_DIR / "protein_a.fa")

    if rna_seq and protein_seq:
        codon_len = get_codon_length(rna_seq, protein_seq)
        if codon_len:
            print(f"\nThe Martian codon length is {codon_len}")
    else:
        print("Error: Failed to read sequence files")


if __name__ == "__main__":
    main()
