"""
Stage 4 - Build the codon-to-amino-acid lookup table.

The mRNA sequence is split into fixed-length codons and each codon is
mapped to the amino acid at the matching position in the reference
protein. If the mRNA is longer than the protein implies, the trailing
codon is recorded as the stop codon.
"""

import re
from collections import Counter
from pathlib import Path

import codon_length
import sequence_analysis

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


# Determine the mRNA -> amino acid table.
# Each codon corresponds to exactly one amino acid.
def lookup_table_mrna_to_amino(rna_seq, protein_seq, codon_len):
    # 1. prepare an empty mRNA -> amino acid table
    codon_table = {}

    # the codon length may arrive as a float from a division; normalize it
    codon_len = int(codon_len)

    # findall returns a list of amino acid symbols (e.g. 'Na', 'Mg')
    amino_acids = re.findall(r'[A-Z][a-z]?', protein_seq)
    counts = Counter(amino_acids)
    print(f"The number of amino acid types : {len(counts)}")

    # record every codon -> amino acid pairing
    for i in range(len(amino_acids)):
        current_aa = amino_acids[i]

        # the i-th amino acid is encoded by the i-th codon_len-sized
        # window of the mRNA, so the indices line up directly
        start_index = i * codon_len
        end_index = start_index + codon_len

        if end_index <= len(rna_seq):
            codon = rna_seq[start_index : end_index]

            if codon not in codon_table:
                codon_table[codon] = current_aa

    # check for a stop codon: an mRNA carries exactly one. If the mRNA is
    # longer than the protein accounts for, the trailing codon is the stop.
    expected_rna_len = len(amino_acids) * codon_len
    if len(rna_seq) > expected_rna_len:
        stop_codon = rna_seq[expected_rna_len : expected_rna_len + codon_len]
        codon_table[stop_codon] = "STOP"

    return codon_table


def main():
    print("Stage 4: Create amino acid and codon lookup table")

    rna_seq = sequence_analysis.read_file_sequence(DATA_DIR / 'rna_a.fa')
    protein_seq = sequence_analysis.read_file_sequence(DATA_DIR / 'protein_a.fa')

    if rna_seq and protein_seq:
        codon_len = codon_length.get_codon_length(rna_seq, protein_seq)
        if codon_len:
            print(f"\nThe Martian codon length is {codon_len}")
            table = lookup_table_mrna_to_amino(rna_seq, protein_seq, codon_len)
            print("codon-amino acid query table")

            for codon, amino in table.items():
                print(f"{codon} : {amino}")
        else:
            print("Error: Failed to read codon length")
    else:
        print("Error: Failed to read sequence files")


if __name__ == "__main__":
    main()
