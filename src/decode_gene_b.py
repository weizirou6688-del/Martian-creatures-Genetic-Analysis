"""
Stage 5 - Decode the unknown Gene B (pipeline entry point).

Gene B has no known mRNA or protein. This module learns the Martian
transcription key, codon length, and codon table from the reference
Gene A (stages 1-4), then applies them to transcribe and translate
Gene B into its mRNA and protein sequences.

Run this file to execute the full pipeline:
    python src/decode_gene_b.py
"""

from pathlib import Path

import codon_length
import codon_table
import sequence_analysis
import transcription_decoder

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


# Derive the transcription key and codon lookup table from Gene A.
def generate_basic_information():

    # 1. read the reference sequences
    gene_a = sequence_analysis.read_file_sequence(DATA_DIR / 'gene_a.fa')
    rna_a = sequence_analysis.read_file_sequence(DATA_DIR / 'rna_a.fa')
    protein_a = sequence_analysis.read_file_sequence(DATA_DIR / 'protein_a.fa')

    if not (gene_a and rna_a and protein_a):
        print("Error: Failed to read Gene A reference files")
        return None, None, None

    # 2. get the transcription key
    print("Deriving transcription key")
    offset = transcription_decoder.get_best_alignment_offset(gene_a, rna_a)
    pair_counts = transcription_decoder.get_aligned_sequences(gene_a, rna_a, offset)
    transcription_key, _ = transcription_decoder.derive_transcription_rules(pair_counts)

    # 3. get the codon length
    print("Calculating codon length")
    codon_len = codon_length.get_codon_length(rna_a, protein_a)

    # 4. get the codon table
    print("Generating the codon lookup table")
    table = codon_table.lookup_table_mrna_to_amino(rna_a, protein_a, codon_len)

    return transcription_key, codon_len, table


def transcribe_dna_to_rna(dna_seq, transcription_key):

    rna_seq = []
    for base in dna_seq:
        rna_base = transcription_key.get(base, '?')
        rna_seq.append(rna_base)

    return "".join(rna_seq)


def translate_rna_to_protein(rna_seq, codon_len, table):
    protein_seq = []

    # iterate through the RNA sequence with step = codon_len
    for i in range(0, len(rna_seq), codon_len):
        start = i
        end = i + codon_len

        if end > len(rna_seq):
            break

        codon = rna_seq[start : end]

        amino_acid = table.get(codon, 'X')

        if amino_acid == "STOP":
            print(f"Termination codon detected at index {start}. Translation stopped.")
            break

        protein_seq.append(amino_acid)

    return "".join(protein_seq)


def main():
    print("Stage 5: Gene B analysis")
    key, length, table = generate_basic_information()

    if not (key and length and table):
        print("Error: could not derive transcription and translation data")
        return

    gene_b_seq = sequence_analysis.read_file_sequence(DATA_DIR / 'gene_b.fa')

    if not gene_b_seq:
        print("Error: could not read gene_b.fa")
        return

    print(f"\nGene B length: {len(gene_b_seq)}\n")

    print("Transcribing Gene B to mRNA")
    mrna_b_seq = transcribe_dna_to_rna(gene_b_seq, key)
    print(f"\nGenerated mRNA sequence: {mrna_b_seq}\n")

    print("Translating mRNA to protein\n")
    protein_b_seq = translate_rna_to_protein(mrna_b_seq, length, table)

    print(protein_b_seq)


if __name__ == "__main__":
    main()
