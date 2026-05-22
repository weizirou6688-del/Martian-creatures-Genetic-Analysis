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


def generate_basic_information() -> (
    tuple[dict[str, str] | None, int | None, dict[str, str] | None]
):
    """Learn the Martian genetic rules from the Gene A reference data.

    Returns:
        A tuple ``(transcription_key, codon_length, codon_table)``, or
        ``(None, None, None)`` if the reference data cannot be decoded.
    """
    gene_a = sequence_analysis.read_file_sequence(DATA_DIR / "gene_a.fa")
    rna_a = sequence_analysis.read_file_sequence(DATA_DIR / "rna_a.fa")
    protein_a = sequence_analysis.read_file_sequence(DATA_DIR / "protein_a.fa")

    if not (gene_a and rna_a and protein_a):
        print("Error: Failed to read Gene A reference files")
        return None, None, None

    print("Deriving transcription key")
    offset = transcription_decoder.get_best_alignment_offset(gene_a, rna_a)
    pair_counts = transcription_decoder.get_aligned_sequences(gene_a, rna_a, offset)
    transcription_key, _ = transcription_decoder.derive_transcription_rules(pair_counts)

    print("Calculating codon length")
    codon_len = codon_length.get_codon_length(rna_a, protein_a)
    if codon_len is None:
        print("Error: could not determine codon length")
        return None, None, None

    print("Generating the codon lookup table")
    table = codon_table.lookup_table_mrna_to_amino(rna_a, protein_a, codon_len)

    return transcription_key, codon_len, table


def transcribe_dna_to_rna(dna_seq: str, transcription_key: dict[str, str]) -> str:
    """Transcribe a DNA sequence to mRNA using the transcription key.

    Bases absent from the key are marked ``'?'``.

    Args:
        dna_seq: The DNA sequence to transcribe.
        transcription_key: DNA -> RNA base mapping.

    Returns:
        The transcribed mRNA sequence.
    """
    rna_seq = []
    for base in dna_seq:
        rna_seq.append(transcription_key.get(base, "?"))
    return "".join(rna_seq)


def translate_rna_to_protein(
    rna_seq: str, codon_len: int, table: dict[str, str]
) -> str:
    """Translate an mRNA sequence into a protein sequence.

    Translation stops at the first stop codon. Codons absent from the
    table are marked ``'X'``.

    Args:
        rna_seq: The mRNA sequence to translate.
        codon_len: Number of nucleotides per codon.
        table: Codon -> amino-acid lookup table.

    Returns:
        The translated protein sequence.
    """
    protein_seq = []
    for i in range(0, len(rna_seq), codon_len):
        codon = rna_seq[i : i + codon_len]
        if len(codon) < codon_len:
            break

        amino_acid = table.get(codon, "X")
        if amino_acid == "STOP":
            print(f"Termination codon detected at index {i}. Translation stopped.")
            break

        protein_seq.append(amino_acid)
    return "".join(protein_seq)


def main() -> None:
    """Run the full Gene B decoding pipeline."""
    print("Stage 5: Gene B analysis")
    key, length, table = generate_basic_information()

    if not (key and length and table):
        print("Error: could not derive transcription and translation data")
        return

    gene_b_seq = sequence_analysis.read_file_sequence(DATA_DIR / "gene_b.fa")
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


if __name__ == "__main__":  # pragma: no cover
    main()
