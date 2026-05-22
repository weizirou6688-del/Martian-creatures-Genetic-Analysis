"""End-to-end test of the full Gene B decoding pipeline on the real data.

This verifies that the pipeline reproduces the documented scientific
result: the Martian genetic rules learned from Gene A, and the protein
sequence decoded for the unknown Gene B.
"""

import re

import decode_gene_b
import sequence_analysis

# The protein produced by decoding Gene B - the documented result of the
# project (see docs/bioinformatic_report.pdf).
EXPECTED_GENE_B_PROTEIN = "AsArAlKNaClNBeBeRnClNaArKNiPbBeSRnNaArPbMgLiNiLiArClAlMgRnArCNClPbNKAlNLiNaClSNaSRnSNaKBeLiPbKMgNCArMgKRnAlNaAlNNClCKNiRnNSAlLiClMgAlPbNiNAlMgNaKCCPbNaArSAlArNNaCCKArNaClZn"


def test_rules_learned_from_gene_a():
    key, codon_len, table = decode_gene_b.generate_basic_information()
    assert key == {"A": "Z", "B": "Y", "C": "X", "T": "U"}
    assert codon_len == 2
    assert len(table) == 16  # 16 codons, one per Martian amino acid


def test_gene_b_decodes_to_expected_protein():
    key, codon_len, table = decode_gene_b.generate_basic_information()

    gene_b = sequence_analysis.read_file_sequence(
        sequence_analysis.DATA_DIR / "gene_b.fa"
    )
    mrna = decode_gene_b.transcribe_dna_to_rna(gene_b, key)
    protein = decode_gene_b.translate_rna_to_protein(mrna, codon_len, table)

    assert protein == EXPECTED_GENE_B_PROTEIN
    # 204-base gene -> 204-base mRNA -> 102 doublet codons
    assert len(re.findall(r"[A-Z][a-z]?", protein)) == 102
