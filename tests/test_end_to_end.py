"""Integration tests: the full decoding pipeline and module entry points.

`test_gene_b_*` verify that the pipeline reproduces the documented
scientific result. The smoke tests confirm that every module's ``main()``
runs end to end on the real data without raising.
"""

import re

import codon_length
import codon_table
import decode_gene_b
import sequence_analysis
import transcription_decoder

# The protein produced by decoding Gene B - the documented result of the
# project (see docs/bioinformatic_report.pdf).
EXPECTED_GENE_B_PROTEIN = "AsArAlKNaClNBeBeRnClNaArKNiPbBeSRnNaArPbMgLiNiLiArClAlMgRnArCNClPbNKAlNLiNaClSNaSRnSNaKBeLiPbKMgNCArMgKRnAlNaAlNNClCKNiRnNSAlLiClMgAlPbNiNAlMgNaKCCPbNaArSAlArNNaCCKArNaClZn"  # noqa: E501

# --- Full pipeline -------------------------------------------------------


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


# --- Smoke tests: every module's main() runs without error ---------------


def test_sequence_analysis_main_runs(tmp_path, monkeypatch):
    monkeypatch.setattr(sequence_analysis, "RESULTS_DIR", tmp_path)
    sequence_analysis.main()


def test_transcription_decoder_main_runs():
    transcription_decoder.main()


def test_codon_length_main_runs():
    codon_length.main()


def test_codon_table_main_runs():
    codon_table.main()


def test_decode_gene_b_main_runs():
    decode_gene_b.main()
