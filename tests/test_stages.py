"""Unit tests for the individual decoding stages."""

import codon_length
import codon_table
import decode_gene_b
import sequence_analysis
import transcription_decoder

# --- Stage 1: sequence cleaning -------------------------------------------


def test_read_file_sequence_cleans_newlines_and_n_padding(tmp_path):
    seq_file = tmp_path / "seq.fa"
    seq_file.write_text("NNNNATCB\nATCB\nATCBNNNN\n")
    assert sequence_analysis.read_file_sequence(seq_file) == "ATCBATCBATCB"


def test_read_file_sequence_missing_file_returns_none(tmp_path):
    assert sequence_analysis.read_file_sequence(tmp_path / "missing.fa") is None


# --- Stage 2: transcription key ------------------------------------------


def test_get_best_alignment_offset_locates_t_u_pairing():
    # the RNA 'U' pairs with the DNA 'T' at index 3
    assert transcription_decoder.get_best_alignment_offset("CCCTACCC", "UA") == 3


def test_get_aligned_sequences_counts_base_pairs():
    pair_counts = transcription_decoder.get_aligned_sequences("ATC", "ZUX", 0)
    assert pair_counts == {("A", "Z"): 1, ("T", "U"): 1, ("C", "X"): 1}


def test_derive_transcription_rules_returns_valid_bijection():
    pair_counts = {("A", "Z"): 10, ("T", "U"): 8, ("C", "X"): 6}
    key, is_valid = transcription_decoder.derive_transcription_rules(pair_counts)
    assert key == {"A": "Z", "C": "X", "T": "U"}
    assert is_valid is True


def test_derive_transcription_rules_rejects_ambiguous_pairing():
    # 'A' maps to 'Z' only 60 % of the time -> not accepted into the key
    pair_counts = {("A", "Z"): 6, ("A", "Y"): 4}
    key, _ = transcription_decoder.derive_transcription_rules(pair_counts)
    assert "A" not in key


# --- Stage 3: codon length -----------------------------------------------


def test_get_codon_length_without_stop_codon():
    # 3 amino acids, mRNA length 6  ->  codon length 2
    assert codon_length.get_codon_length("ABCDEF", "NaMgCl") == 2


def test_get_codon_length_with_stop_codon():
    # 3 amino acids, mRNA length 8  ->  8 / (3 + 1) = codon length 2
    assert codon_length.get_codon_length("ABCDEFGH", "NaMgCl") == 2


def test_get_codon_length_empty_protein_returns_none():
    assert codon_length.get_codon_length("ABCDEF", "") is None


# --- Stage 4: codon table ------------------------------------------------


def test_lookup_table_maps_codons_to_amino_acids():
    table = codon_table.lookup_table_mrna_to_amino("UUZX", "AsBe", 2)
    assert table == {"UU": "As", "ZX": "Be"}


def test_lookup_table_records_trailing_stop_codon():
    # the mRNA carries one extra codon ('XX') beyond the protein
    table = codon_table.lookup_table_mrna_to_amino("UUZXXX", "AsBe", 2)
    assert table["XX"] == "STOP"


# --- Stage 5: transcription & translation --------------------------------


def test_transcribe_dna_to_rna_applies_key():
    key = {"A": "Z", "B": "Y", "C": "X", "T": "U"}
    assert decode_gene_b.transcribe_dna_to_rna("ATCB", key) == "ZUXY"


def test_transcribe_marks_unknown_base():
    assert decode_gene_b.transcribe_dna_to_rna("AG", {"A": "Z"}) == "Z?"


def test_translate_rna_to_protein_builds_sequence():
    table = {"UU": "As", "ZX": "Be"}
    assert decode_gene_b.translate_rna_to_protein("UUZX", 2, table) == "AsBe"


def test_translate_halts_at_stop_codon():
    table = {"UU": "As", "XX": "STOP", "ZZ": "Be"}
    assert decode_gene_b.translate_rna_to_protein("UUXXZZ", 2, table) == "As"


def test_translate_marks_unknown_codon():
    assert decode_gene_b.translate_rna_to_protein("UUZZ", 2, {"UU": "As"}) == "AsX"
