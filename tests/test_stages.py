"""Unit tests for the individual decoding stages."""

import codon_length
import codon_table
import decode_gene_b
import sequence_analysis
import transcription_decoder

# --- Stage 1: sequence cleaning ------------------------------------------


def test_read_file_sequence_cleans_newlines_and_n_padding(tmp_path):
    seq_file = tmp_path / "seq.fa"
    seq_file.write_text("NNNNATCB\nATCB\nATCBNNNN\n")
    assert sequence_analysis.read_file_sequence(seq_file) == "ATCBATCBATCB"


def test_read_file_sequence_missing_file_returns_none(tmp_path):
    assert sequence_analysis.read_file_sequence(tmp_path / "missing.fa") is None


def test_read_file_sequence_all_padding_returns_empty(tmp_path):
    seq_file = tmp_path / "all_n.fa"
    seq_file.write_text("NNNNNNNN\n")
    assert sequence_analysis.read_file_sequence(seq_file) == ""


def test_read_file_sequence_empty_file_returns_empty(tmp_path):
    seq_file = tmp_path / "empty.fa"
    seq_file.write_text("")
    assert sequence_analysis.read_file_sequence(seq_file) == ""


def test_analyze_sequence_warns_on_empty_input(capsys):
    sequence_analysis.analyze_sequence("", "Empty Gene")
    assert "empty" in capsys.readouterr().out.lower()


def test_visualize_counts_writes_chart(tmp_path, monkeypatch):
    monkeypatch.setattr(sequence_analysis, "RESULTS_DIR", tmp_path)
    sequence_analysis.visualize_counts("Test Chart", {"A": 6, "B": 4}, 10)
    assert (tmp_path / "test_chart_distribution.png").exists()


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


def test_derive_transcription_rules_flags_non_bijective_key():
    # 'A' and 'T' both map 100 % to 'Z' -> not a clean one-to-one mapping
    pair_counts = {("A", "Z"): 5, ("T", "Z"): 5}
    _, is_valid = transcription_decoder.derive_transcription_rules(pair_counts)
    assert is_valid is False


def test_check_transcription_accepts_clean_bijection():
    key = {"A": "Z", "T": "U"}
    reverse = {"Z": ["A"], "U": ["T"]}
    assert transcription_decoder.check_transcription(key, reverse) is True


def test_check_transcription_rejects_many_to_one_mapping():
    key = {"A": "Z", "T": "Z"}
    reverse = {"Z": ["A", "T"]}
    assert transcription_decoder.check_transcription(key, reverse) is False


# --- Stage 3: codon length -----------------------------------------------


def test_get_codon_length_without_stop_codon():
    # 3 amino acids, mRNA length 6  ->  codon length 2
    assert codon_length.get_codon_length("ABCDEF", "NaMgCl") == 2


def test_get_codon_length_with_stop_codon():
    # 3 amino acids, mRNA length 8  ->  8 / (3 + 1) = codon length 2
    assert codon_length.get_codon_length("ABCDEFGH", "NaMgCl") == 2


def test_get_codon_length_empty_protein_returns_none():
    assert codon_length.get_codon_length("ABCDEF", "") is None


def test_get_codon_length_returns_none_when_no_integer_length():
    # mRNA length 5, 2 amino acids: 5/2 and 5/3 are both non-integers
    assert codon_length.get_codon_length("ABCDE", "NaMg") is None


# --- Stage 4: codon table ------------------------------------------------


def test_lookup_table_maps_codons_to_amino_acids():
    table = codon_table.lookup_table_mrna_to_amino("UUZX", "AsBe", 2)
    assert table == {"UU": "As", "ZX": "Be"}


def test_lookup_table_records_trailing_stop_codon():
    # the mRNA carries one extra codon ('XX') beyond the protein
    table = codon_table.lookup_table_mrna_to_amino("UUZXXX", "AsBe", 2)
    assert table["XX"] == "STOP"


def test_lookup_table_accepts_float_codon_length():
    # a codon length arriving as a float (e.g. from a division) is normalized
    table = codon_table.lookup_table_mrna_to_amino("UUZX", "AsBe", 2.0)
    assert table == {"UU": "As", "ZX": "Be"}


# --- Stage 5: transcription & translation --------------------------------


def test_transcribe_dna_to_rna_applies_key():
    key = {"A": "Z", "B": "Y", "C": "X", "T": "U"}
    assert decode_gene_b.transcribe_dna_to_rna("ATCB", key) == "ZUXY"


def test_transcribe_marks_unknown_base():
    assert decode_gene_b.transcribe_dna_to_rna("AG", {"A": "Z"}) == "Z?"


def test_transcribe_empty_sequence_returns_empty():
    assert decode_gene_b.transcribe_dna_to_rna("", {"A": "Z"}) == ""


def test_translate_rna_to_protein_builds_sequence():
    table = {"UU": "As", "ZX": "Be"}
    assert decode_gene_b.translate_rna_to_protein("UUZX", 2, table) == "AsBe"


def test_translate_halts_at_stop_codon():
    table = {"UU": "As", "XX": "STOP", "ZZ": "Be"}
    assert decode_gene_b.translate_rna_to_protein("UUXXZZ", 2, table) == "As"


def test_translate_marks_unknown_codon():
    assert decode_gene_b.translate_rna_to_protein("UUZZ", 2, {"UU": "As"}) == "AsX"


def test_translate_empty_sequence_returns_empty():
    assert decode_gene_b.translate_rna_to_protein("", 2, {"UU": "As"}) == ""
