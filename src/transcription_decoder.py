"""
Stage 2 - Discover the Martian transcription key.

A reference DNA gene and its mRNA transcript are aligned, the observed
base pairs are counted, and a candidate DNA -> RNA mapping is derived.
The mapping is then validated to confirm it is strictly bijective
(one-to-one), as required by the central dogma of transcription.
"""

import collections
from pathlib import Path

import sequence_analysis

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


# Slide the RNA along the DNA and return the offset where the conserved
# T-U pairing is best satisfied. In nature an mRNA binds only a fragment
# of a much longer DNA strand, so the alignment start must be discovered.
def get_best_alignment_offset(dna_seq, rna_seq):
    dna_len = len(dna_seq)
    rna_len = len(rna_seq)

    best_offset = 0
    max_matches = -1

    for i in range(dna_len - rna_len + 1):
        # extract the segment of DNA to compare
        dna_segment = dna_seq[i : i + rna_len]

        current_matches = 0
        for d, r in zip(dna_segment, rna_seq):
            if d == 'T' and r == 'U':
                current_matches += 1

        if current_matches > max_matches:
            max_matches = current_matches
            best_offset = i

    print(f"Alignment Found! The starting point of DNA segment begin with {best_offset}")
    return best_offset


# 1. Align the DNA segment to the RNA at the given offset.
# 2. Return pair_counts: {(dna_base, rna_base): pair count}.
def get_aligned_sequences(dna_seq, rna_seq, offset):
    # Line breaks and 'N' padding have already been removed from the data,
    # so no further cleaning is needed here.
    pair_counts = {}
    rna_len = len(rna_seq)

    aligned_dna = dna_seq[offset : offset + rna_len]

    for d, r in zip(aligned_dna, rna_seq):
        pair = (d, r)
        if pair in pair_counts:
            pair_counts[pair] += 1
        else:
            pair_counts[pair] = 1

    return pair_counts


def check_transcription(transcription_key, rna_reverse_pairing):

    is_one_to_one = True

    # each RNA base must pair with exactly one DNA base
    for rna_base, dna_list in rna_reverse_pairing.items():
        if len(dna_list) != 1:
            is_one_to_one = False
            print(f"error:RNA[{rna_base}] is paired with {len(dna_list)} DNA bases: {dna_list} ")

    # check that one DNA base is not pointed to by multiple RNA bases
    dna_to_rna = collections.defaultdict(list)
    for rna_base, dna_list in rna_reverse_pairing.items():
        if len(dna_list) == 1:
            dna_base = dna_list[0]
            dna_to_rna[dna_base].append(rna_base)

    for dna_base, rna_sources in dna_to_rna.items():
        if len(rna_sources) > 1:
            is_one_to_one = False
            print(f"error: DNA{dna_base} is pointed by multiple RNA bases: {rna_sources}")

    # check the round trip A -> B -> A
    for dna_origin, rna_goal in transcription_key.items():
        sources = rna_reverse_pairing.get(rna_goal, [])

        if len(sources) == 1 and sources[0] == dna_origin:
            print(f"DNA[{dna_origin}] -> RNA[{rna_goal}] -> DNA[{sources[0]}]")
        else:
            is_one_to_one = False
            print(f"error: DNA[{dna_origin}] -> RNA[{rna_goal}] -> {sources} (Inconsistency or missing sources)")

    return is_one_to_one


# 1. dna_base_mappings = {dna_base: {rna_base: count}} -> all observed pairings
# 2. transcription_key  = {dna_base: most_frequent_rna} -> the final mapping
# 3. rna_reverse_pairing = {most_frequent_rna: [dna_base]} -> reverse check
def derive_transcription_rules(pair_counts):

    dna_base_mappings = {}
    for (dna_base, rna_base), count in pair_counts.items():
        if dna_base not in dna_base_mappings:
            dna_base_mappings[dna_base] = {}
        dna_base_mappings[dna_base][rna_base] = count

    transcription_key = {}
    rna_reverse_pairing = {}

    # select one-to-one matching bases
    for dna_base in sorted(dna_base_mappings.keys()):
        mappings = dna_base_mappings[dna_base]
        total_count = sum(mappings.values())

        most_frequent_rna = max(mappings, key=mappings.get)
        count = mappings[most_frequent_rna]
        matching_ratio = (count / total_count) * 100

        if count / total_count == 1.0:
            # organize into a simplified dictionary
            transcription_key[dna_base] = most_frequent_rna

        print(f"Transcribed DNA [{dna_base}] into RNA [{most_frequent_rna}] "
              f"(matching_ratio: {matching_ratio:.2f}%)")

        # build the RNA -> DNA reverse pairing for the validity check
        if most_frequent_rna not in rna_reverse_pairing:
            rna_reverse_pairing[most_frequent_rna] = []
        rna_reverse_pairing[most_frequent_rna].append(dna_base)

    is_key_valid = check_transcription(transcription_key, rna_reverse_pairing)

    return transcription_key, is_key_valid


def main():
    dna_seq = sequence_analysis.read_file_sequence(DATA_DIR / 'gene_a.fa')
    rna_seq = sequence_analysis.read_file_sequence(DATA_DIR / 'rna_a.fa')

    if dna_seq and rna_seq:
        best_offset = get_best_alignment_offset(dna_seq, rna_seq)
        pair_counts = get_aligned_sequences(dna_seq, rna_seq, best_offset)
        transcription_key, is_valid = derive_transcription_rules(pair_counts)
        print("\nThe base pairing result is:")
        print(transcription_key)
        print(f"Bijective key valid: {is_valid}")
    else:
        print("error, Fail to read sequence file")


if __name__ == "__main__":
    main()
