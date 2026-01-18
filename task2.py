# Task2 : discover transcription key

import task1
import collections

'''
We’re running a sliding matching process here to 
locate the sequence start.
 In practice, RNA binds to a fragment of a long DNA
 sequence—but since there are no known pairing rules to follow,
 I don’t know how to code this logic.
'''
#find the initial position of sliding match
def get_best_alignment_offset(dna_seq, rna_seq):
    dna_len = len(dna_seq)
    rna_len = len(rna_seq)
    
    best_offset = 0
    max_matches = -1

    for i in range(dna_len - rna_len +1):
        #Extract the segment of DNA to compare
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
        
#1.sequence alignment ——》 do the length matching
#2.return pair_counts: {(dna_base, rna_base): pair counts}
def get_aligned_sequences(dna_seq, rna_seq, offset):
    #Since line breaks and 'N' in the sequence have already been removed from the
    #data,we will not handling sequence here.
       
    pair_counts = {}
    rna_len= len(rna_seq)
    
    aligned_dna = dna_seq[offset : offset + rna_len]
    
    for d, r in zip(aligned_dna, rna_seq):
        pair = (d,r)
        if pair in pair_counts:
            pair_counts[pair] += 1
        else:
            pair_counts[pair] =1
    
    return pair_counts

def check_transcription(transcription_key, rna_reverse_pairing):
    
    is_one_to_one = True
    
    #Each RNA base pairs only one DNA base
    for rna_base, dna_list in rna_reverse_pairing.items():
        if len(dna_list) != 1:
            is_one_to_one = False
            print(f"error:RNA[{rna_base}] is paired with {len(dna_list)} DNA bases: {dna_list} ")
    
    #Check one RNA base only point to one type of DNA base
    dna_to_rna = collections.defaultdict(list)
    for rna_base, dna_list in rna_reverse_pairing.items():
        if len(dna_list) == 1:
            dna_base = dna_list[0]
            dna_to_rna[dna_base].append(rna_base)
            
    for dna_base, rna_sources in dna_to_rna.items():
        if len(rna_sources) > 1:
            is_one_to_one = False
            print(f"error: DNA{dna_base} is pointed by multiple RNA bases: {rna_sources}")
    
    #check A->B, B->A
    for dna_origin, rna_goal in transcription_key.items():
        sources = rna_reverse_pairing.get(rna_goal, [])
        
        if len(sources) == 1 and sources[0] == dna_origin:
            print(f"DNA[{dna_origin}] -> RNA[{rna_goal}] -> DNA[{sources[0]}]")
        else:
            is_one_to_one = False
            print(f"error: DNA[{dna_origin}] -> RNA[{rna_goal}] -> {sources} (Inconsistency or missing sources)")
     
    return is_one_to_one

#1.dna_base_mappings = {dna_base : {rna_base : count}} -> identified the paired bases, mappings (inside dna_base_mappings)= {rna_base : count}
#2.transcription_key = {dna_base: most_frequent_rna} -> sort out the final paired bases
#3.rna_reverse_pairing = {most_frequence_rna : dna_base} -> RNA's reverse pairing check with DNA
def derive_transcription_rules(pair_counts):
    
    dna_base_mappings = {}
    for (dna_base, rna_base), count in pair_counts.items():
        if dna_base not in dna_base_mappings:
            dna_base_mappings[dna_base] = {}
        dna_base_mappings[dna_base][rna_base] = count

    transcription_key = {}
    rna_reverse_pairing = {}

    #Select one-on-one matching bases
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
        
        #2.rna reverse complement dictionary
        #rna_reverse_pairing{most_frequent_rna : dna_base}
        if most_frequent_rna not in rna_reverse_pairing:
            rna_reverse_pairing[most_frequent_rna] = []
        rna_reverse_pairing[most_frequent_rna].append(dna_base)
    
    #call the check function
    is_key_valid = check_transcription(transcription_key, rna_reverse_pairing)
    
    return transcription_key, is_key_valid
    


def main():
    dna_file_path = 'gene_a.fa'
    rna_file_path = 'rna_a.fa'
    
    dna_seq = task1.read_file_sequence(dna_file_path)
    rna_seq = task1.read_file_sequence(rna_file_path)
    
    if dna_seq and rna_seq:
        best_offset = get_best_alignment_offset(dna_seq, rna_seq)
        pair_counts = get_aligned_sequences(dna_seq, rna_seq, best_offset)
        final_key = derive_transcription_rules(pair_counts)
        print("The base pairing result is")
        print(final_key)
    else:
        print("error, Fail to read sequence file")
        
if __name__ ==  "__main__":
    main()