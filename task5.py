# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 06:55:06 2025

"""

"""
Task5: transcription and translation of gene b
"""

import task1
import task2
import task3
import task4

# derive the information of transcription key and codon table lookup table from gene_a
def generate_basic_information():
    
    #1.read sequences
    gene_a_path = 'gene_a.fa'
    rna_a_path = 'rna_a.fa'
    protein_a_path = 'protein_a.fa'
    
   
    gene_a = task1.read_file_sequence(gene_a_path)
    rna_a = task1.read_file_sequence(rna_a_path)
    protein_a = task1.read_file_sequence(protein_a_path)
    
    if not (gene_a and rna_a and protein_a):
        print("Error:Failed to read Gene B files")
        return None, None, None
    
    #2.get transcription_key
    print("Deriving transcription key")
    offset = task2.get_best_alignment_offset(gene_a, rna_a)
    pair_counts = task2.get_aligned_sequences(gene_a, rna_a, offset)
    transcription_key, _ = task2.derive_transcription_rules(pair_counts)
    
    #3.get codon length
    print("calculating Codon length")
    codon_len = task3.get_codon_length(rna_a, protein_a)
    
    #4.get codon table
    print("generating the codon lookup table")
    codon_table = task4.lookup_table_mrna_to_amino(rna_a, protein_a, codon_len)
    
    return transcription_key, codon_len, codon_table

def transcribe_dna_to_rna(dna_seq, transcription_key):
    
    rna_seq = []
    for base in dna_seq:
        rna_base = transcription_key.get(base, '?')
        rna_seq.append(rna_base)
        
    return "".join(rna_seq)

def translate_rna_to_protein(rna_seq, codon_len, codon_table):
    protein_seq = []
    
    #iterate through the RNA sequence with step = codon_len
    for i in range(0, len(rna_seq), codon_len):
        start = i
        end = i + codon_len
        
        if end > len(rna_seq):
            break
        
        codon = rna_seq[start : end]
        
        amino_acid = codon_table.get(codon, 'X')
        
        if amino_acid == "STOP":
            print(f" Termination codon detected at index {start}。 Translation Stop")
            break
        
        protein_seq.append(amino_acid)
        
    return "".join(protein_seq)

def main():
    print("task5 : Gene B analysis")
    key, length, table = generate_basic_information()
    
    if not (key and length and table):
        print("error: could not derive transcripted and translated data now")
        return
    
    gene_b_path = 'gene_b.fa'
    gene_b_seq = task1.read_file_sequence(gene_b_path)

    if not gene_b_seq:
        print(f"error: couldo not read {gene_b_path}")
        return

    print(f"\nGene b length: {len(gene_b_seq)}\n")
    
    print("transcribing Gene B to mRNA")
    mrna_b_seq = transcribe_dna_to_rna(gene_b_seq, key)
    print(f"\n Generated mRNA sequence {mrna_b_seq}\n")
    
    print("Translate mrna to protein\n")
    protein_b_seq = translate_rna_to_protein(mrna_b_seq, length, table)
    
    print(protein_b_seq)
    
if __name__ == "__main__":
    main()