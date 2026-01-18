# -*- coding: utf-8 -*-
"""
Created on Sun Nov 23 23:08:48 2025

@author: weicq
"""
#task3 : Determine the length of the codon
import re
import task1

def get_codon_length(rna_seq, protein_seq):
   
    rna_len = len(rna_seq)
    
    amino_acids = re.findall(r'[A-Z][a-z]?', protein_seq)
    protein_len = len(amino_acids)
    
    if protein_len == 0:
        print("Error: Protein sequence is empty")
        return None
    
    print(f"RNA length: {rna_len}, Protein length {protein_len}")
    
    calculated_codon_len_a = rna_len / protein_len
    calculated_codon_len_b = rna_len / (protein_len + 1)
    
    if calculated_codon_len_a % 1 == 0:
        codon_len = int(calculated_codon_len_a)
        print(f"Result: No termination codon")
        print(f"Codon length is {codon_len}")
        return codon_len
    
    elif calculated_codon_len_b % 1 == 0:
        codon_len = int(calculated_codon_len_b)
        print(f"Result: Have termination codon")
        print(f"Codon length is {codon_len}")
        return codon_len
    
    else:
        print(f"Error: Could not find integer codon length")
        return None
    
def main():
    rna_file_path = 'rna_a.fa'
    protein_file_path = 'protein_a.fa'
    
    print("Task 3: determine codon length")
    
    rna_seq = task1.read_file_sequence(rna_file_path)
    protein_seq = task1.read_file_sequence(protein_file_path)
    
    if rna_seq and protein_seq:
        codon_len = get_codon_length(rna_seq, protein_seq)
        if codon_len:
            print(f"\nThe Martian codon length is {codon_len}")
    else:
        print("Error: Failed to read sequence files")
        
if __name__ == "__main__":
    main()