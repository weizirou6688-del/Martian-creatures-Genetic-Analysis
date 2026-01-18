# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 05:27:39 2025

"""
import re
import task3 #codon_length
import task1
from collections import Counter 

# Determine mrna - amino acids table
# Each codon corresponds to one amino acid
def lookup_table_mrna_to_amino (rna_seq, protein_seq, codon_len):
    #1.prepare empty mrna - amino acids table
    codon_table = {}    
    
    #prepare parameter
    #because the result obtained from division may be a floating-pint numeber
    codon_len = int(codon_len)
    
    #findall function outcome is a list
    amino_acids  =  re.findall(r'[A-Z][a-z]?', protein_seq)
    counts = Counter(amino_acids)
    print(f"The number of amino acid types : {len(counts)}")
    
    #len(amino_acids):
    #due to the findall function return list, so through len() to 
    #visit list, every combining charactor  is correponding index
    expected_rna_len = len(amino_acids) * codon_len
    
    #record every type of amino_acids
    for i in range(len(amino_acids)):
        #obtain the current amino_aicds
        current_aa = amino_acids[i]
        
        #caculate the starting and ending position of every amino_acids
        # +codon_len means to break down amino acids into small nucleotides
        # these correspond exactly to the nucleotides of mrna, so the index match
        start_index =  i * codon_len
        end_index = start_index + codon_len
        
        if end_index <= len(rna_seq):
            codon = rna_seq[start_index : end_index]
            
            if codon not in codon_table:
                codon_table[codon] = current_aa

    
    #check if there is a stop codon
    expected_rna_len = len(amino_acids) * codon_len
    
    #An RNA sequence contains only one stop codon
    #Add a stop to the amino acid;when a stop codon is detected in the mrna,
    #add a stop to the amino acid
    if len(rna_seq) > expected_rna_len:
        #take out the last section
        stop_codon = rna_seq[expected_rna_len : expected_rna_len + codon_len]
        codon_table[stop_codon] = "STOP"
    
    return codon_table

def main():
    print("Task4: Create amino acids and codon lookup table")
    
    rna_file_path = 'rna_a.fa'
    protein_file_path = 'protein_a.fa'
    
    rna_seq = task1.read_file_sequence(rna_file_path)
    protein_seq =  task1.read_file_sequence(protein_file_path)
    
    if rna_seq and protein_seq:
        codon_len = task3.get_codon_length(rna_seq, protein_seq)
        if codon_len:
            print(f"\nThe Martian codon length is {codon_len}")
            table = lookup_table_mrna_to_amino(rna_seq, protein_seq, codon_len)
            print("codon-amino acid query table")
            
            for codon, amino in table.items():
                print(f"{codon} : {amino}")
        else:
            print("Error: Failed to read codon length from task3")
    else:
        print("Error: Failed to read sequence files")

if __name__ == "__main__":
    main()
    
     
     
     
     
     