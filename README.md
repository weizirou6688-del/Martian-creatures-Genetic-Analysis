# Martian Genome Decoder 🧬

## Project Overview
This project is a computational bioinformatics assignment designed to decode the genetic language of a Martian species. Unlike Earth-based biology, the Martian genome utilizes a unique set of nucleotides and a doublet codon structure. 

The goal was to reverse-engineer the transcription and translation rules from a reference gene (Gene A) to decode an unknown sequence (Gene B), utilizing statistical analysis and algorithmic alignment.

## Key Features
- **Dynamic Nucleotide Discovery:** Automatically identifies and quantifies unknown base types (A, B, C, T for DNA; U, X, Y, Z for RNA).
- **Transcription Logic Discovery:** Uses frequency ratio analysis to derive a 100% bijective base-pairing rule.
- **Sliding Window Algorithm:** Implemented a custom alignment algorithm to synchronize DNA and RNA fragments based on conserved T-U pairings.
- **Codon Architecture Analysis:** Mathematically verified a doublet (2-base) codon structure using integer divisibility checks.
- **Data Visualization:** Generates frequency distribution charts using `matplotlib`.

## Tech Stack
- **Language:** Python 3.x
- **Libraries:** `matplotlib`, `re` (Regular Expressions), `collections`
- **Concepts:** Bio-algorithms, Data Cleaning, Statistical Analysis

## Results Summary
| Feature | Earth System | Martian System |
| :--- | :--- | :--- |
| **DNA Bases** | A, T, C, G | A, T, C, B |
| **RNA Bases** | A, U, C, G | U, X, Y, Z |
| **Codon Length** | 3 (Triplet) | 2 (Doublet) |
| **Protein Unit** | Amino Acids | Chemical Elements |

## Project Structure
- `task1-task5.py`: Main script containing the 5-task pipeline.
- `bioinformatic_report.pdf`: Detailed analysis report including biological comparison.

## Author
**Ava Taylor**  
Computational Bioinformatics Student
