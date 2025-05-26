"""
Calculates the nucleotide diversity (pi) from a FASTA alignment file.

Nucleotide diversity is defined as the average number of nucleotide
differences per site between any two DNA sequences chosen randomly from the
sample population.

pi = (sum of pairwise differences) / (number of pairs)
Number of pairs = n * (n - 1) / 2, where n is the number of sequences.

Arguments:
  fasta_file (str): Path to the input FASTA alignment file.

Output:
  Prints the calculated nucleotide diversity (pi) to standard output.
"""

import argparse
from Bio import SeqIO
import itertools # Will be used for pairwise combinations

def calculate_pairwise_differences(seq1, seq2):
    """Calculates the number of nucleotide differences between two sequences."""
    differences = 0
    # Ensure sequences are of the same length for a meaningful comparison
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of the same length for comparison.")
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            differences += 1
    return differences

def main():
    parser = argparse.ArgumentParser(description="Calculate nucleotide diversity (pi) from a FASTA alignment.")
    parser.add_argument("fasta_file", help="Path to the input FASTA alignment file.")
    args = parser.parse_args()

    sequences = []
    try:
        for record in SeqIO.parse(args.fasta_file, "fasta"):
            sequences.append(record.seq)
    except FileNotFoundError:
        print(f"Error: File not found at {args.fasta_file}")
        return
    except Exception as e: # Catch other potential parsing errors
        print(f"Error parsing FASTA file: {e}")
        return

    if not sequences:
        print("No sequences found in the FASTA file.")
        return
    
    num_sequences = len(sequences)
    if num_sequences < 2:
        print("Nucleotide diversity calculation requires at least two sequences.")
        # Pi would be 0 or undefined, depending on interpretation.
        # For n=0 or n=1, n*(n-1)/2 is 0.
        # It's clearer to state it requires at least two.
        return

    # Basic check for alignment: are all sequences of the same length?
    first_seq_len = len(sequences[0])
    if not all(len(seq) == first_seq_len for seq in sequences):
        print("Error: Sequences in the FASTA file are not all of the same length. Please provide an aligned FASTA file.")
        return

    total_pairwise_differences = 0
    # Generate all unique pairs of sequences
    for seq1, seq2 in itertools.combinations(sequences, 2):
        total_pairwise_differences += calculate_pairwise_differences(seq1, seq2)

    # Calculate number of pairs
    # n * (n - 1) / 2
    num_pairs = num_sequences * (num_sequences - 1) / 2

    if num_pairs == 0: # Should be caught by num_sequences < 2, but as a safeguard
        nucleotide_diversity = 0.0
    else:
        # The sum of pairwise differences is for the entire length of sequences.
        # Nucleotide diversity (pi) is usually expressed per site.
        # So, we need to divide by the length of the alignment (first_seq_len).
        if first_seq_len == 0: # Avoid division by zero if sequences are empty
             print("Error: Sequences are empty (length 0). Cannot calculate diversity per site.")
             return
        nucleotide_diversity = total_pairwise_differences / num_pairs / first_seq_len


    print(f"Number of sequences: {num_sequences}")
    print(f"Alignment length: {first_seq_len}")
    print(f"Total pairwise differences: {total_pairwise_differences}")
    print(f"Number of pairs: {int(num_pairs)}")
    print(f"Nucleotide diversity (pi): {nucleotide_diversity:.6f}")

if __name__ == "__main__":
    main()
