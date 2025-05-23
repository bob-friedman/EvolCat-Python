#!/usr/bin/env python3

import argparse
import sys
# import re # No longer needed directly if is_valid_dna_strict_regex handles it
from collections import Counter
from pylib.utils import seq_parser # For reading FASTA files
from pylib.utils.sequence_utils import generate_kmers, is_valid_dna_strict_regex

def main():
    """
    Counts k-mer frequencies in sequences from a FASTA file.
    """
    parser = argparse.ArgumentParser(
        description="Count k-mer frequencies in sequences from a FASTA file. "
                    "Sequences with non-DNA characters (A,T,C,G) are skipped."
    )
    parser.add_argument(
        "input_fasta_file",
        help="Path to the input FASTA file."
    )
    parser.add_argument(
        "--kmer_len",
        type=int,
        default=5,
        help="Length of the k-mers to count (default: 5)."
    )
    args = parser.parse_args()

    if args.kmer_len <= 0:
        print("Error: --kmer_len must be a positive integer.", file=sys.stderr)
        sys.exit(1)

    kmer_counts = Counter()

    try:
        for record in seq_parser.parse_fasta_file(args.input_fasta_file):
            sequence_str = str(record.seq).upper()
            
            # Validate sequence: Check for non-DNA characters
            if not is_valid_dna_strict_regex(sequence_str):
                print(f"Warning: Sequence {record.id} contains non-DNA characters (other than A, T, C, G). Skipping.", 
                      file=sys.stderr)
                continue
            
            # Extract and count k-mers using the utility function
            # The generate_kmers function already handles the case where sequence_str is shorter than kmer_len.
            for kmer in generate_kmers(sequence_str, args.kmer_len):
                kmer_counts[kmer] += 1
            
            # Optional: Warning if a sequence is shorter than kmer_len
            # if len(sequence_str) < args.kmer_len:
                # print(f"Warning: Sequence {record.id} (length {len(sequence_str)}) is shorter than k-mer length ({args.kmer_len}). No k-mers extracted.", 
                #       file=sys.stderr)


    except FileNotFoundError:
        print(f"Error: Input file not found at {args.input_fasta_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e: # Catches other errors like parsing issues
        print(f"An error occurred while processing {args.input_fasta_file}: {e}", file=sys.stderr)
        sys.exit(1)

    # Sort k-mers: primarily by frequency (descending), secondarily by k-mer (alphabetically)
    # Counter.most_common() sorts by frequency (desc) then by insertion order for ties.
    # To get alphabetical secondary sort for ties, we need a custom sort.
    # Step 1: Get items
    # Step 2: Sort them. Python's default sort is stable.
    # Sort by k-mer alphabetically first (ascending).
    # Then sort by frequency (descending).
    sorted_kmer_counts = sorted(kmer_counts.items(), key=lambda item: item[0]) # Sort by k-mer
    sorted_kmer_counts.sort(key=lambda item: item[1], reverse=True) # Sort by frequency (desc)

    # Print k-mer and count
    for kmer, count in sorted_kmer_counts:
        print(f"{kmer}:{count}")

if __name__ == "__main__":
    main()
