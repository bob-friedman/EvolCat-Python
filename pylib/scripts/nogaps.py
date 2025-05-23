#!/usr/bin/env python3

import argparse
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO # For writing to stdout
from pylib.utils import seq_parser # For reading

def main():
    """
    Removes columns containing non-alphabetic characters (A-Z) from an aligned FASTA file.
    Requires all input sequences to be of the same length.
    Output is to standard output in FASTA format.
    """
    parser = argparse.ArgumentParser(
        description="Remove columns with non-alphabetic characters (gaps, etc.) from an aligned FASTA file. "
                    "All sequences must be of the same length. Output is to standard output."
    )
    parser.add_argument(
        "input_fasta_file",
        help="Path to the input aligned FASTA file."
    )
    args = parser.parse_args()

    records = []
    try:
        for record in seq_parser.parse_fasta_file(args.input_fasta_file):
            records.append(record)
    except FileNotFoundError:
        print(f"Error: Input file not found at {args.input_fasta_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e: # Catches other parsing errors
        print(f"Error reading or parsing FASTA file {args.input_fasta_file}: {e}", file=sys.stderr)
        sys.exit(1)

    # Validate Input: a. If no sequences are found
    if not records:
        print("Error: No sequences found in the input file.", file=sys.stderr)
        sys.exit(1)

    # Validate Input: b. Check if all sequences have the same length
    first_seq_len = len(records[0].seq)
    for i, record in enumerate(records):
        if len(record.seq) != first_seq_len:
            print(f"Error: Sequences are not of equal length. "
                  f"Sequence '{records[0].id}' has length {first_seq_len}, "
                  f"but sequence '{record.id}' (index {i}) has length {len(record.seq)}. "
                  "Aligned sequences of equal length are required.", file=sys.stderr)
            sys.exit(1)
    
    alignment_length = first_seq_len
    if alignment_length == 0: # All sequences are empty, nothing to do
        print("Warning: Input sequences are empty. Output will also be empty.", file=sys.stderr)
        SeqIO.write([], sys.stdout, "fasta") # Write an empty fasta
        sys.exit(0)


    # Identify Gap Columns
    is_valid_column = [True] * alignment_length
    for j in range(alignment_length): # Iterate columns
        for record in records: # Iterate sequences for the current column
            char = str(record.seq[j]).upper() # Get char and convert to uppercase
            if not ('A' <= char <= 'Z'): # If char is not an uppercase letter
                is_valid_column[j] = False
                break # Column j is invalid, no need to check other sequences for this column

    # Construct and Output New Sequences
    new_records = []
    for original_record in records:
        new_sequence_chars = []
        for j in range(alignment_length):
            if is_valid_column[j]:
                new_sequence_chars.append(str(original_record.seq[j]))
        
        new_sequence_str = "".join(new_sequence_chars)
        
        # Create a new SeqRecord. Preserve original ID and description.
        new_record = SeqRecord(Seq(new_sequence_str), 
                               id=original_record.id, 
                               description=original_record.description)
        new_records.append(new_record)

    try:
        SeqIO.write(new_records, sys.stdout, "fasta")
    except Exception as e:
        print(f"An error occurred while writing FASTA to stdout: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
