#!/usr/bin/env python3

import argparse
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO # For writing to stdout if needed
from pylib.utils import seq_parser # For reading and possibly writing

def main():
    """
    Extracts a specified region from sequences in a FASTA file.
    """
    parser = argparse.ArgumentParser(
        description="Extract a region from sequences in a FASTA file."
    )
    parser.add_argument(
        "input_fasta_file",
        help="Path to the input FASTA file."
    )
    parser.add_argument(
        "start_bp",
        type=int,
        help="Start base pair of the region to extract (1-based)."
    )
    parser.add_argument(
        "end_bp",
        type=int,
        help="End base pair of the region to extract (1-based)."
    )
    parser.add_argument(
        "--output_file",
        help="Optional path to the output FASTA file. If not provided, output is to stdout."
    )
    args = parser.parse_args()

    if args.start_bp <= 0 or args.end_bp <= 0:
        print("Error: Start and end base pairs must be positive.", file=sys.stderr)
        sys.exit(1)
    
    if args.start_bp > args.end_bp:
        print(f"Error: Start base pair ({args.start_bp}) cannot be greater than end base pair ({args.end_bp}).", file=sys.stderr)
        sys.exit(1)

    new_records = []

    try:
        for record in seq_parser.parse_fasta_file(args.input_fasta_file):
            seq_len = len(record.seq)
            
            # Adjust from 1-based (user) to 0-based (Python slice)
            py_start = args.start_bp - 1
            py_end = args.end_bp # Slicing is exclusive at the end, so end_bp directly works

            # Validate coordinates for the current sequence
            if not (0 <= py_start < seq_len and py_start < py_end and py_end <= seq_len):
                print(f"Warning: Region {args.start_bp}-{args.end_bp} is out of bounds "
                      f"for sequence {record.id} (length {seq_len}). Skipping this record.", file=sys.stderr)
                continue

            sub_sequence = record.seq[py_start:py_end]
            
            new_id = f"{record.id}_region_{args.start_bp}-{args.end_bp}"
            
            original_desc = record.description if record.description and record.description.strip() else record.id
            # Ensure we don't duplicate the ID in the description if record.description was just the ID
            if original_desc == record.id:
                 new_desc = f"{record.id} [extracted {args.start_bp}-{args.end_bp}]"
            elif record.description: # record.description is present and different from record.id
                 new_desc = f"{record.description} [extracted {args.start_bp}-{args.end_bp}]"
            else: # Only record.id was present originally
                 new_desc = f"{record.id} [extracted {args.start_bp}-{args.end_bp}]"


            new_record = SeqRecord(sub_sequence, id=new_id, description=new_desc)
            new_records.append(new_record)

    except FileNotFoundError:
        print(f"Error: Input file not found at {args.input_fasta_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while processing {args.input_fasta_file}: {e}", file=sys.stderr)
        sys.exit(1)

    if not new_records:
        # This can happen if the input file was empty, not FASTA, or all records had invalid regions.
        # If the input file was valid FASTA but no regions could be extracted, a warning might be useful.
        # Check if the input file exists and is not empty to differentiate.
        try:
            with open(args.input_fasta_file, 'r') as f_check:
                if f_check.read(1): # File exists and is not empty
                    print("Warning: No valid regions extracted from any sequence.", file=sys.stderr)
                else: # File is empty
                     print(f"Warning: Input file {args.input_fasta_file} is empty.", file=sys.stderr)
        except FileNotFoundError:
            pass # Already handled
        except Exception:
            pass # Other read errors also handled
        # No sys.exit(1) here, as creating an empty output might be valid if no regions matched.
        
    try:
        if args.output_file:
            seq_parser.write_fasta_file(new_records, args.output_file)
        else:
            SeqIO.write(new_records, sys.stdout, "fasta")
    except Exception as e:
        print(f"An error occurred while writing output: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
