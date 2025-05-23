#!/usr/bin/env python3

import argparse
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO # For writing to stdout
from pylib.utils import seq_parser # For reading

def main():
    """
    Generates the reverse complement of sequences from a FASTA file.
    Output is to standard output in FASTA format.
    """
    parser = argparse.ArgumentParser(
        description="Generate the reverse complement of sequences from a FASTA file. "
                    "Output is to standard output in FASTA format."
    )
    parser.add_argument(
        "input_fasta_file",
        help="Path to the input FASTA file."
    )
    args = parser.parse_args()

    new_records = []

    try:
        for original_record in seq_parser.parse_fasta_file(args.input_fasta_file):
            rev_comp_seq = original_record.seq.reverse_complement()
            
            new_id = f"{original_record.id}_rc"
            
            # Construct new description
            # Use original_record.description if it exists and is not just the ID itself.
            # Otherwise, base the description on the new ID.
            original_desc_text = original_record.description
            if original_desc_text and original_desc_text.strip() and original_desc_text.strip() != original_record.id:
                new_desc = f"{original_desc_text} (reverse complement)"
            else: # Original description was empty, whitespace, or same as ID
                new_desc = f"{new_id} (reverse complement)"

            new_record = SeqRecord(rev_comp_seq, id=new_id, description=new_desc)
            new_records.append(new_record)

    except FileNotFoundError:
        print(f"Error: Input file not found at {args.input_fasta_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e: # Catches errors from parsing or reverse_complement
        print(f"An error occurred while processing {args.input_fasta_file}: {e}", file=sys.stderr)
        sys.exit(1)

    if not new_records and args.input_fasta_file:
        try:
            with open(args.input_fasta_file, 'r') as f:
                if not f.read(1):
                    print(f"Warning: Input file {args.input_fasta_file} is empty.", file=sys.stderr)
        except: 
            pass # File not found already handled, other read errors also handled
    
    try:
        SeqIO.write(new_records, sys.stdout, "fasta")
    except Exception as e:
        print(f"An error occurred while writing FASTA to stdout: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
