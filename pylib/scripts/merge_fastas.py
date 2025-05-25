#!/usr/bin/env python3

import argparse
import sys
import os
import textwrap
from collections import OrderedDict # To preserve the order of sequences as they are first encountered

# Add the project root directory to the Python path
# This allows imports like 'from pylib.scripts...'
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.insert(0, project_root)
from pylib.utils import seq_parser # For reading FASTA files

def main():
    """
    Merges multiple FASTA files. Sequences with identical headers are concatenated.
    Output is to standard output.
    """
    parser = argparse.ArgumentParser(
        description="Merge multiple FASTA files. Sequences with identical headers are concatenated. "
                    "The order of headers in the output is based on their first appearance."
    )
    parser.add_argument(
        "input_fasta_files",
        nargs='+', # One or more input files
        help="Paths to the input FASTA files."
    )
    args = parser.parse_args()

    merged_sequences = OrderedDict()

    for fasta_file_path in args.input_fasta_files:
        try:
            for record in seq_parser.parse_fasta_file(fasta_file_path):
                # The Perl script uses the full line after '>' as the key.
                # Biopython's record.description contains this, while record.id is usually the first word.
                header_key = record.description if record.description and record.description.strip() else record.id
                
                sequence_part = str(record.seq)
                
                if header_key in merged_sequences:
                    merged_sequences[header_key] += sequence_part
                else:
                    merged_sequences[header_key] = sequence_part
        except FileNotFoundError:
            print(f"Error: Input file not found at {fasta_file_path}", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f"An error occurred while processing {fasta_file_path}: {e}", file=sys.stderr)
            # Decide if we should exit or continue with other files.
            # The Perl one-liner might try to continue, but for robustness, exiting on error is safer.
            sys.exit(1)

    if not merged_sequences:
        # This can happen if all input files were empty or not valid FASTA.
        # Check if at least one file was attempted to differentiate from no files given (handled by nargs='+')
        if args.input_fasta_files:
             print("Warning: No sequences found in any input file or all files were empty/invalid.", file=sys.stderr)
        # No sys.exit(1) here, as outputting nothing for empty inputs might be desired.

    # Print the merged FASTA to standard output
    try:
        for header, sequence in merged_sequences.items():
            print(f">{header}")
            wrapped_sequence = textwrap.wrap(sequence, width=70)
            for line in wrapped_sequence:
                print(line)
    except Exception as e:
        print(f"An error occurred while writing merged FASTA to stdout: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
