#!/usr/bin/env python3

import argparse
from pylib.utils import seq_parser # Assuming pylib is in PYTHONPATH or script is run from the project root

def main():
    """
    Converts a FASTA file to a CSV file (name,sequence).
    """
    parser = argparse.ArgumentParser(
        description="Convert a FASTA file to a CSV file, with each line containing 'name,sequence'."
    )
    parser.add_argument(
        "input_fasta_file",
        help="Path to the input FASTA file."
    )
    args = parser.parse_args()

    try:
        for record in seq_parser.parse_fasta_file(args.input_fasta_file):
            # Determine the name:
            # The Perl script uses everything after '>' in the header line as the name.
            # Biopython's record.description usually captures the part of the header after the ID.
            # If record.description is empty, record.id is a fallback.
            # A more direct equivalent to the Perl script is to use record.description if available,
            # as it often contains the full title line without the ID part.
            # If record.description is available and non-empty, use it.
            # Otherwise, use record.id.
            # The Bio.SeqIO.FastaIO parser sets record.id to the first word after '>',
            # and record.description to the entire line after '>', including the id.
            # So, record.description is closer to the Perl script's behavior.

            name = record.description if record.description else record.id
            
            sequence_string = str(record.seq)
            
            # Print in "name,sequence" format
            print(f"{name},{sequence_string}")

    except FileNotFoundError:
        print(f"Error: Input file not found at {args.input_fasta_file}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()
