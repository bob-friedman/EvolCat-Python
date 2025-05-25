#!/usr/bin/env python3

import argparse
import os
import sys

# Add the project root directory to the Python path
# This allows imports like 'from pylib.scripts...'
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.insert(0, project_root)
from pylib.utils import seq_parser # Assuming pylib is in PYTHONPATH or script is run from the project root

def main():
    """
    Converts a FASTA file to MEGA (.meg) format.
    """
    parser = argparse.ArgumentParser(
        description="Convert a FASTA file to MEGA (.meg) format."
    )
    parser.add_argument(
        "input_fasta_file",
        help="Path to the input FASTA file."
    )
    args = parser.parse_args()

    # Print MEGA header
    print("#mega")
    print("title:")
    print() # Extra newline as per format

    try:
        for record in seq_parser.parse_fasta_file(args.input_fasta_file):
            # Determine the name:
            # Biopython's record.description usually captures the part of the header after the ID.
            # If record.description is empty, record.id is a fallback.
            # The Bio.SeqIO.FastaIO parser sets record.id to the first word after '>',
            # and record.description to the entire line after '>', including the id.
            # So, record.description is closer to the Perl script's behavior of taking the whole line.
            name = record.description if record.description else record.id
            
            # Clean the name
            cleaned_name = name.strip()      # Remove leading/trailing whitespace
            cleaned_name = cleaned_name.replace('.', '') # Remove dots
            cleaned_name = cleaned_name.replace(',', '') # Remove commas
            
            sequence_string = str(record.seq) # Biopython Seq object to string, already no whitespace
            
            # Print sequence in MEGA format
            print(f"#{cleaned_name}")
            print(sequence_string)
            print() # Newline after each sequence block in MEGA format

    except FileNotFoundError:
        print(f"Error: Input file not found at {args.input_fasta_file}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()
