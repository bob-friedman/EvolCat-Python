#!/usr/bin/env python3

import argparse
import re
import sys
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO # For writing to stdout

# Add the project root directory to the Python path
# This allows imports like 'from pylib.scripts...'
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.insert(0, project_root)
from pylib.utils import seq_parser # For reading

def main():
    """
    Cleans FASTA headers: replaces whitespace and periods with underscores,
    and converts to uppercase. Writes modified FASTA to standard output.
    """
    parser = argparse.ArgumentParser(
        description="Clean FASTA headers: replace whitespace and periods with underscores, "
                    "and convert to uppercase. Output to stdout."
    )
    parser.add_argument(
        "input_fasta_file",
        help="Path to the input FASTA file."
    )
    args = parser.parse_args()

    modified_records = []

    try:
        for record in seq_parser.parse_fasta_file(args.input_fasta_file):
            # Start with the record.description. If it's empty or None, use record.id.
            # Biopython's FastaIO parser puts the whole line after '>' into record.description,
            # and the first word into record.id. Using record.description is closer to the Perl script.
            header_base = record.description if record.description and record.description.strip() else record.id
            
            # Apply cleaning transformations
            # 1. Replace all occurrences of one or more whitespace characters with a single underscore.
            cleaned_header = re.sub(r'\s+', '_', header_base)
            # 2. Replace all period characters (.) with underscores (_).
            cleaned_header = cleaned_header.replace('.', '_')
            # 3. Convert the entire string to uppercase.
            cleaned_header = cleaned_header.upper()
            
            # Create a new SeqRecord.
            # The id of this new record should be the cleaned string.
            # The description field of the new record can be set to an empty string.
            # The sequence (record.seq) remains unchanged.
            new_record = SeqRecord(record.seq, id=cleaned_header, description="")
            modified_records.append(new_record)

    except FileNotFoundError:
        print(f"Error: Input file not found at {args.input_fasta_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while processing {args.input_fasta_file}: {e}", file=sys.stderr)
        sys.exit(1)

    if not modified_records and args.input_fasta_file: # Check if file was processed but no records found
        try: # Check if file exists but is empty or not FASTA
            with open(args.input_fasta_file, 'r') as f:
                if not f.read(1):
                    print(f"Warning: Input file {args.input_fasta_file} is empty.", file=sys.stderr)
        except: # If file not found, already handled
            pass
    
    # Use Bio.SeqIO.write to write all modified SeqRecord objects to standard output.
    # The `write_fasta_file` in `pylib.utils.seq_parser` takes a file path.
    # For writing to stdout, it's cleaner to use SeqIO.write directly here.
    try:
        SeqIO.write(modified_records, sys.stdout, "fasta")
    except Exception as e:
        # This might happen if there's an issue with the records list or stdout, though less common.
        print(f"An error occurred while writing FASTA to stdout: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
