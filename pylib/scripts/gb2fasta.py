#!/usr/bin/env python3

import argparse
import sys
import os

# Add the project root directory to the Python path
# This allows imports like 'from pylib.scripts...'
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.insert(0, project_root)
from pylib.utils import seq_parser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main():
    """
    Converts a GenBank file to FASTA format.
    """
    parser = argparse.ArgumentParser(description="Convert a GenBank file to FASTA format.")
    parser.add_argument("genbank_file", help="Path to the input GenBank file.")
    args = parser.parse_args()

    try:
        for record in seq_parser.parse_genbank_file(args.genbank_file):
            # Construct FASTA header
            # Biopython's SeqRecord.format("fasta") already creates a header like ">id description"
            # so we can directly use that.
            # However, the Perl script uses ">gb:LOCUS_NAME DEFINITION"
            # record.id in Biopython usually corresponds to LOCUS.
            # record.description corresponds to DEFINITION.
            # The "gb:" prefix seems specific to the Perl script's output.
            # We will stick to ">record.id record.description" as per instruction.

            # Format the record as FASTA. This includes the header and wrapped sequence.
            # SeqRecord.format("fasta") handles wrapping at 60 characters by default.
            # To ensure 70 characters, we'll have to do it manually or use a Bio.SeqIO.FastaIO.FastaWriter
            # For simplicity and to match the 70 char requirement, let's format manually.

            print(f">{record.id} {record.description}")
            seq_str = str(record.seq)
            for i in range(0, len(seq_str), 70):
                print(seq_str[i:i+70])

    except FileNotFoundError:
        print(f"Error: File not found at {args.genbank_file}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
