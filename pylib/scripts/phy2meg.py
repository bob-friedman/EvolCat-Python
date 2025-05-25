#!/usr/bin/env python3

import argparse
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def main():
    """
    Converts a PHYLIP file (interleaved or sequential) to MEGA (.meg) format.
    """
    parser = argparse.ArgumentParser(
        description="Convert a PHYLIP file to MEGA (.meg) format. "
                    "Handles both interleaved and sequential PHYLIP."
    )
    parser.add_argument(
        "input_phylip_file",
        help="Path to the input PHYLIP file."
    )
    args = parser.parse_args()

    # Print MEGA header
    print("#mega")
    print("title:")
    print() # Extra newline as per format

    try:
        with open(args.input_phylip_file, "r") as handle:
            for record in SeqIO.parse(handle, "phylip"):
                # Biopython's SeqIO.parse for "phylip" format reads the sequence name (typically up to 10 chars)
                # into record.id. The Perl script uses the first part of the line as $name[$_].
                # record.id from Biopython is the direct equivalent.
                # The Perl script does no further cleaning of this name for phy2meg.pl.
                
                name = record.id
                sequence_string = str(record.seq) # Biopython already concatenates interleaved sequences.
                
                # Print sequence in MEGA format
                # The Perl script adds a newline after each sequence block.
                print(f"#{name}")
                print(sequence_string)
                print() # Newline after each sequence block

    except FileNotFoundError:
        print(f"Error: Input file not found at {args.input_phylip_file}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e: # Biopython often raises ValueError for format issues
        print(f"Error parsing PHYLIP file {args.input_phylip_file}: {e}", file=sys.stderr)
        print("Please ensure the file is a valid PHYLIP format.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
