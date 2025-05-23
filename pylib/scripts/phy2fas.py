#!/usr/bin/env python3

import argparse
import sys
import textwrap
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def main():
    """
    Converts a PHYLIP file (interleaved or sequential) to FASTA format.
    """
    parser = argparse.ArgumentParser(
        description="Convert a PHYLIP file to FASTA format. "
                    "Handles both interleaved and sequential PHYLIP."
    )
    parser.add_argument(
        "input_phylip_file",
        help="Path to the input PHYLIP file."
    )
    args = parser.parse_args()

    try:
        with open(args.input_phylip_file, "r") as handle:
            for record in SeqIO.parse(handle, "phylip"):
                # Construct FASTA header
                header = f">{record.id}"
                if record.description and record.description.strip():
                    # Add description if it's not just whitespace or empty
                    # SeqIO.parse for phylip might put the original name (up to 10 chars) in .id
                    # and potentially other info (if any was parsed from a more complex phylip name line) in .description
                    # Or, .description might just be a copy of .id or empty.
                    # We'll include .description if it's distinct and not empty.
                    # A common case is that record.description will be identical to record.id for PHYLIP.
                    # Let's avoid printing ">ID ID" by checking if description is different from ID.
                    if record.description.strip() != record.id:
                         header += f" {record.description.strip()}"
                
                print(header)
                
                # Print sequence wrapped at 70 characters
                sequence_string = str(record.seq)
                wrapped_sequence = textwrap.wrap(sequence_string, width=70)
                for line in wrapped_sequence:
                    print(line)
                
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
