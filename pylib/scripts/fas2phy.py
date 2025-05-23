#!/usr/bin/env python3

import argparse
import sys
import io
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pylib.utils import seq_parser # Assuming pylib is in PYTHONPATH or script is run from the project root

def main():
    """
    Converts a FASTA file to PHYLIP (sequential) format.
    """
    parser = argparse.ArgumentParser(
        description="Convert a FASTA file to PHYLIP (sequential) format. "
                    "All sequences in the input FASTA file must be of the same length."
    )
    parser.add_argument(
        "input_fasta_file",
        help="Path to the input FASTA file."
    )
    args = parser.parse_args()

    records = []
    try:
        for record in seq_parser.parse_fasta_file(args.input_fasta_file):
            records.append(record)
    except FileNotFoundError:
        print(f"Error: Input file not found at {args.input_fasta_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while reading {args.input_fasta_file}: {e}", file=sys.stderr)
        sys.exit(1)

    if not records:
        print(f"Error: No sequences found in {args.input_fasta_file}", file=sys.stderr)
        sys.exit(1)

    # Biopython's SeqIO.write for phylip expects SeqRecord objects.
    # The name (ID) of the SeqRecord will be used for the PHYLIP name,
    # truncated or padded to 10 characters by Biopython.
    # Let's ensure our records are suitable. seq_parser already returns SeqRecord objects.

    try:
        # Use a StringIO buffer to capture the output, then print to stdout.
        # This can sometimes be more robust than writing directly to sys.stdout
        # with certain library functions or in specific environments.
        output_buffer = io.StringIO()
        SeqIO.write(records, output_buffer, "phylip-sequential")
        print(output_buffer.getvalue(), end='') # end='' because getvalue() will include necessary newlines
        output_buffer.close()
    except ValueError as e:
        # This exception is often raised by Biopython if sequences have different lengths
        # or if names are problematic beyond simple truncation.
        print(f"Error during PHYLIP conversion: {e}", file=sys.stderr)
        print("Please ensure all sequences in the FASTA file are of the same length "
              "and have valid identifiers.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred during PHYLIP conversion: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
