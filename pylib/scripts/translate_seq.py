#!/usr/bin/env python3

import argparse
from Bio.Seq import Seq
from Bio.Data import CodonTable
from pylib.utils import seq_parser # Assuming pylib is in PYTHONPATH or script is run from parent dir

def main():
    """
    Translates nucleotide sequences from a FASTA file to protein sequences.
    """
    parser = argparse.ArgumentParser(
        description="Translate nucleotide sequences from a FASTA file to protein sequences."
    )
    parser.add_argument(
        "--input_file",
        required=True,
        help="Path to the input FASTA file."
    )
    parser.add_argument(
        "--frame",
        type=int,
        default=1,
        choices=[1, 2, 3],
        help="Translation frame (1, 2, or 3). Default is 1."
    )
    parser.add_argument(
        "--table",
        type=int,
        default=1,
        help="NCBI translation table ID (integer). Default is 1 (standard genetic code)."
    )
    args = parser.parse_args()

    try:
        # Check if the NCBI table ID is valid before parsing the file
        # This will raise an exception if the table is not found
        CodonTable.ambiguous_dna_by_id[args.table]
    except KeyError:
        print(f"Error: Invalid NCBI translation table ID: {args.table}")
        print(f"Available table IDs can be found at https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi")
        return
    except Exception as e:
        print(f"An unexpected error occurred with table ID {args.table}: {e}")
        return


    try:
        for record in seq_parser.parse_fasta_file(args.input_file):
            if not record.seq:
                print(f"Warning: Record {record.id} is empty. Skipping.")
                continue

            # Adjust sequence for the frame. Frame is 1, 2, or 3.
            # Slicing is 0-indexed, so frame 1 is seq[0:], frame 2 is seq[1:], etc.
            framed_seq = record.seq[args.frame - 1:]

            if not framed_seq:
                print(f"Warning: Record {record.id} with frame {args.frame} results in an empty sequence. Skipping.")
                continue
            
            try:
                # Translate the sequence. to_stop=True stops at the first stop codon.
                protein_seq = framed_seq.translate(table=args.table, to_stop=True)
            except CodonTable.TranslationError as e:
                print(f"Warning: Could not translate record {record.id} (frame {args.frame}, table {args.table}): {e}. Skipping.")
                continue


            # Output FASTA format
            output_id = f"{record.id}_translated_frame_{args.frame}_table_{args.table}"
            print(f">{output_id}")
            
            # Print sequence wrapped at a standard width (e.g., 70 characters)
            seq_str = str(protein_seq)
            for i in range(0, len(seq_str), 70):
                print(seq_str[i:i+70])

    except FileNotFoundError:
        print(f"Error: Input file not found at {args.input_file}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()
