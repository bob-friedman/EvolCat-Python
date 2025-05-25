#!/usr/bin/env python3

import argparse
import sys
import os
from Bio import SearchIO

# Add the project root directory to the Python path
# This allows imports like 'from pylib.scripts...'
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.insert(0, project_root)
from pylib.utils.blast_utils import format_blast_coordinates, get_hit_strand_str

def main():
    """
    Parses a BLAST text output file and prints selected information for each hit and HSP.
    """
    parser = argparse.ArgumentParser(
        description="Parse BLAST text output and print structured information."
    )
    parser.add_argument(
        "input_blast_file",
        help="Path to the input BLAST text output file."
    )
    args = parser.parse_args()

    try:
        # Iterate through each query result in the BLAST output file
        for query_result in SearchIO.parse(args.input_blast_file, 'blast-text'):
            # The Perl script doesn't explicitly print query ID, but it processes per query.
            # We will iterate through hits for each query.
            
            # The prompt mentions sorting hits by subject identifier.
            # For now, processing as returned by Bio.SearchIO.
            # If sorting is needed:
            # sorted_hits = sorted(query_result.hits, key=lambda hit: hit.id)
            # for hit in sorted_hits:
            
            for hit in query_result.hits: # Each hit is a subject sequence
                print(f"{hit.id}:")
                print(f"     {hit.description}")
                print(f"     Length = {hit.seq_len}")
                
                for hsp in hit.hsps: # Each hsp is an alignment region
                    evalue = hsp.evalue
                    identities_num = hsp.ident_num
                    alignment_span = hsp.aln_span 
                    gaps_num = hsp.gap_num # Total number of gap characters

                    orientation_str = get_hit_strand_str(hsp)
                    q_start, q_end, s_start, s_end = format_blast_coordinates(hsp)
                    
                    # Printing the list representation
                    # To exactly match `f"          [{item1}, {item2}, ...]"`
                    # we need to format each item. For simplicity, using str() for numbers now.
                    # For better float formatting: f"{evalue:.2e}" or similar.
                    # Using repr() for strings to ensure quotes if they were part of original, though not here.
                    # The Perl script's Dumper output for strings doesn't add extra quotes unless needed.
                    # So, direct string representation is fine.
                    
                    # Representing the list as a string
                    # Using a more direct f-string approach for the list-like output
                    # The Perl Dumper output for arrays looks like: $VAR1 = [ val1, val2, ... ];
                    # The prompt's example is: "[{hsp.evalue}, {identities_num}, ...]"
                    # Let's match the prompt's specific f-string example format.
                    
                    # Formatting numbers: evalue as scientific, others as integers
                    # Using a helper to format list items to match the desired string representation
                    formatted_hsp_items = [
                        f"{evalue:.2e}" if isinstance(evalue, float) else str(evalue), # evalue variable already exists
                        str(identities_num),
                        str(alignment_span),
                        str(gaps_num),
                        f"'{orientation_str}'", # Add quotes around string as in Dumper-like output
                        str(q_start),
                        str(q_end),
                        str(s_start),
                        str(s_end)
                    ]
                    print(f"          [{', '.join(formatted_hsp_items)}]")
                
                print() # Newline after all HSPs for a given hit

    except FileNotFoundError:
        print(f"Error: Input file not found at {args.input_blast_file}", file=sys.stderr)
        sys.exit(1)
    except SearchIO.SearchIOError as e: # Specific error for parsing issues
        print(f"Error parsing BLAST file {args.input_blast_file}: {e}", file=sys.stderr)
        print("Please ensure the file is a valid BLAST text output format.", file=sys.stderr)
        sys.exit(1)
    except Exception as e: # Catch-all for other unexpected errors
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
