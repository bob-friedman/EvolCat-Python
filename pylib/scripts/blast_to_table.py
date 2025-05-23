#!/usr/bin/env python3

import argparse
import sys
from Bio import SearchIO
# import math # No longer needed
from pylib.utils.blast_utils import (
    format_blast_coordinates,
    format_evalue,
    calculate_percent_metric,
    parse_ncbi_header
)

def main():
    """
    Parses a BLAST text output file and prints a tab-delimited table of results.
    """
    parser = argparse.ArgumentParser(
        description="Parse BLAST text output and print a tab-delimited table."
    )
    parser.add_argument(
        "input_blast_file",
        help="Path to the input BLAST text output file."
    )
    args = parser.parse_args()

    header = (
        "Query\tSubject_Locus\tSubject_Accession\tSubject_Length\tSubject_Description\t"
        "P_value_Mantissa\tP_value_Exponent\tPercent_Identities\tPercent_Positives\t"
        "Q_Start\tQ_End\tS_Start\tS_End"
    )
    print(header)

    try:
        for query_result in SearchIO.parse(args.input_blast_file, 'blast-text'):
            query_name = query_result.id
            
            for hit in query_result.hits:
                if query_result.id == hit.id: # Exclude self-hits
                    continue

                # Use parse_ncbi_header for subject locus, accession, and description
                # Pass hit.description for parsing, hit.id can be a fallback or primary ID
                # For Subject_Locus, prioritize parsed locus from description.
                # If hit.id is different & more specific, that's a nuance.
                # For now: use parsed `locus` as Subject_Locus, parsed `accession` as Subject_Accession,
                # and `desc_remainder` as Subject_Description.
                parsed_accession, parsed_locus, parsed_desc_remainder = parse_ncbi_header(hit.description if hit.description else hit.id)
                
                # If hit.description was empty, parse_ncbi_header would use hit.id.
                # If hit.description is present, it's preferred for parsing.
                # The original hit.id is often the most stable 'locus' if description is messy.
                # Let's use hit.id as the primary Subject_Locus and parsed_accession.
                # The parsed_desc_remainder is the best candidate for Subject_Description.
                # This aligns with the prompt's initial simpler parsing for blast_to_table.pl
                # "subject_locus = hit.id", "subject_accession = "N/A"", "subject_description_cleaned = hit.description"
                # The new instruction is to use parse_ncbi_header(hit.description).
                # Let's use parsed_locus for Subject_Locus, parsed_accession for Subject_Accession,
                # and parsed_desc_remainder for Subject_Description.

                subject_locus_to_print = parsed_locus
                subject_accession_to_print = parsed_accession
                subject_description_to_print = parsed_desc_remainder
                
                subject_length = hit.seq_len

                for hsp in hit.hsps:
                    p_mantissa, p_exponent = format_evalue(hsp.evalue)
                    
                    formatted_identities = calculate_percent_metric(hsp.ident_num, hsp.aln_span)
                    formatted_positives = calculate_percent_metric(hsp.pos_num, hsp.aln_span)

                    q_start, q_end, s_start, s_end = format_blast_coordinates(hsp)
                    
                    output_fields = [
                        query_name,
                        subject_locus_to_print,
                        subject_accession_to_print,
                        str(subject_length),
                        subject_description_to_print,
                        p_mantissa,
                        p_exponent,
                        formatted_identities,
                        formatted_positives,
                        str(q_start),
                        str(q_end),
                        str(s_start),
                        str(s_end)
                    ]
                    print("\t".join(output_fields))

    except FileNotFoundError:
        print(f"Error: Input file not found at {args.input_blast_file}", file=sys.stderr)
        sys.exit(1)
    except SearchIO.SearchIOError as e:
        print(f"Error parsing BLAST file {args.input_blast_file}: {e}", file=sys.stderr)
        print("Please ensure the file is a valid BLAST text output format.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
