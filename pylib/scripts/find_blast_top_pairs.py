#!/usr/bin/env python3

import argparse
import sys

def main():
    """
    Processes a tab-delimited input file to find unique top-scoring pairs,
    based on specified columns for QueryID, SubjectID, Score, and a FilterValue.
    """
    parser = argparse.ArgumentParser(
        description="Find unique top-scoring query-subject pairs from tab-delimited BLAST-like results."
    )
    parser.add_argument(
        "input_file",
        nargs='?', # Optional: if not provided, reads from stdin
        help="Path to the input tab-delimited file. Reads from stdin if not specified."
    )
    parser.add_argument(
        "--score_col_idx",
        type=int,
        default=2,
        help="0-based index of the score column (default: 2)."
    )
    parser.add_argument(
        "--filter_col_idx",
        type=int,
        default=3,
        help="0-based index of the column used for filtering (default: 3)."
    )
    parser.add_argument(
        "--filter_threshold",
        type=float,
        default=50.0,
        help="Minimum value for the filter column (default: 50.0)."
    )
    parser.add_argument(
        "--query_id_col_idx",
        type=int,
        default=0,
        help="0-based index of the Query ID column (default: 0)."
    )
    parser.add_argument(
        "--subject_id_col_idx",
        type=int,
        default=1,
        help="0-based index of the Subject ID column (default: 1)."
    )
    args = parser.parse_args()

    hits_data = []
    line_number = 0

    try:
        input_source = open(args.input_file, 'r') if args.input_file else sys.stdin
        
        for line in input_source:
            line_number += 1
            fields = line.rstrip('\n').split('\t')
            
            try:
                query_id = fields[args.query_id_col_idx]
                subject_id = fields[args.subject_id_col_idx]
                
                # Exclude self-hits
                if query_id == subject_id:
                    continue
                
                score_str = fields[args.score_col_idx]
                filter_value_str = fields[args.filter_col_idx]
                
                score = float(score_str)
                filter_value = float(filter_value_str)
                
                # Store enough original fields for potential output (first 4 as per Perl)
                # and the necessary parsed values.
                # The Perl script stores the first 4 fields.
                # We need query_id, subject_id, score, filter_value for logic.
                # We also need the original representation of the first 4 fields for printing.
                original_fields_to_print = fields[:max(4, args.query_id_col_idx+1, args.subject_id_col_idx+1, args.score_col_idx+1, args.filter_col_idx+1)]
                # Ensure at least 4 fields if available, or up to the max index used.
                # For simplicity, if the script is meant to output a fixed number of fields from the original line,
                # this could be simplified. The Perl script takes $fields[0..3].

                hits_data.append({
                    "query_id": query_id,
                    "subject_id": subject_id,
                    "score": score,
                    "filter_value": filter_value,
                    "original_line_fields": fields # Store all original fields
                })

            except IndexError:
                print(f"Warning: Line {line_number} has fewer columns than expected. Skipping.", file=sys.stderr)
                continue
            except ValueError as e:
                print(f"Warning: Could not convert score/filter value on line {line_number} ('{score_str}', '{filter_value_str}'): {e}. Skipping.", file=sys.stderr)
                continue

        if args.input_file:
            input_source.close()

    except FileNotFoundError:
        print(f"Error: Input file not found at {args.input_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred during input processing: {e}", file=sys.stderr)
        sys.exit(1)

    # Sort hits_data in descending order based on Score
    hits_data.sort(key=lambda x: x["score"], reverse=True)

    processed_queries = set()
    processed_subjects = set()

    for hit in hits_data:
        query_id = hit["query_id"]
        subject_id = hit["subject_id"]
        
        if query_id in processed_queries or subject_id in processed_subjects:
            continue
            
        # Apply filter threshold
        if hit["filter_value"] > args.filter_threshold:
            # The Perl script prints the first four original fields.
            # We stored all original fields, so we can reconstruct this.
            # Ensure we don't go out of bounds if original line had < 4 fields.
            fields_to_print = hit["original_line_fields"][:4] 
            # If original line had fewer than 4 fields, this will take what's available.
            # To strictly match Perl's potential for undefined values if original had <4:
            # output_str_fields = []
            # for i in range(4):
            #     output_str_fields.append(hit["original_line_fields"][i] if i < len(hit["original_line_fields"]) else "")
            # print("\t".join(output_str_fields))
            
            # Simpler: print the fields we know exist and were used.
            # The prompt says "print the relevant fields... (e.g. QueryID, SubjectID, Score, FilterValue)"
            # The Perl script prints $line[0..3]. Let's match the Perl script.
            
            # Get the string representation of score and filter_value as they were in the file,
            # or use our float conversion. Perl's output would be from original string.
            # We stored the full original_line_fields.
            
            # Ensure we have at least 4 fields to print, padding with empty strings if not
            # This more closely matches Perl's behavior with `print "$line[0]\t$line[1]\t$line[2]\t$line[3]\n";`
            # where if $line[2] or $line[3] were not there, it would print just a tab or an empty field.
            output_display_fields = []
            for i in range(max(4, args.query_id_col_idx + 1, args.subject_id_col_idx + 1, args.score_col_idx + 1, args.filter_col_idx + 1)):
                 if i < len(hit["original_line_fields"]):
                     output_display_fields.append(hit["original_line_fields"][i])
                 else:
                     output_display_fields.append("") # Pad if original line was too short
            
            # The Perl script specifically prints the first four fields of the original line.
            print("\t".join(output_display_fields[:4]))


            processed_queries.add(query_id)
            processed_subjects.add(subject_id)

if __name__ == "__main__":
    main()
