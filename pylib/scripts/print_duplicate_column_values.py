#!/usr/bin/env python3

import argparse
import sys
from collections import defaultdict

def process_lines(lines_iterator):
    """
    Processes lines from an iterator, extracts values from the first two columns,
    and prints values that appear more than once in either column (cumulatively).
    """
    column1_values = []
    column2_values = []

    for line in lines_iterator:
        stripped_line = line.rstrip('\n')
        fields = stripped_line.split('\t')
        
        if len(fields) >= 1:
            column1_values.append(fields[0])
        if len(fields) >= 2:
            column2_values.append(fields[1])

    item_counts = defaultdict(int)

    # Process column 1 values
    for item in column1_values:
        item_counts[item] += 1
        if item_counts[item] > 1:
            print(item)
            
    # Process column 2 values, using the same item_counts
    for item in column2_values:
        item_counts[item] += 1
        if item_counts[item] > 1:
            print(item)

def main():
    """
    Reads from input files (or stdin) and prints values from the first two columns
    that appear more than once (cumulatively across both columns).
    """
    parser = argparse.ArgumentParser(
        description="Print values from the first two columns of input file(s) that appear more than once. "
                    "Counts are cumulative across both columns. Reads from stdin if no files are specified."
    )
    parser.add_argument(
        "input_files",
        nargs='*', # Zero or more input files
        help="Paths to the input tab-delimited files. If not specified, reads from stdin."
    )
    args = parser.parse_args()

    if not args.input_files: # No files specified, read from stdin
        if sys.stdin.isatty():
            # Avoid hanging if stdin is a TTY and no data is piped
            # However, the original Perl script would hang, so for direct mimicry,
            # we might allow hanging. For better UX, we could print a help message.
            # Let's stick to closer mimicry for now.
            pass
        process_lines(sys.stdin)
    else:
        for file_path in args.input_files:
            try:
                with open(file_path, 'r') as f:
                    process_lines(f)
            except FileNotFoundError:
                print(f"Error: Input file not found at {file_path}", file=sys.stderr)
                # Decide on behavior: exit or continue with other files.
                # The Perl script would likely die on the first error if files are explicit.
                # If processing multiple files, it's common to report error and continue,
                # but the Perl script doesn't explicitly handle multiple files in its loop.
                # Given it's `while(<>)`, it handles multiple files by concatenating them.
                # So, if one is missing, it might die before processing others.
                # Let's exit on first error for explicit files.
                sys.exit(1)
            except Exception as e:
                print(f"An error occurred while processing {file_path}: {e}", file=sys.stderr)
                sys.exit(1)

if __name__ == "__main__":
    main()
