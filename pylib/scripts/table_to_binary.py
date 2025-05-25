#!/usr/bin/env python3

import argparse
import sys

def process_lines_for_binary_conversion(lines_iterator):
    """
    Reads lines from an iterator, converts each element in a tab-separated line
    to "1" if it's a number >= 1.0, otherwise to "0".
    Prints the resulting binary row.
    """
    # warning_issued = False # Optional: for printing warning only once

    for line in lines_iterator:
        stripped_line = line.rstrip('\n')
        elements = stripped_line.split('\t')
        binary_row = []
        
        for element in elements:
            try:
                value = float(element)
                if value >= 1.0:
                    binary_row.append("1")
                else:
                    binary_row.append("0")
            except ValueError:
                # Treat non-numeric values as 0
                binary_row.append("0")
                # Optional: Print warning for non-numeric data
                # if not warning_issued:
                #     print("Warning: Non-numeric data encountered, treating as 0.", file=sys.stderr)
                #     warning_issued = True
        
        print("\t".join(binary_row))

def main():
    """
    Reads tab-delimited data from input file(s) (or stdin), converts values
    to binary (1 if >= 1.0, else 0), and prints the binary rows.
    """
    parser = argparse.ArgumentParser(
        description="Convert tab-delimited table data to binary (1 if element >= 1.0, else 0). "
                    "Non-numeric elements are treated as 0. Reads from stdin if no files are specified."
    )
    parser.add_argument(
        "input_files",
        nargs='*', # Zero or more input files
        help="Paths to the input tab-delimited files. If not specified, reads from stdin."
    )
    args = parser.parse_args()

    if not args.input_files: # No files specified, read from stdin
        process_lines_for_binary_conversion(sys.stdin)
    else:
        for file_path in args.input_files:
            try:
                with open(file_path, 'r') as f:
                    process_lines_for_binary_conversion(f)
            except FileNotFoundError:
                print(f"Error: Input file not found at {file_path}", file=sys.stderr)
                sys.exit(1)
            except Exception as e:
                print(f"An error occurred while processing {file_path}: {e}", file=sys.stderr)
                sys.exit(1)

if __name__ == "__main__":
    main()
