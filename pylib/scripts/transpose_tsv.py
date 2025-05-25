#!/usr/bin/env python3

import argparse
import sys

def read_data_matrix(lines_iterator):
    """
    Reads lines from an iterator and constructs a matrix (list of lists).
    Each line is split by tabs.
    """
    data_matrix = []
    for line in lines_iterator:
        stripped_line = line.rstrip('\n')
        row = stripped_line.split('\t')
        data_matrix.append(row)
    return data_matrix

def main():
    """
    Transposes a tab-delimited matrix from input file(s) (or stdin).
    """
    parser = argparse.ArgumentParser(
        description="Transpose a tab-delimited matrix. Reads from stdin if no files are specified."
    )
    parser.add_argument(
        "input_files",
        nargs='*', # Zero or more input files
        help="Paths to the input tab-delimited files. If not specified, reads from stdin."
    )
    args = parser.parse_args()

    data_matrix = []

    if not args.input_files: # No files specified, read from stdin
        data_matrix = read_data_matrix(sys.stdin)
    else:
        for file_path in args.input_files:
            try:
                with open(file_path, 'r') as f:
                    # If multiple files, we append rows from each file to the same matrix
                    data_matrix.extend(read_data_matrix(f)) 
            except FileNotFoundError:
                print(f"Error: Input file not found at {file_path}", file=sys.stderr)
                sys.exit(1)
            except Exception as e:
                print(f"An error occurred while processing {file_path}: {e}", file=sys.stderr)
                sys.exit(1)
    
    # Validate Matrix: a. If data_matrix is empty
    if not data_matrix:
        # print("No data provided or all input files were empty.", file=sys.stderr) # Optional message
        sys.exit(0) # Exit cleanly, printing nothing, similar to Perl script with no input.

    # Validate Matrix: b. Determine num_cols and c. Verify consistent column count
    num_cols = len(data_matrix[0])
    for i, row in enumerate(data_matrix):
        if len(row) != num_cols:
            print(f"Error: Row {i+1} (0-indexed {i}) has {len(row)} columns, "
                  f"but expected {num_cols} columns (based on the first row). "
                  "Cannot transpose.", file=sys.stderr)
            sys.exit(1)

    # Transpose the Matrix
    transposed_matrix = []
    if num_cols > 0: # Only proceed if there are columns to transpose
        for j in range(num_cols): # Iterate original columns
            new_row = []
            for i in range(len(data_matrix)): # Iterate original rows
                new_row.append(data_matrix[i][j])
            transposed_matrix.append(new_row)
    
    # Print Transposed Matrix
    for row in transposed_matrix:
        print("\t".join(row))

if __name__ == "__main__":
    main()
