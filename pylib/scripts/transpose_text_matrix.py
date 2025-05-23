#!/usr/bin/env python3

import argparse
import sys

def read_character_matrix(lines_iterator):
    """
    Reads lines from an iterator and constructs a character matrix (list of lists of characters).
    Each line is stripped and then converted into a list of its characters.
    """
    character_matrix = []
    for line in lines_iterator:
        line_content = line.rstrip('\n')
        character_matrix.append(list(line_content)) # Convert line to list of characters
    return character_matrix

def main():
    """
    Transposes a matrix of characters from input file(s) (or stdin).
    Each line in the input is a row, each character an element.
    """
    parser = argparse.ArgumentParser(
        description="Transpose a character matrix (each line is a row, each char an element). "
                    "Reads from stdin if no files are specified."
    )
    parser.add_argument(
        "input_files",
        nargs='*', # Zero or more input files
        help="Paths to the input files. If not specified, reads from stdin."
    )
    args = parser.parse_args()

    character_matrix = []

    if not args.input_files: # No files specified, read from stdin
        character_matrix = read_character_matrix(sys.stdin)
    else:
        for file_path in args.input_files:
            try:
                with open(file_path, 'r') as f:
                    # If multiple files, append rows from each file to the same matrix
                    character_matrix.extend(read_character_matrix(f))
            except FileNotFoundError:
                print(f"Error: Input file not found at {file_path}", file=sys.stderr)
                sys.exit(1)
            except Exception as e:
                print(f"An error occurred while processing {file_path}: {e}", file=sys.stderr)
                sys.exit(1)
    
    # Validate Matrix: a. If character_matrix is empty
    if not character_matrix:
        # print("No data provided or all input files were empty.", file=sys.stderr) # Optional
        sys.exit(0) # Exit cleanly, printing nothing.

    # Validate Matrix: b. Determine num_cols and c. Verify consistent line lengths
    num_cols = len(character_matrix[0])
    for i, row_chars in enumerate(character_matrix):
        if len(row_chars) != num_cols:
            print(f"Error: Line {i+1} (0-indexed {i}) has {len(row_chars)} characters, "
                  f"but expected {num_cols} characters (based on the first line). "
                  "Lines must have consistent lengths to transpose character matrix.", file=sys.stderr)
            sys.exit(1)

    # Transpose the Matrix
    transposed_matrix = []
    if num_cols > 0: # Only proceed if there are characters in lines
        for j in range(num_cols): # Iterate original character positions (columns)
            new_row_chars = []
            for i in range(len(character_matrix)): # Iterate original lines (rows)
                new_row_chars.append(character_matrix[i][j])
            transposed_matrix.append(new_row_chars)
    
    # Print Transposed Matrix
    for row_chars in transposed_matrix:
        # The Perl script prints " $LoL[$i][$j]", which means a leading space for each char.
        # " ".join(row_chars) will put a space *between* chars.
        # To match the Perl script's " $char1 $char2 $char3" format for a row:
        # print(" " + " ".join(row_chars))
        # However, the prompt asks for "characters in row joined by a single space"
        # which implies "char1 char2 char3". Let's stick to the prompt's specific instruction for Python.
        print(" ".join(row_chars))


if __name__ == "__main__":
    main()
