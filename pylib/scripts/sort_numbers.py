#!/usr/bin/env python3

import argparse
import sys

def process_lines_for_sorting(lines_iterator, numbers_list):
    """
    Reads lines from an iterator, attempts to convert them to floats,
    and adds them to the numbers_list. Prints a warning for non-numeric lines.
    """
    for line in lines_iterator:
        line_content = line.rstrip('\n')
        try:
            number = float(line_content)
            numbers_list.append(number)
        except ValueError:
            print(f"Warning: Could not convert '{line_content}' to a number. Skipping.", file=sys.stderr)

def main():
    """
    Reads numbers from input file(s) (or stdin), sorts them numerically,
    and prints them to standard output.
    """
    parser = argparse.ArgumentParser(
        description="Sort numbers from input file(s) or stdin and print them. "
                    "Non-numeric lines are skipped with a warning."
    )
    parser.add_argument(
        "input_files",
        nargs='*', # Zero or more input files
        help="Paths to the input files containing numbers. If not specified, reads from stdin."
    )
    args = parser.parse_args()

    numbers_to_sort = []

    if not args.input_files: # No files specified, read from stdin
        # If stdin is a TTY and no data is piped, script will wait for input.
        # This matches the behavior of the Perl script `while(<>)`.
        process_lines_for_sorting(sys.stdin, numbers_to_sort)
    else:
        for file_path in args.input_files:
            try:
                with open(file_path, 'r') as f:
                    process_lines_for_sorting(f, numbers_to_sort)
            except FileNotFoundError:
                print(f"Error: Input file not found at {file_path}", file=sys.stderr)
                # The Perl script would likely die. For robustness and user feedback,
                # we exit here. If processing should continue for other files,
                # this sys.exit(1) could be removed and error reported differently.
                sys.exit(1)
            except Exception as e: # Catch other potential errors during file processing
                print(f"An error occurred while processing {file_path}: {e}", file=sys.stderr)
                sys.exit(1)

    # Sort the collected numbers numerically
    numbers_to_sort.sort() # Default sort for numbers is numerical

    # Print each sorted number
    for number in numbers_to_sort:
        # To match Perl's default output for numbers (e.g., no unnecessary ".0"),
        # we can format integers specially.
        if number == int(number):
            print(int(number))
        else:
            print(number)

if __name__ == "__main__":
    main()
