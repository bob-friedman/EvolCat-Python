#!/usr/bin/env python3

import argparse

def main():
    """
    Joins two tab-delimited files based on the first column.
    Mimics the behavior of the Perl script tables/join.pl.
    """
    parser = argparse.ArgumentParser(
        description="Join two tab-delimited files based on the first column."
    )
    parser.add_argument("file1", help="Path to the first input file (left file for the join).")
    parser.add_argument("file2", help="Path to the second input file (right file for the join).")
    args = parser.parse_args()

    data_file2 = {}

    try:
        with open(args.file2, 'r') as f2:
            for line in f2:
                stripped_line = line.rstrip('\n')
                parts = stripped_line.split('\t', 1) # Split only on the first tab
                if parts: # Ensure the line is not empty
                    key = parts[0]
                    data_file2[key] = stripped_line
    except FileNotFoundError:
        print(f"Error: File not found at {args.file2}")
        return
    except Exception as e:
        print(f"An error occurred while reading {args.file2}: {e}")
        return

    try:
        with open(args.file1, 'r') as f1:
            for line1 in f1:
                stripped_line1 = line1.rstrip('\n')
                key1_parts = stripped_line1.split('\t', 1) # Split only on the first tab
                
                if not key1_parts: # Handle empty lines in file1
                    print(f"{stripped_line1}\t")
                    continue

                key1 = key1_parts[0]
                
                # Mimic Perl's behavior: if key is not found, $hash{$key} is undef,
                # which prints as an empty string in a print statement.
                line_from_file2 = data_file2.get(key1, "") 
                
                print(f"{stripped_line1}\t{line_from_file2}")

    except FileNotFoundError:
        print(f"Error: File not found at {args.file1}")
        return
    except Exception as e:
        print(f"An error occurred while reading {args.file1} or processing data: {e}")
        return

if __name__ == "__main__":
    main()
