#!/usr/bin/env python3

import argparse
import sys
import collections
from pylib.utils import seq_parser

def main():
    parser = argparse.ArgumentParser(
        description="Find tandem repeats of a pattern in a DNA sequence."
    )
    parser.add_argument(
        "--sequence",
        help="The DNA sequence string."
    )
    parser.add_argument(
        "--sequence_file",
        help="Path to a file containing the DNA sequence. If FASTA, first record is used. Otherwise, first non-empty line."
    )
    parser.add_argument(
        "--pattern",
        required=True,
        help="The repeat pattern string."
    )
    parser.add_argument(
        "--fudge",
        type=int,
        default=1,
        help="Maximum allowed separation (fudge factor) between repeats (default: 1)."
    )
    parser.add_argument(
        "--verbose",
        action='store_true',
        help="Print messages during discovery."
    )
    args = parser.parse_args()

    sequence = ""
    if args.sequence_file:
        try:
            # Attempt to parse as FASTA first
            try:
                records = list(seq_parser.parse_fasta_file(args.sequence_file))
                if records:
                    sequence = str(records[0].seq)
                else:
                    # If not FASTA or empty FASTA, try reading as plain text
                    with open(args.sequence_file, 'r') as f:
                        for line in f:
                            stripped_line = line.strip()
                            if stripped_line:
                                sequence = stripped_line
                                break # Use first non-empty line
            except Exception: # Broad exception if FASTA parsing fails
                 with open(args.sequence_file, 'r') as f:
                    for line in f:
                        stripped_line = line.strip()
                        if stripped_line:
                            sequence = stripped_line
                            break 
            if not sequence:
                print(f"Error: Could not read a sequence from file: {args.sequence_file}", file=sys.stderr)
                sys.exit(1)
        except FileNotFoundError:
            print(f"Error: Sequence file not found: {args.sequence_file}", file=sys.stderr)
            sys.exit(1)
    elif args.sequence:
        sequence = args.sequence
    else:
        print("Error: Either --sequence or --sequence_file must be provided.", file=sys.stderr)
        sys.exit(1)

    if not args.pattern:
        print("Error: --pattern must be provided.", file=sys.stderr) # Should be caught by argparse required=True
        sys.exit(1)
    
    pattern = args.pattern
    pattern_len = len(pattern)
    if pattern_len == 0:
        print("Error: Pattern cannot be empty.", file=sys.stderr)
        sys.exit(1)

    sequence = sequence.upper()
    pattern = pattern.upper()

    results_store = collections.defaultdict(list)
    
    search_start_idx = 0
    current_block_start_offset = -1
    repeats_in_current_block = 0
    last_known_match_start_offset = -1 

    while search_start_idx < len(sequence):
        match_start_offset = sequence.find(pattern, search_start_idx)

        if match_start_offset != -1:
            if current_block_start_offset == -1: # First match of a new block
                current_block_start_offset = match_start_offset
                repeats_in_current_block = 1
                if args.verbose:
                    print(f"Found at {match_start_offset}, repeats {repeats_in_current_block}")
            # Condition for starting a new block vs continuing existing one
            # A new block starts if this match is too far from the *start* of the previous match
            elif match_start_offset > last_known_match_start_offset + pattern_len + args.fudge:
                # Record previous block
                if repeats_in_current_block > 0: # Ensure there was a block to record
                    block_end_to_store = last_known_match_start_offset + pattern_len - 1
                    results_store[repeats_in_current_block].append((current_block_start_offset, block_end_to_store))
                
                # Reset for new block
                current_block_start_offset = match_start_offset
                repeats_in_current_block = 1
                if args.verbose:
                    print(f"Found at {match_start_offset}, repeats {repeats_in_current_block}")
            else: # Continues current block
                repeats_in_current_block += 1
                if args.verbose:
                    print(f"Found at {match_start_offset}, repeats {repeats_in_current_block}")
            
            last_known_match_start_offset = match_start_offset
            search_start_idx = match_start_offset + 1 # Start next search after the beginning of current match
        else: # No more matches
            break
            
    # After the loop, record the last block if it exists
    if repeats_in_current_block > 0:
        block_end_to_store = last_known_match_start_offset + pattern_len - 1
        results_store[repeats_in_current_block].append((current_block_start_offset, block_end_to_store))

    # Print Summary Output
    if not results_store:
        print("No repeats found.")
    else:
        for repeat_count in sorted(results_store.keys()):
            blocks = results_store[repeat_count]
            plural_s = "s" if repeat_count > 1 else ""
            found_plural_s = "s" if len(blocks) > 1 else "" # For "1 found" vs "2 founds"
            print(f"{repeat_count} repeat{plural_s}")
            # Correcting "founds"
            print(f"\t{len(blocks)} instance{'' if len(blocks) == 1 else 's'} found")


            for block_start, block_end in blocks:
                block_len = block_end - block_start + 1
                print(f"\t\tOffset {block_start} - {block_end} ({block_len})")

if __name__ == "__main__":
    main()
