#!/usr/bin/env python3

import argparse
import sys

def main():
    """
    Finds the best approximate match of a pattern in a text using edit distance.
    """
    parser = argparse.ArgumentParser(
        description="Find the best approximate match of a pattern in a text using edit distance."
    )
    parser.add_argument(
        "--pattern",
        required=True,
        help="The pattern string."
    )
    parser.add_argument(
        "--text",
        required=True,
        help="The text string in which to search for the pattern."
    )
    parser.add_argument(
        "--print_matrix",
        action='store_true',
        help="Print the full distance matrix."
    )
    args = parser.parse_args()

    pattern = args.pattern
    text = args.text

    PLEN = len(pattern)
    TLEN = len(text)

    # Initialize the distance matrix D
    # D will be (TLEN + 1) rows x (PLEN + 1) columns
    D = [[0 for _ in range(PLEN + 1)] for _ in range(TLEN + 1)]

    # Set initial conditions
    # D[t][0] = 0 for finding pattern anywhere in text (cost of deletions from text is 0 at start)
    for t in range(TLEN + 1):
        D[t][0] = 0
    
    # D[0][p] = p (cost of insertions to match pattern prefix if text is empty)
    for p in range(PLEN + 1):
        D[0][p] = p

    # Fill the rest of the matrix
    for t in range(1, TLEN + 1):
        for p in range(1, PLEN + 1):
            cost = 0 if text[t-1] == pattern[p-1] else 1
            D[t][p] = min(
                D[t-1][p-1] + cost,  # Match/Mismatch
                D[t-1][p] + 1,       # Deletion from text
                D[t][p-1] + 1        # Insertion into text
            )

    # Print the distance matrix if requested
    if args.print_matrix:
        print("Distance Matrix (D):")
        # Header row for pattern
        header = "    " + " ".join(f"{pat_char:2}" for pat_char in pattern)
        print(header)
        print("   " + "-" * (len(header)-3))

        for t in range(TLEN + 1):
            row_str = f"{text[t-1] if t > 0 else ' ':1} |"
            for p in range(PLEN + 1):
                row_str += f" {D[t][p]:2}"
            print(row_str)
        print("-" * (len(header) + 1))


    # Find the best match(es)
    best_score = float('inf')
    matches_end_positions = [] # Will store 1-based end positions in text

    # Iterate through the last column of D (where pattern is fully consumed)
    # t represents the end position in text (1-based due to loop range and matrix structure)
    for t in range(1, TLEN + 1): 
        score = D[t][PLEN]
        if score < best_score:
            best_score = score
            matches_end_positions = [t]
        elif score == best_score:
            matches_end_positions.append(t)
    
    # If PLEN is 0 (empty pattern), best_score would be 0 and all positions are matches.
    # The Perl script likely assumes non-empty pattern.
    # If pattern is empty, D[t][0] = 0, so best_score is 0, matches_end_positions are all t.
    # If pattern is empty (PLEN=0), D[t][0] is 0. The loop for p won't run for D[0][p]=p.
    # D[t][PLEN] where PLEN=0 means D[t][0], which is always 0.
    # So best_score will be 0, and matches_end_positions will be 1..TLEN.
    # This seems like a reasonable interpretation for an empty pattern.

    # Print results
    print(f"\nPattern: {pattern}")
    print(f"Text:    {text}")
    print(f"Best edit distance: {best_score}")
    if matches_end_positions:
        # Convert 0-based list of end positions to 1-based for reporting if necessary.
        # However, since our 't' loop for finding best match goes from 1 to TLEN,
        # and D[t] corresponds to text ending at text[t-1], 't' is already 1-based.
        print(f"Ending position(s) in text (1-based): {', '.join(map(str, matches_end_positions))}")
    else:
        # This case should ideally not happen if TLEN > 0, as there will always be a score.
        # Could happen if TLEN is 0 and PLEN > 0, then best_score remains float('inf').
        if TLEN == 0 and PLEN > 0:
            print("Ending position(s) in text: N/A (empty text, non-empty pattern)")
            # best_score would be D[0][PLEN] which is PLEN.
            # The loop `for t in range(1, TLEN + 1)` won't run if TLEN=0
            # So best_score remains float('inf') and matches_end_positions is empty.
            # Let's correct this part for empty text:
            if TLEN == 0 and PLEN > 0:
                best_score = D[0][PLEN] # Should be PLEN
                print(f"Best edit distance: {best_score}") # Re-print corrected score
                print("Ending position(s) in text: N/A (empty text)")

if __name__ == "__main__":
    main()
