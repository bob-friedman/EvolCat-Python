#!/usr/bin/env python3

import argparse
import sys

IUPAC_TO_REGEX = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': '[GA]', 'Y': '[CT]', 'M': '[AC]', 'K': '[GT]',
    'S': '[GC]', 'W': '[AT]', 'B': '[CGT]',
    'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]',
    'N': '[ACGT]', # As per Perl script, N maps to [ACGT]
    # U is sometimes used for T, but not in the original Perl script's map
}

def main():
    """
    Converts IUPAC strings (one per line from input file) to regular expressions.
    """
    parser = argparse.ArgumentParser(
        description="Convert IUPAC ambiguity strings to their equivalent regular expressions."
    )
    parser.add_argument(
        "input_file",
        help="Path to the input file containing IUPAC strings, one per line."
    )
    args = parser.parse_args()

    unknown_char_warning_issued = False

    try:
        with open(args.input_file, 'r') as f:
            for line in f:
                iupac_string = line.strip()
                if not iupac_string: # Skip empty lines
                    continue
                
                # Remove all '^' characters
                iupac_string_cleaned = iupac_string.replace('^', '')
                
                regex_output = [] # Use a list for efficiency, then join
                
                for char_code in iupac_string_cleaned:
                    char_upper = char_code.upper() # Ensure lookup is case-insensitive for IUPAC codes
                    regex_component = IUPAC_TO_REGEX.get(char_upper)
                    
                    if regex_component:
                        regex_output.append(regex_component)
                    else:
                        # Character not in IUPAC_TO_REGEX map
                        regex_output.append(char_code) # Append the original character as a literal
                        if not unknown_char_warning_issued:
                            print(f"Warning: Unknown IUPAC code(s) encountered (e.g., '{char_code}'). "
                                  "Unknown characters are being appended as literals to the regex.", file=sys.stderr)
                            unknown_char_warning_issued = True
                
                print("".join(regex_output))

    except FileNotFoundError:
        print(f"Error: Input file not found at {args.input_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
