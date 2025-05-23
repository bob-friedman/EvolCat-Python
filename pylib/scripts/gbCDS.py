#!/usr/bin/env python3

import argparse
import sys
from pylib.utils import seq_parser # Assuming pylib is in PYTHONPATH or script is run from the project root
from Bio.SeqFeature import FeatureLocation # For type hinting if needed, not strictly necessary for feature.extract

def main():
    """
    Extracts information about the first CDS feature from a GenBank file.
    For each record, prints Locus ID, CDS nucleotide sequence, and CDS translation.
    """
    parser = argparse.ArgumentParser(
        description="Extract information (Locus ID, nucleotide sequence, translation) "
                    "for the first CDS feature in each record of a GenBank file."
    )
    parser.add_argument(
        "input_genbank_file",
        help="Path to the input GenBank file."
    )
    args = parser.parse_args()

    try:
        for record in seq_parser.parse_genbank_file(args.input_genbank_file):
            found_cds = False
            for feature in record.features:
                if feature.type == "CDS":
                    # Extract nucleotide sequence of this CDS feature
                    # feature.extract() correctly handles complex locations (joins, etc.)
                    cds_nucleotide_sequence = feature.extract(record.seq)
                    
                    # Get protein translation
                    # The qualifier value is a list, take the first element.
                    # Default to "N/A" if "translation" qualifier is missing or empty.
                    translation = feature.qualifiers.get("translation", ["N/A"])[0]
                    if not translation.strip(): # Handle empty string in list
                        translation = "N/A"
                        
                    print(record.id)
                    print(str(cds_nucleotide_sequence))
                    print(translation)
                    
                    found_cds = True
                    break # Process only the first CDS feature

            if not found_cds:
                print(record.id)
                print("No CDS feature found")
                print("N/A") # Maintain consistent output structure

    except FileNotFoundError:
        print(f"Error: Input file not found at {args.input_genbank_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e: # Broad exception for parsing errors or other issues
        print(f"An error occurred while processing {args.input_genbank_file}: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
