#!/usr/bin/env python3

"""
Processes GenBank files to extract and analyze coding sequence (CDS) information.

This script provides two main functionalities:
1. Range Extraction Mode: Extracts nucleotide sequences from a specified
   genomic range. If this range overlaps with CDS features, it translates
   the overlapping CDS portions into protein sequences.
2. Single Position Mode: Identifies the codon and corresponding amino acid
   at a specific nucleotide position if it falls within a CDS feature.

The script uses BioPython for GenBank parsing and sequence operations.
It handles various CDS feature complexities, including forward/reverse strands,
compound locations (exons), and `codon_start` qualifiers.

Usage via command line:
  For range extraction:
    python extract_cds_region.py <genbank_file> --start_pos <start> --end_pos <end>
  For single position analysis:
    python extract_cds_region.py <genbank_file> --single_position <position>

See pylib/scripts/README.md for detailed documentation.
"""

import argparse
import sys
import os

# Add the project root directory to the Python path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.insert(0, project_root)

from pylib.utils import seq_parser
from Bio.Seq import Seq # Removed UnknownSeq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, CompoundLocation, ExactPosition
from Bio.Data import CodonTable

def extract_and_translate_cds_region(genbank_file, start_pos, end_pos):
    """
    Extracts a nucleotide sequence from a GenBank record based on 1-based start and end
    positions, and translates it to a protein sequence if it overlaps with any CDS feature.

    This function iterates through records in the provided GenBank file. For each record,
    it extracts the nucleotide sequence corresponding to the given `start_pos` and `end_pos`.
    It then identifies all CDS features within that record and checks for overlaps
    with the queried genomic range.

    For each overlapping CDS:
    - It determines the exact nucleotide sequence from the CDS that falls within the query.
    - It considers the CDS's strand, `codon_start` qualifier, and `transl_table`.
    - It translates this portion into an amino acid sequence.
    - Results, including nucleotide and protein sequences, are printed to standard output.

    Args:
        genbank_file (str): Path to the input GenBank file.
        start_pos (int): The 1-based start position of the genomic region to query.
        end_pos (int): The 1-based end position of the genomic region to query.

    Prints:
        Detailed information about the processed records, extracted sequences,
        CDS features found, and translated protein sequences to standard output.
        Error messages are printed to standard error if issues occur (e.g.,
        positions out of bounds, file not found).

    Returns:
        None
    """
    if start_pos > end_pos:
        print(f"Error: Start position ({start_pos}) cannot be greater than end position ({end_pos}).", file=sys.stderr)
        return

    if start_pos <= 0 or end_pos <= 0:
        print(f"Error: Start and end positions must be positive integers.", file=sys.stderr)
        return

    # User inputs are 1-based, convert to 0-based for internal use
    query_start_0based = start_pos - 1
    query_end_0based = end_pos # Biopython slicing [start:end] means end is exclusive

    try:
        for record in seq_parser.parse_genbank_file(genbank_file):
            print(f"Processing record: {record.id} (length {len(record.seq)})")

            if not (0 <= query_start_0based < len(record.seq) and query_end_0based <= len(record.seq)):
                print(f"Error: Query range {start_pos}-{end_pos} is out of bounds for record {record.id} (length {len(record.seq)}). Skipping record.", file=sys.stderr)
                continue

            # Extract the raw nucleotide sequence from the record based on the user's query range.
            queried_nucleotide_sequence_full_range = record.seq[query_start_0based:query_end_0based]
            print(f"Nucleotide sequence in query range ({start_pos}-{end_pos}): {queried_nucleotide_sequence_full_range}")

            any_cds_produced_sequence_in_record = False # Renamed and will be used correctly

            # Iterate through all features in the record to find CDS features.
            for feature in record.features:
                if feature.type == "CDS":
                    cds_id = feature.qualifiers.get('protein_id', feature.qualifiers.get('locus_tag', ['Unknown_CDS']))[0]
                    print(f"  Found CDS: {cds_id} located at {feature.location}")

                    # Extract the full coding sequence for this CDS.
                    # BioPython's feature.extract() correctly handles strand orientation (forward/reverse)
                    # and compound locations (e.g., joins for exons), returning the 5'-3' coding sequence.
                    full_cds_sequence = feature.extract(record.seq)
                    # Removed UnknownSeq check: if isinstance(full_cds_sequence, UnknownSeq):
                    # The Seq object itself can handle unknown characters.
                    # A general check for 'N's or other non-standard characters can be added if needed,
                    # but Biopython's translate handles them by typically producing 'X'.
                    if "N" in str(full_cds_sequence).upper(): # Check for unknown bases
                        print(f"    Warning: CDS {cds_id} sequence contains 'N' characters. Translation may produce 'X' or fail depending on table.")
                    if not full_cds_sequence:
                        print(f"    Warning: CDS {cds_id} feature extracted an empty sequence. Skipping.")
                        continue

                    target_sub_sequence = Seq("") # Initialize as an empty Bio.Seq object

                    # Create a map from genomic 0-based indices to 0-based indices within the full_cds_sequence.
                    # This map is crucial for correctly piecing together the CDS portion that overlaps the user's query,
                    # especially for CDS with compound locations (exons) or on the reverse strand.
                    genomic_to_cds_index_map = {}
                    map_idx_counter = 0 # Counter for position in the full_cds_sequence
                    if feature.location.strand == 1: # Forward strand
                        for part in feature.location.parts:
                            # For each genomic position in this part, map it to its corresponding index in full_cds_sequence
                            for i in range(int(part.start), int(part.end)):
                                genomic_to_cds_index_map[i] = map_idx_counter
                                map_idx_counter += 1
                    else: # Reverse strand
                        # For reverse strand, iterate parts and bases in transcription order (genomic high to low for each part, parts high to low)
                        for part in reversed(feature.location.parts):
                            for i in reversed(range(int(part.start), int(part.end))):
                                genomic_to_cds_index_map[i] = map_idx_counter
                                map_idx_counter += 1

                    # List to store the 0-based indices from full_cds_sequence that fall within the user's genomic query range.
                    indices_in_full_cds_for_query = []
                    for record_idx_in_query in range(query_start_0based, query_end_0based):
                        if record_idx_in_query in genomic_to_cds_index_map:
                            indices_in_full_cds_for_query.append(genomic_to_cds_index_map[record_idx_in_query])

                    if indices_in_full_cds_for_query:
                        # Sort and get unique indices to reconstruct the CDS portion correctly.
                        # Sorting is vital if the query spans multiple exons that are not contiguous in the genome
                        # but are contiguous in the coding sequence. Set ensures uniqueness.
                        sorted_indices_in_full_cds = sorted(list(set(indices_in_full_cds_for_query)))

                        temp_target_seq_list = []
                        for cds_idx in sorted_indices_in_full_cds:
                            try:
                               temp_target_seq_list.append(str(full_cds_sequence[cds_idx]))
                            except IndexError:
                                print(f"    Warning: Index {cds_idx} out of bounds for full_cds_sequence (len {len(full_cds_sequence)}) for CDS {cds_id}. This might indicate an issue with location parsing or mapping.")
                                continue

                        if temp_target_seq_list:
                            target_sub_sequence = Seq("".join(temp_target_seq_list))
                            print(f"    Nucleotide sequence from CDS ({cds_id}) overlapping query: {target_sub_sequence}")
                            # found_cds_in_query_region = True # This CDS produced sequence
                        # Removed the else block that set found_cds_in_query_region to False here
                    else: # No indices from the query range map to the CDS for this feature
                        # This print is fine, it's specific to this CDS feature
                        print(f"    CDS {cds_id} does not have exonic bases within the query range {start_pos}-{end_pos}.")
                        # Deliberately not changing any_cds_produced_sequence_in_record here

                    if target_sub_sequence: # If this particular CDS feature yielded a sequence
                        any_cds_produced_sequence_in_record = True # Mark that at least one CDS in record was processed
                        if not indices_in_full_cds_for_query: # Should ideally not happen if target_sub_sequence is true
                             print(f"    Internal logic error: target_sub_sequence is present but its source indices (indices_in_full_cds_for_query) are not. Skipping translation for {cds_id}.")
                             continue

                        # The 0-based index of the first base of the overlapping CDS region within the *full* CDS sequence.
                        # This is used to calculate the correct frame for translation.
                        first_base_offset_in_cds = sorted_indices_in_full_cds[0]

                        # Get the 'codon_start' qualifier (1, 2, or 3), defaulting to 1 if not present.
                        # This indicates the reading frame of the CDS.
                        codon_start_qualifier = int(feature.qualifiers.get("codon_start", [1])[0])

                        # Calculate the effective frame for the *extracted segment* (target_sub_sequence).
                        # (codon_start_qualifier - 1) is the 0-based frame of the full CDS.
                        # first_base_offset_in_cds is how many bases into the full CDS our extracted segment begins.
                        # The sum modulo 3 gives the 0-based starting frame for our specific segment.
                        effective_frame_start = ( (codon_start_qualifier - 1) + first_base_offset_in_cds ) % 3

                        # Slice the extracted CDS segment to begin at the correct frame for translation.
                        sequence_for_translation = target_sub_sequence[effective_frame_start:]

                        table_id = int(feature.qualifiers.get("transl_table", [1])[0])

                        if len(sequence_for_translation) < 3:
                            protein_sequence = "Error: Not enough nucleotides for a codon after frame adjustment."
                            print(f"    {protein_sequence} (CDS {cds_id}, length {len(sequence_for_translation)}, frame adj {effective_frame_start})")
                        else:
                            try:
                                protein_sequence = sequence_for_translation.translate(table=table_id, to_stop=True, cds=False)
                                print(f"    Protein sequence from CDS {cds_id} (frame adj: {effective_frame_start}, table: {table_id}): {protein_sequence}")
                            except CodonTable.TranslationError as te:
                                protein_sequence = f"Error during translation: {te}"
                                print(f"    {protein_sequence} (CDS {cds_id})")
                            except Exception as e:
                                protein_sequence = f"Unexpected error during translation: {e}"
                                print(f"    {protein_sequence} (CDS {cds_id})")
                    # Removed: elif found_cds_in_query_region and not target_sub_sequence :
                    # This condition is covered because target_sub_sequence would be empty.
                    # A message for a specific CDS not yielding sequence is already printed if indices_in_full_cds_for_query is empty.


            if not any_cds_produced_sequence_in_record: # Check after iterating all features in the record
                print(f"  No CDS features found whose exonic parts overlap with the query range {start_pos}-{end_pos}.")

    except FileNotFoundError:
        print(f"Error: GenBank file not found at {genbank_file}", file=sys.stderr)
        sys.exit(1)
    except CodonTable.NCBICodonTableNotFoundError as e: # Ensure this is the correct exception type
        print(f"Error: Invalid NCBI translation table ID used in a feature: {e}", file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1) # Exit if an unexpected error occurs during range processing

# New function placeholder
def process_single_position_mode(genbank_file, target_nucleotide_pos_1based):
    """
    Identifies and analyzes the codon containing a single specified nucleotide position
    within CDS features of a GenBank file.

    This function iterates through records in the GenBank file. For each record, it
    checks if the `target_nucleotide_pos_1based` falls within any CDS feature.

    If the position is found within a CDS:
    - It determines the specific codon containing this nucleotide.
    - It identifies the position of the target nucleotide within that codon (1st, 2nd, or 3rd).
    - It translates the codon to its corresponding amino acid using the `transl_table`
      and `codon_start` qualifier from the CDS feature.
    - Detailed information, including record ID, CDS ID, genomic location,
      CDS nucleotide index, codon sequence, and translated amino acid, is printed
      in a formatted block to standard output.

    Args:
        genbank_file (str): Path to the input GenBank file.
        target_nucleotide_pos_1based (int): The 1-based genomic nucleotide position to analyze.

    Prints:
        Formatted information about the codon and amino acid at the target position
        if it's found within a CDS.
        Prints messages if the position is not within a CDS, is before the translation
        start, or is part of an incomplete codon.
        Error messages are printed to standard error for issues like file not found
        or invalid positions.

    Returns:
        None
    """
    if target_nucleotide_pos_1based <= 0:
        print(f"Error: Target nucleotide position must be positive. Got {target_nucleotide_pos_1based}", file=sys.stderr)
        sys.exit(1)

    # User input is 1-based, convert to 0-based for internal calculations.
    target_pos_0based = target_nucleotide_pos_1based - 1
    found_target_in_any_cds = False

    try:
        for record in seq_parser.parse_genbank_file(genbank_file):
            print(f"Processing record: {record.id} (length {len(record.seq)}) for position {target_nucleotide_pos_1based}")

            if not (0 <= target_pos_0based < len(record.seq)):
                print(f"  Error: Target position {target_nucleotide_pos_1based} is out of bounds for record {record.id} (length {len(record.seq)}). Skipping record.", file=sys.stderr)
                continue

            for feature in record.features:
                if feature.type == "CDS":
                    cds_id = feature.qualifiers.get('protein_id', feature.qualifiers.get('locus_tag', ['Unknown_CDS']))[0]
                    # print(f"  Checking CDS: {cds_id} located at {feature.location}") # Debug

                    # This will store the 0-based index of the target nucleotide within the full CDS sequence, if found.
                    pos_in_cds_0based = None
                    # Counter for the current position in the conceptual full_cds_sequence as we iterate through parts.
                    current_map_idx = 0

                    # Determine iteration order for parts based on strand
                    parts_to_iterate = feature.location.parts
                    if feature.location.strand == -1:
                        # For reverse strand, BioPython's feature.extract() effectively processes parts
                        # as if they were ordered from "highest coordinate part" to "lowest coordinate part"
                        # to construct the 5'-3' coding sequence.
                        # So, when we map genomic coordinates to CDS coordinates, we iterate parts in this order.
                        parts_to_iterate = reversed(feature.location.parts)

                    for part in parts_to_iterate:
                        part_start_0based = int(part.start)
                        part_end_0based = int(part.end) # Exclusive end for length

                        if target_pos_0based >= part_start_0based and target_pos_0based < part_end_0based:
                            # Target nucleotide is within this part (e.g., exon) of the CDS feature.
                            if feature.location.strand == 1: # Forward strand
                                # Offset within this specific part.
                                offset_in_part = target_pos_0based - part_start_0based
                                # Position in full_cds_sequence is sum of lengths of prior parts plus offset in current part.
                                pos_in_cds_0based = current_map_idx + offset_in_part
                            else: # Reverse strand
                                # For reverse strand, the 5' end of this exon (in coding sense) is part_end_0based-1 on the genome.
                                # Offset is calculated from this 5' end.
                                offset_in_part = (part_end_0based - 1) - target_pos_0based
                                # Position in full_cds_sequence is sum of lengths of prior parts (iterated in reverse genomic order) plus offset in current part.
                                pos_in_cds_0based = current_map_idx + offset_in_part
                            break # Found the part containing the target nucleotide

                        current_map_idx += len(part) # Add length of this part to map_idx for next iteration.

                    if pos_in_cds_0based is not None:
                        found_target_in_any_cds = True
                        # cds_id is already defined from: feature.qualifiers.get('protein_id', feature.qualifiers.get('locus_tag', ['Unknown_CDS']))[0]
                        print(f"  Target nucleotide {target_nucleotide_pos_1based} (0-based genomic: {target_pos_0based}) found in CDS: {cds_id}")
                        print(f"    Original CDS location: {feature.location}")
                        print(f"    Target is at 0-based index {pos_in_cds_0based} within the conceptual full CDS sequence.")

                        full_cds_sequence = feature.extract(record.seq)
                        if not isinstance(full_cds_sequence, Seq): # Should be a Seq object
                            full_cds_sequence = Seq(str(full_cds_sequence))

                        if pos_in_cds_0based >= len(full_cds_sequence):
                            print(f"    Error: Calculated pos_in_cds_0based ({pos_in_cds_0based}) is out of bounds for extracted full_cds_sequence (length {len(full_cds_sequence)}). This might indicate an issue with location parsing or mapping. Skipping this CDS.")
                            # Reset found_target_in_any_cds if this was the only potential find and it's problematic
                            # This depends on whether we want to report the CDS if we can't get the nucleotide
                            # For now, let's consider it a failed attempt for this CDS.
                            # If we 'continue', found_target_in_any_cds should not be reset here, rather it should not have been set true for THIS cds yet.
                            # Let's adjust: only set found_target_in_any_cds = True AFTER all checks pass for this CDS.
                            # However, the current logic sets it if pos_in_cds_0based is not None.
                            # For this step, let's keep found_target_in_any_cds as True and print error.
                            continue # To the next feature

                        # The actual nucleotide base at the calculated 0-based index within the full CDS sequence.
                        target_nucleotide_in_cds = full_cds_sequence[pos_in_cds_0based]
                        print(f"    Nucleotide at this CDS index: {target_nucleotide_in_cds}")

                        codon_start_qualifier = int(feature.qualifiers.get("codon_start", [1])[0])
                        # The 0-based index on the full_cds_sequence where translation actually begins, accounting for codon_start.
                        cds_translation_start_offset = codon_start_qualifier - 1

                        # The 0-based index of the target nucleotide relative to the actual start of translation within the CDS.
                        # This is used to determine its position within a codon.
                        target_pos_relative_to_translation_start = pos_in_cds_0based - cds_translation_start_offset

                        if target_pos_relative_to_translation_start < 0:
                            print(f"    Position {target_nucleotide_pos_1based} is in CDS {cds_id} (at CDS index {pos_in_cds_0based}) but occurs *before* the translation start indicated by codon_start={codon_start_qualifier}. It is not part of a translated codon.")
                            break # Stop processing this specific CDS here, move to next record or finish.

                        # Calculate the 0-based start index of the codon that contains the target nucleotide.
                        # This index is relative to the start of the translated portion of the CDS.
                        codon_internal_start_0based = (target_pos_relative_to_translation_start // 3) * 3

                        # Calculate the actual 0-based start index of this codon within the *full_cds_sequence*.
                        codon_start_in_full_cds = cds_translation_start_offset + codon_internal_start_0based

                        # Extract the three-nucleotide codon sequence from the full CDS sequence.
                        the_codon_seq = full_cds_sequence[codon_start_in_full_cds : codon_start_in_full_cds + 3]

                        if len(the_codon_seq) == 3:
                            # Determine the 1-based position (1, 2, or 3) of the target nucleotide within its codon.
                            base_in_codon_1based = (target_pos_relative_to_translation_start % 3) + 1

                            # Translate the codon
                            table_id = int(feature.qualifiers.get("transl_table", [1])[0])
                            protein_letter_str = "" # Initialize to empty string
                            try:
                                # Ensure the_codon_seq is a Seq object for translation
                                if not isinstance(the_codon_seq, Seq):
                                    codon_to_translate = Seq(str(the_codon_seq))
                                else:
                                    codon_to_translate = the_codon_seq

                                # Check for 'N' or other non-DNA characters in codon before translating
                                if any(c not in "ATGCatgc" for c in str(codon_to_translate)):
                                    protein_letter_str = "X (contains non-standard base)"
                                else:
                                    protein_letter = codon_to_translate.translate(table=table_id)
                                    protein_letter_str = str(protein_letter)
                                    # Standard tables return '*' for stop. No need to manually strip if that's desired.

                            except CodonTable.TranslationError as te:
                                protein_letter_str = f"TranslationError ({te})"
                            except Exception as e: # Catch any other unexpected error during translation
                                protein_letter_str = f"Error translating ({e})"

                            # Output the collected information
                            print(f"  ---------------------------------------------------")
                            print(f"  Record ID:                    {record.id}")
                            print(f"  CDS ID:                       {cds_id}")
                            print(f"  CDS Location:                 {feature.location}")
                            print(f"  Target Nucleotide Position:   {target_nucleotide_pos_1based} (genomic, 1-based)")
                            print(f"  CDS Nucleotide Index:         {pos_in_cds_0based} (0-based, within full CDS sequence)")
                            print(f"  Target Base in Codon:       {base_in_codon_1based}")
                            print(f"  Codon Sequence:               {the_codon_seq}")
                            print(f"  Translation Table ID:         {table_id}")
                            print(f"  Translated Codon:             {protein_letter_str}")
                            print(f"  ---------------------------------------------------")

                        else: # Incomplete codon
                            print(f"  ---------------------------------------------------")
                            print(f"  Record ID:                    {record.id}")
                            print(f"  CDS ID:                       {cds_id}")
                            print(f"  CDS Location:                 {feature.location}")
                            print(f"  Target Nucleotide Position:   {target_nucleotide_pos_1based} (genomic, 1-based)")
                            print(f"  CDS Nucleotide Index:         {pos_in_cds_0based} (0-based, within full CDS sequence)")
                            print(f"  Codon Sequence:               Fragment '{the_codon_seq}' (length {len(the_codon_seq)})")
                            print(f"  Note:                         Target nucleotide is part of an incomplete codon at the end of CDS {cds_id}.")
                            print(f"  ---------------------------------------------------")

                        break # Stop checking other features in this record once a relevant one is processed.

            if found_target_in_any_cds: # If found in current record
                 # break # Uncomment this to stop processing further records after first hit.
                 # For now, it will process all records and report first hit in each.
                 pass

        if not found_target_in_any_cds: # This message will print if the loop finishes for all records and target was never found.
                                      # This needs to be inside the record loop if we want a per-record "not found" message,
                                      # or outside if it's an overall "not found in any record".
                                      # The current logic will print it once if never found in *any* record.
                                      # Let's adjust for per-record status.
            # This 'if not found_target_in_any_cds' would be better inside the record loop for per-record message.
            # However, the 'found_target_in_any_cds' flag as defined is global across records.
            # For this step, let's keep it as is; refinement of per-record/global messages can come later.
            pass # The final "not found" message below handles the global case.

    except FileNotFoundError:
        print(f"Error: GenBank file not found at {genbank_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

    if not found_target_in_any_cds: # Final check after all records.
        print(f"Target nucleotide position {target_nucleotide_pos_1based} was not found within any CDS feature in the processed GenBank file.")

def main():
    """
    Main function to parse command-line arguments and call the appropriate processing function.

    It sets up `argparse` to handle command-line inputs, distinguishing between
    range extraction mode (using `--start_pos` and `--end_pos`) and single
    position mode (using `--single_position`). Based on the provided arguments,
    it validates them and then invokes either `extract_and_translate_cds_region`
    or `process_single_position_mode`.
    """
    parser = argparse.ArgumentParser(
        description="Processes GenBank records to either extract and translate a nucleotide range "
                    "overlapping CDS features, or identify and translate the specific codon "
                    "for a single nucleotide position."
    )
    parser.add_argument(
        "genbank_file",
        help="Path to the input GenBank file."
    )

    parser.add_argument(
        "--start_pos",
        type=int,
        help="Start base pair of the region to extract (1-based). Required for range "
             "extraction mode; ignored if --single_position is used."
    )
    parser.add_argument(
        "--end_pos",
        type=int,
        help="End base pair of the region to extract (1-based). Required for range "
             "extraction mode; ignored if --single_position is used."
    )
    parser.add_argument(
        "--single_position",
        type=int,
        help="Specify a single 1-based nucleotide position to find its codon. "
             "Activates single-position mode; start_pos and end_pos will be ignored if also provided."
    )

    args = parser.parse_args()

    if args.single_position is not None:
        # Single position mode
        if args.start_pos is not None or args.end_pos is not None:
            print("Warning: --start_pos and --end_pos are ignored when --single_position is used.", file=sys.stderr)

        if args.single_position <= 0:
            parser.error(f"--single_position ({args.single_position}) must be a positive integer.")

        process_single_position_mode(args.genbank_file, args.single_position)

    elif args.start_pos is not None and args.end_pos is not None:
        # Range extraction mode
        if args.start_pos <= 0:
            parser.error(f"--start_pos ({args.start_pos}) must be a positive integer for range mode.")
        if args.end_pos <= 0:
            parser.error(f"--end_pos ({args.end_pos}) must be a positive integer for range mode.")
        if args.start_pos > args.end_pos:
            parser.error(f"--start_pos ({args.start_pos}) cannot be greater than --end_pos ({args.end_pos}) for range mode.")

        extract_and_translate_cds_region(args.genbank_file, args.start_pos, args.end_pos)
    else:
        # Neither mode's required arguments were satisfactorily provided
        parser.print_help(sys.stderr)
        # parser.error("You must specify EITHER --single_position OR (BOTH --start_pos AND --end_pos).")
        # Using print and sys.exit(1) for more control over message formatting if needed outside of parser.error
        print("\nError: You must specify EITHER --single_position OR (BOTH --start_pos AND --end_pos).", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
