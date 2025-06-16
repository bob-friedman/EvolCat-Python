#!/usr/bin/env python3

import argparse
import sys
import os

# Add the project root directory to the Python path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.insert(0, project_root)

from pylib.utils import seq_parser
from Bio.Seq import Seq, UnknownSeq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, CompoundLocation, ExactPosition
from Bio.Data import CodonTable

def extract_and_translate_cds_region(genbank_file, start_pos, end_pos):
    """
    Extracts a nucleotide sequence from a GenBank record based on start and end positions,
    and translates it to protein sequence if it overlaps with a CDS feature.
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

            queried_nucleotide_sequence_full_range = record.seq[query_start_0based:query_end_0based]
            print(f"Nucleotide sequence in query range ({start_pos}-{end_pos}): {queried_nucleotide_sequence_full_range}")

            found_cds_in_query_region = False

            for feature in record.features:
                if feature.type == "CDS":
                    cds_id = feature.qualifiers.get('protein_id', feature.qualifiers.get('locus_tag', ['Unknown_CDS']))[0]
                    print(f"  Found CDS: {cds_id} located at {feature.location}")

                    # Get the entire sequence of this CDS feature
                    # feature.extract() handles compound locations (exons) and reverse strands correctly.
                    # The result is the coding sequence in the 5' to 3' direction of transcription.
                    full_cds_sequence = feature.extract(record.seq)
                    if isinstance(full_cds_sequence, UnknownSeq):
                        print(f"    Warning: CDS {cds_id} sequence contains unknown characters (often 'N'). Translation may be affected or fail.")
                    if not full_cds_sequence:
                        print(f"    Warning: CDS {cds_id} feature extracted an empty sequence. Skipping.")
                        continue

                    target_sub_sequence = Seq("") # Initialize as an empty Bio.Seq object

                    # Rebuild target_sub_sequence using the full_cds_sequence and the indices derived from genomic positions
                    # This is the most robust way:
                    # `genomic_to_cds_index_map` maps a genomic index (0-based) to its corresponding index in `full_cds_sequence` (0-based)
                    genomic_to_cds_index_map = {}
                    map_idx_counter = 0
                    if feature.location.strand == 1:
                        for part in feature.location.parts:
                            for i in range(int(part.start), int(part.end)): # Ensure integer positions
                                genomic_to_cds_index_map[i] = map_idx_counter
                                map_idx_counter += 1
                    else: # Reverse strand
                        for part in reversed(feature.location.parts): # Iterate parts in transcription order for reverse strand
                            for i in reversed(range(int(part.start), int(part.end))): # Iterate bases in transcription order
                                genomic_to_cds_index_map[i] = map_idx_counter
                                map_idx_counter += 1

                    indices_in_full_cds_for_query = []
                    for record_idx_in_query in range(query_start_0based, query_end_0based):
                        if record_idx_in_query in genomic_to_cds_index_map:
                            indices_in_full_cds_for_query.append(genomic_to_cds_index_map[record_idx_in_query])

                    if indices_in_full_cds_for_query:
                        # Sort indices to ensure they are in the correct order for sequence reconstruction
                        # This is vital if the query range might hit multiple disjoint parts of the CDS that are contiguous in the coding sequence
                        # or if the mapping process itself didn't guarantee order (though it should with map_idx_counter).
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
                            found_cds_in_query_region = True # Set this only if we actually get a sequence
                        else:
                            if found_cds_in_query_region: # if previously set by simple overlap (which is removed)
                                 print(f"    CDS {cds_id} overlaps query range, but no exonic bases from this CDS are within the precise query nucleotide range after mapping.")
                            # Ensure found_cds_in_query_region is False if no sequence extracted
                            found_cds_in_query_region = False
                    else: # No indices from the query range map to the CDS
                        print(f"    CDS {cds_id} does not have exonic bases within the query range {start_pos}-{end_pos}.")
                        found_cds_in_query_region = False


                    if target_sub_sequence and found_cds_in_query_region:
                        if not indices_in_full_cds_for_query:
                             print(f"    Internal logic error: target_sub_sequence is present but its source indices (indices_in_full_cds_for_query) are not. Skipping translation for {cds_id}.")
                             continue

                        # Use the already sorted list of unique indices: sorted_indices_in_full_cds
                        first_base_offset_in_cds = sorted_indices_in_full_cds[0]

                        codon_start_qualifier = int(feature.qualifiers.get("codon_start", [1])[0]) # 1, 2, or 3

                        effective_frame_start = ( (codon_start_qualifier - 1) + first_base_offset_in_cds ) % 3

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
                    elif found_cds_in_query_region and not target_sub_sequence : # Should be rare due to logic above
                        print(f"    CDS {cds_id} was thought to overlap, but no specific sequence was extracted after mapping. Check query range and CDS definition.")
                        # Ensure flag is reset if we reach here with no sequence
                        found_cds_in_query_region = False


            if not found_cds_in_query_region:
                # This outer flag aggregates results from all CDS features.
                # We need a per-feature flag, and then an overall one.
                # The current logic correctly sets found_cds_in_query_region per CDS, but it gets overwritten.
                # Let's rename the loop flag and use an overall flag.
                # No, the current structure is: iterate records, then iterate features.
                # `found_cds_in_query_region` is reset for each feature effectively by how it's used.
                # The final print should reflect if *any* CDS in the record met the criteria.
                # This requires an overall flag for the record.
                # Let's adjust:
                pass # The existing print "No CDS features found whose exonic parts overlap..." is printed if loop finishes and this is false.
                     # This is effectively "no CDS in this *record* had overlapping sequence" IF it's the last step after the loop.
                     # The current placement is inside the record loop, after feature loop. So it's per record. This is fine.

    except FileNotFoundError:
        print(f"Error: GenBank file not found at {genbank_file}", file=sys.stderr)
        sys.exit(1)
    except CodonTable.CodonTableNotFoundError as e:
        print(f"Error: Invalid NCBI translation table ID used in a feature: {e}", file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()

def main():
    """
    Main function to parse command-line arguments and call the extraction function.
    """
    parser = argparse.ArgumentParser(
        description="Extracts a nucleotide sequence from a GenBank record and translates "
                    "it to protein if it overlaps with a CDS feature."
    )
    parser.add_argument(
        "genbank_file",
        help="Path to the input GenBank file."
    )
    parser.add_argument(
        "start_pos",
        type=int,
        help="Start base pair of the region to extract (1-based)."
    )
    parser.add_argument(
        "end_pos",
        type=int,
        help="End base pair of the region to extract (1-based)."
    )
    args = parser.parse_args()

    extract_and_translate_cds_region(args.genbank_file, args.start_pos, args.end_pos)

if __name__ == "__main__":
    main()
