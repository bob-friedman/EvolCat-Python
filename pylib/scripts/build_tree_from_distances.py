#!/usr/bin/env python3

"""
Constructs a phylogenetic tree from a distance matrix using Neighbor-Joining (NJ)
or UPGMA methods.
"""

import argparse
import sys
# import os # Not strictly used, but often useful.
from Bio import Phylo
from Bio.Nexus import Nexus
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

def parse_phylip_distance_matrix(filepath):
    """
    Parses a PHYLIP-formatted distance matrix (lower triangle without diagonal).
    Format:
    N
    Taxon1
    Taxon2  d(2,1)
    Taxon3  d(3,1)  d(3,2)
    ...

    Returns:
        tuple: (list_of_names, dict_of_distances) or (None, None) if parsing fails.
               The dict_of_distances is {(name_i, name_j): distance}
               where name_i is the taxon of the current row and name_j is a taxon from a previous column.
    """
    names = []
    distances_dict = {}
    try:
        with open(filepath, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]

            if not lines:
                # print("DEBUG PHYLIP: File is empty or contains only whitespace.", file=sys.stderr)
                return None, None

            try:
                num_taxa_declared = int(lines[0])
            except ValueError:
                # print(f"DEBUG PHYLIP: First line '{lines[0]}' is not a valid integer for number of taxa.", file=sys.stderr)
                return None, None

            # Expected number of lines: 1 (count) + num_taxa_declared (taxa lines)
            if len(lines) < num_taxa_declared + 1:
                # print(f"DEBUG PHYLIP: Declared {num_taxa_declared} taxa, but found only {len(lines)-1} data lines after count.", file=sys.stderr)
                return None, None

            actual_taxa_lines = lines[1 : num_taxa_declared + 1]
            if len(actual_taxa_lines) != num_taxa_declared: # Should be redundant due to above check, but good for safety
                # print(f"DEBUG PHYLIP: Expected {num_taxa_declared} taxon lines, found {len(actual_taxa_lines)}.", file=sys.stderr)
                return None, None

            for i, line_content in enumerate(actual_taxa_lines):
                parts = line_content.split()
                if not parts:
                    # print(f"DEBUG PHYLIP: Taxon line {i+1} is empty.", file=sys.stderr)
                    return None, None

                taxon_name = parts[0]
                names.append(taxon_name)
                
                expected_distances_on_line = i # 0 for 1st taxon (names[0]), 1 for 2nd (names[1]), etc.
                actual_distances_on_line = len(parts) - 1

                if actual_distances_on_line != expected_distances_on_line:
                    # print(f"DEBUG PHYLIP: Taxon '{taxon_name}': expected {expected_distances_on_line} distances, found {actual_distances_on_line}.", file=sys.stderr)
                    return None, None

                for j in range(actual_distances_on_line):
                    try:
                        dist_val = float(parts[j+1])
                    except ValueError:
                        # print(f"DEBUG PHYLIP: Invalid distance value '{parts[j+1]}' for pair involving '{taxon_name}' and '{names[j]}'.", file=sys.stderr)
                        return None, None
                    
                    # names[i] is the current taxon_name. names[j] is the taxon from a previous column.
                    # The key order (taxon_in_row, taxon_in_column) is (names[i], names[j])
                    distances_dict[(taxon_name, names[j])] = dist_val
            
            if len(names) == num_taxa_declared:
                # print(f"DEBUG PHYLIP: Successfully parsed {num_taxa_declared} taxa and their distances.", file=sys.stderr)
                return names, distances_dict
            else:
                # This case should ideally not be reached if logic is sound.
                # print(f"DEBUG PHYLIP: Parsing completed, but number of names found ({len(names)}) does not match declared ({num_taxa_declared}).", file=sys.stderr)
                return None, None

    except Exception as e: # Catches general errors during PHYLIP parsing
        # print(f"DEBUG PHYLIP: An unexpected error occurred during PHYLIP parsing: {e}", file=sys.stderr)
        return None, None


def main():
    """Main function to build tree from distances."""
    parser = argparse.ArgumentParser(
        description="Constructs a phylogenetic tree from a distance matrix.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "distance_matrix_file",
        help="Path to the input distance matrix file."
    )
    parser.add_argument(
        "--method",
        choices=["nj", "upgma"],
        default="nj",
        help="Tree construction method: 'nj' (Neighbor-Joining) or 'upgma' (UPGMA).\nDefault is 'nj'."
    )
    parser.add_argument(
        "--outfile",
        default="phylogenetic_tree.nwk",
        help="Output file name for the tree in Newick format.\nDefault is 'phylogenetic_tree.nwk'."
    )
    parser.add_argument(
        "--informat",
        choices=["nexus"], # Currently only one choice, implies Bio.Nexus.Nexus parser + PHYLIP fallback
        default="nexus",
        help="Input distance matrix file format.\nCurrently supports 'nexus' (for PHYLIP-style matrices, standalone or within a Nexus DISTANCES block).\nDefault is 'nexus'."
    )

    args = parser.parse_args()

    names = None
    nexus_matrix_values = None 

    try:
        # Attempt 1: Parse as Nexus using Bio.Nexus.Nexus
        # print(f"DEBUG: Attempting to parse '{args.distance_matrix_file}' with Bio.Nexus.Nexus", file=sys.stderr)
        nex_obj = Nexus.Nexus(args.distance_matrix_file)

        # --- Debugging prints (can be minimized later) ---
        # print(f"DEBUG: dir(nex_obj) = {dir(nex_obj)}", file=sys.stderr)
        # if hasattr(nex_obj, 'taxlabels'): print(f"DEBUG: nex_obj.taxlabels = {nex_obj.taxlabels}", file=sys.stderr)
        # if hasattr(nex_obj, 'matrix'): print(f"DEBUG: nex_obj.matrix (type) = {type(nex_obj.matrix)}, value = {nex_obj.matrix}", file=sys.stderr)
        # if hasattr(nex_obj, 'distances') and nex_obj.distances:
        #     dist_block_debug = nex_obj.distances
        #     if hasattr(dist_block_debug, 'names'): print(f"DEBUG: nex_obj.distances.names = {dist_block_debug.names}", file=sys.stderr)
        #     if hasattr(dist_block_debug, 'matrix'): print(f"DEBUG: nex_obj.distances.matrix (type) = {type(dist_block_debug.matrix)}", file=sys.stderr)
        # else:
        #     print(f"DEBUG: nex_obj.distances is None or does not exist.", file=sys.stderr)
        # --- End debugging prints ---

        # Try to get data from NEXUS DISTANCES block
        if hasattr(nex_obj, 'distances') and nex_obj.distances is not None:
            dist_block = nex_obj.distances
            if hasattr(dist_block, 'names') and dist_block.names:
                names = dist_block.names
            elif hasattr(nex_obj, 'taxlabels') and nex_obj.taxlabels: # Fallback for names
                names = nex_obj.taxlabels
            
            if hasattr(dist_block, 'matrix') and dist_block.matrix is not None:
                nexus_matrix_values = dist_block.matrix
            
            if names and nexus_matrix_values:
                 print(f"INFO: Successfully read data from NEXUS DISTANCES block.", file=sys.stderr)

        # If not found in DISTANCES, try PHYLIP-like interpretation from top-level nex_obj attributes
        if not names or nexus_matrix_values is None:
            # print(f"DEBUG: Data not fully extracted from DISTANCES block. Attempting PHYLIP-like interpretation via nex_obj.taxlabels/matrix.", file=sys.stderr)
            current_names = None
            if hasattr(nex_obj, 'taxlabels') and nex_obj.taxlabels:
                current_names = nex_obj.taxlabels
            
            current_matrix_data = None
            if hasattr(nex_obj, 'matrix') and nex_obj.matrix is not None:
                m = nex_obj.matrix
                if isinstance(m, list) and m and all(isinstance(row, list) for row in m):
                    if current_names and len(m) == len(current_names) and all(len(row) == len(current_names) for row in m):
                        current_matrix_data = m
            
            if current_names and current_matrix_data:
                if not names: names = current_names
                if nexus_matrix_values is None: nexus_matrix_values = current_matrix_data
                print(f"INFO: Successfully read data using Bio.Nexus.Nexus PHYLIP-like interpretation (nex_obj.taxlabels, nex_obj.matrix).", file=sys.stderr)

    except FileNotFoundError:
        print(f"Error: Input file '{args.distance_matrix_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except Nexus.NexusError as e: 
        print(f"Warning: Could not parse file '{args.distance_matrix_file}' as Nexus using Bio.Nexus.Nexus. Details: {e}", file=sys.stderr)
        # Don't exit yet, will try PHYLIP parser
    except Exception as e: 
        print(f"An unexpected error occurred while attempting Bio.Nexus.Nexus parsing of '{args.distance_matrix_file}': {e}", file=sys.stderr)
        # Don't exit yet, will try PHYLIP parser

    # Attempt 2: If Nexus parsing failed or yielded no data, try dedicated PHYLIP parser
    if not names or nexus_matrix_values is None:
        print(f"INFO: Bio.Nexus.Nexus parsing did not yield complete data. Attempting dedicated PHYLIP distance matrix parsing for '{args.distance_matrix_file}'.", file=sys.stderr)
        phylip_names, phylip_matrix_dict = parse_phylip_distance_matrix(args.distance_matrix_file)
        
        if phylip_names and phylip_matrix_dict:
            names = phylip_names
            nexus_matrix_values = phylip_matrix_dict # This is a dict, handled by downstream code
            print(f"INFO: Successfully parsed '{args.distance_matrix_file}' as PHYLIP distance matrix.", file=sys.stderr)
        else:
            print(f"Error: Failed to extract taxon names and/or matrix values from '{args.distance_matrix_file}' using either Nexus or PHYLIP parsing methods.", file=sys.stderr)
            # print(f"DEBUG FINAL: Extracted names: {names}", file=sys.stderr)
            # print(f"DEBUG FINAL: Extracted matrix data (type): {type(nexus_matrix_values)}", file=sys.stderr)
            sys.exit(1)


    if not names or nexus_matrix_values is None: # Should be caught by the above, but as a final check
        print(f"Error: Critical failure in data extraction. Names: {names}, Matrix Type: {type(nexus_matrix_values)}", file=sys.stderr)
        sys.exit(1)

    num_taxa = len(names)
    lower_tri_matrix = []

    try:
        if isinstance(nexus_matrix_values, dict):
            # Case 1: Matrix is a dictionary {('TaxonA', 'TaxonB'): distance}
            # This is expected from nex_obj.distances.matrix or our PHYLIP parser
            for i in range(num_taxa):
                row = []
                for j in range(i + 1): 
                    if i == j:
                        row.append(0.0)
                    else:
                        name_i, name_j = names[i], names[j]
                        dist = nexus_matrix_values.get((name_i, name_j))
                        if dist is None:
                            dist = nexus_matrix_values.get((name_j, name_i))
                        
                        if dist is None: # Index-based fallback (less likely needed now)
                            dist = nexus_matrix_values.get((i,j))
                        if dist is None:
                            dist = nexus_matrix_values.get((j,i))

                        if dist is None:
                            raise ValueError(f"Distance not found between '{name_i}' (idx {i}) and '{name_j}' (idx {j})")
                        row.append(float(dist))
                lower_tri_matrix.append(row)
        
        elif isinstance(nexus_matrix_values, list) and all(isinstance(r, list) for r in nexus_matrix_values):
            # Case 2: Matrix is a list of lists (e.g. from nex_obj.matrix for some PHYLIP-like files)
            if len(nexus_matrix_values) == num_taxa and all(len(r) == (idx + 1) for idx, r in enumerate(nexus_matrix_values)):
                lower_tri_matrix = [[float(d) for d in r] for r in nexus_matrix_values]
            elif len(nexus_matrix_values) == num_taxa and all(len(r) == num_taxa for r in nexus_matrix_values):
                for i in range(num_taxa):
                    row = []
                    for j in range(i + 1):
                        row.append(float(nexus_matrix_values[i][j]))
                    lower_tri_matrix.append(row)
            else:
                raise ValueError("Matrix format (list of lists) is not a recognized full square or lower-triangular with diagonal matching taxon count.")
        else:
            raise ValueError(f"Unsupported matrix data type from input: {type(nexus_matrix_values)}. Expected dict or list of lists.")

    except ValueError as e:
        print(f"Error processing distance matrix values: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred during matrix conversion: {e}", file=sys.stderr)
        sys.exit(1)
    
    try:
        # print(f"DEBUG: Final names for DistanceMatrix: {names}", file=sys.stderr)
        # print(f"DEBUG: Final lower_tri_matrix for DistanceMatrix: {lower_tri_matrix}", file=sys.stderr)
        dm = DistanceMatrix(names, lower_tri_matrix)
    except Exception as e:
        print(f"Error creating Bio.Phylo.TreeConstruction.DistanceMatrix object: {e}", file=sys.stderr)
        sys.exit(1)

    constructor = DistanceTreeConstructor()
    tree = None
    try:
        if args.method == "nj":
            tree = constructor.nj(dm)
        elif args.method == "upgma":
            tree = constructor.upgma(dm)
    except Exception as e:
        print(f"Error during tree construction ({args.method}): {e}", file=sys.stderr)
        sys.exit(1)

    if tree is None:
        print(f"Error: Tree construction failed for an unknown reason.", file=sys.stderr)
        sys.exit(1)
        
    try:
        Phylo.write(tree, args.outfile, "newick")
    except Exception as e:
        print(f"Error writing tree to '{args.outfile}': {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Successfully built tree using '{args.method}' method.")
    print(f"Tree written to '{args.outfile}'.")

if __name__ == '__main__':
    main()
