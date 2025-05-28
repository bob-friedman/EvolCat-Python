#!/usr/bin/env python3

"""
Constructs a phylogenetic tree from a distance matrix using Neighbor-Joining (NJ)
or UPGMA methods.
"""

import argparse
import sys
import os
from Bio import Phylo
from Bio.Nexus import Nexus
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

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
        choices=["nexus"],
        default="nexus",
        help="Input distance matrix file format.\nCurrently supports 'nexus' (for PHYLIP-style matrices, standalone or within a Nexus DISTANCES block).\nDefault is 'nexus'."
    )

    args = parser.parse_args()

    # Read and parse the distance matrix
    names = None
    nexus_matrix_values = None 
    try:
        nex_obj = Nexus.Nexus(args.distance_matrix_file)

        # --- Temporary debugging prints ---
        print(f"DEBUG: dir(nex_obj) = {dir(nex_obj)}", file=sys.stderr)
        if hasattr(nex_obj, 'names'):
            print(f"DEBUG: nex_obj.names = {nex_obj.names}", file=sys.stderr)
        else:
            print("DEBUG: nex_obj has no 'names' attribute", file=sys.stderr)
        if hasattr(nex_obj, 'taxlabels'):
            print(f"DEBUG: nex_obj.taxlabels = {nex_obj.taxlabels}", file=sys.stderr)
        else:
            print("DEBUG: nex_obj has no 'taxlabels' attribute", file=sys.stderr)
        if hasattr(nex_obj, 'matrix'):
            print(f"DEBUG: nex_obj.matrix (type) = {type(nex_obj.matrix)}", file=sys.stderr)
            if nex_obj.matrix is not None:
                 print(f"DEBUG: nex_obj.matrix = {nex_obj.matrix}", file=sys.stderr)
        else:
            print("DEBUG: nex_obj has no 'matrix' attribute", file=sys.stderr)
        if hasattr(nex_obj, 'unaltered_matrix'):
            print(f"DEBUG: nex_obj.unaltered_matrix (type) = {type(nex_obj.unaltered_matrix)}", file=sys.stderr)
            if nex_obj.unaltered_matrix is not None:
                print(f"DEBUG: nex_obj.unaltered_matrix = {nex_obj.unaltered_matrix}", file=sys.stderr)
        else:
            print("DEBUG: nex_obj has no 'unaltered_matrix' attribute", file=sys.stderr)
        # --- End temporary debugging prints ---

        # Prioritize taxlabels and matrix for PHYLIP, then names and unaltered_matrix for Nexus DISTANCES
        if hasattr(nex_obj, 'taxlabels') and nex_obj.taxlabels and \
           hasattr(nex_obj, 'matrix') and nex_obj.matrix is not None:
            names = nex_obj.taxlabels
            nexus_matrix_values = nex_obj.matrix 
            # print(f"DEBUG: Used nex_obj.taxlabels and nex_obj.matrix", file=sys.stderr)
        elif hasattr(nex_obj, 'names') and nex_obj.names and \
             hasattr(nex_obj, 'unaltered_matrix') and nex_obj.unaltered_matrix is not None:
            names = nex_obj.names
            nexus_matrix_values = nex_obj.unaltered_matrix
            # print(f"DEBUG: Used nex_obj.names and nex_obj.unaltered_matrix", file=sys.stderr)
        else:
            print(f"Error: Could not extract matrix and/or taxon names from '{args.distance_matrix_file}' using common Bio.Nexus attributes.", file=sys.stderr)
            sys.exit(1)

    except FileNotFoundError:
        print(f"Error: Input file '{args.distance_matrix_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except Nexus.NexusError as e: 
        print(f"Error: Could not parse file '{args.distance_matrix_file}' as Nexus/PHYLIP. Details: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e: 
        print(f"An unexpected error occurred while reading '{args.distance_matrix_file}': {e}", file=sys.stderr)
        sys.exit(1)

    if not names or nexus_matrix_values is None:
        print(f"Error: Failed to extract names or matrix values from '{args.distance_matrix_file}'. Names: {names}, Matrix Value Type: {type(nexus_matrix_values)}", file=sys.stderr)
        sys.exit(1)

    num_taxa = len(names)
    lower_tri_matrix = []

    # The 'nexus_matrix_values' can be a dictionary (taxon_name_pairs -> dist) 
    # or a list of lists (raw matrix for PHYLIP-like).
    # We need to construct the lower_tri_matrix [[0], [d_BA, 0], [d_CA, d_CB, 0]]

    try:
        if isinstance(nexus_matrix_values, dict):
            # Case 1: Matrix is a dictionary {('TaxonA', 'TaxonB'): distance}
            # Or it could be {(idxA, idxB): distance} if Bio.Nexus processes it that way.
            # We need to ensure we map names to indices correctly.
            name_to_idx = {name: i for i, name in enumerate(names)}
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
                        # Also check for index-based keys if name-based keys fail
                        if dist is None:
                            dist = nexus_matrix_values.get((i,j))
                        if dist is None:
                            dist = nexus_matrix_values.get((j,i))

                        if dist is None:
                            raise ValueError(f"Distance not found between '{name_i}' and '{name_j}'")
                        row.append(float(dist))
                lower_tri_matrix.append(row)
        elif isinstance(nexus_matrix_values, list) and all(isinstance(r, list) for r in nexus_matrix_values):
            # Case 2: Matrix is already a list of lists (expected for PHYLIP-like input that Bio.Nexus might return as nex_obj.matrix)
            # Ensure it's what DistanceMatrix expects: lower triangular with diagonal.
            if len(nexus_matrix_values) == num_taxa and all(len(r) == (idx + 1) for idx, r in enumerate(nexus_matrix_values)):
                # Already in the correct [[0], [d_BA, 0], ...] format
                lower_tri_matrix = [[float(d) for d in r] for r in nexus_matrix_values]
            elif len(nexus_matrix_values) == num_taxa and all(len(r) == num_taxa for r in nexus_matrix_values):
                # It's a full matrix, extract lower triangle
                for i in range(num_taxa):
                    row = []
                    for j in range(i + 1):
                        row.append(float(nexus_matrix_values[i][j]))
                    lower_tri_matrix.append(row)
            else: # PHYLIP lower triangle without diagonal (as in the input file)
                  # This is what nex_obj.unaltered_matrix usually gives for PHYLIP input.
                  # Example: [[], [0.2], [0.5, 0.4], [0.8, 0.7, 0.6]] for a 4-taxa matrix
                  # The first element is the number of taxa (e.g., 4)
                  # The subsequent elements are names, then the matrix rows.
                  # This structure is for PHYLIP files *directly*, not after Nexus parsing usually.
                  # The Bio.Nexus.Nexus object should abstract this away.
                  # If nexus_matrix_values is from nex_obj.unaltered_matrix for a PHYLIP file,
                  # it might be a list of lists where each inner list is a row of the lower triangle *without* diagonal.
                  # e.g. [[0.20], [0.50, 0.40], [0.80, 0.70, 0.60]] for a 4-taxa matrix (after names)
                  # Let's assume nexus_matrix_values from Bio.Nexus is already somewhat processed.
                  # If it's the direct unaltered_matrix from a PHYLIP, it's tricky.
                  # The current structure of sample_distance_matrix.phy is NEXUS.
                  # So, nex_obj.unaltered_matrix or nex_obj.distances.matrix should provide a dict or list of lists.
                  # The most likely for a NEXUS DISTANCES block is a dictionary.
                raise ValueError("Matrix format from Nexus object is not a recognized list of lists (full or lower-triangular with diagonal).")

        else:
            raise ValueError(f"Unsupported matrix format from Nexus object: {type(nexus_matrix_values)}")

    except ValueError as e:
        print(f"Error processing distance matrix: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred during matrix conversion: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Create DistanceMatrix object
    try:
        dm = DistanceMatrix(names, lower_tri_matrix)
    except Exception as e:
        print(f"Error creating DistanceMatrix object: {e}", file=sys.stderr)
        sys.exit(1)

    # Create DistanceTreeConstructor object
    constructor = DistanceTreeConstructor()

    # Build the tree
    tree = None
    try:
        if args.method == "nj":
            tree = constructor.nj(dm)
        elif args.method == "upgma":
            tree = constructor.upgma(dm)
        else:
            # Should not happen due to argparse choices
            print(f"Error: Unknown tree construction method '{args.method}'.", file=sys.stderr)
            sys.exit(1)
    except Exception as e:
        print(f"Error during tree construction ({args.method}): {e}", file=sys.stderr)
        sys.exit(1)

    if tree is None:
        print(f"Error: Tree construction failed for an unknown reason.", file=sys.stderr)
        sys.exit(1)
        
    # Write the tree to the output file
    try:
        Phylo.write(tree, args.outfile, "newick")
    except Exception as e:
        print(f"Error writing tree to '{args.outfile}': {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Successfully built tree using '{args.method}' method.")
    print(f"Tree written to '{args.outfile}'.")

if __name__ == '__main__':
    main()
