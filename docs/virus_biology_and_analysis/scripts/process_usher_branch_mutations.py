# The following are Bash shell commands that should be run before this script.
# Ensure the necessary files (usher_tree.nwk, usher_branch_mutations.tsv) are present.

# Download the latest MAT file (if not already done):
# wget http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.all.masked.pb.gz
# gunzip public-latest.all.masked.pb.gz # This creates public-latest.all.masked.pb

# Export the Newick tree (using the decompressed .pb file):
# matUtils extract -i public-latest.all.masked.pb -t usher_tree.nwk

# Dump the branch mutations (using the decompressed .pb file):
# matUtils dump -i public-latest.all.masked.pb -m usher_branch_mutations.tsv

# The rest of this file is Python code.

# Uncomment the following lines and run if you don't have these libraries installed
# !pip install dendropy pandas

import dendropy
import pandas as pd
from collections import defaultdict

# Define your input file names
tree_file = "usher_tree.nwk"
mutations_file = "usher_branch_mutations.tsv"

print(f"Loading tree from: {tree_file}")
tree = dendropy.Tree.get_from_path(tree_file, schema="newick", preserve_underscores=True)
# UShER trees are usually rooted, but if yours is unrooted, you might need to root it:
# tree.reroot_at_midpoint(update_bipartitions=False)
print(f"Tree loaded with {len(tree.nodes())} nodes and {len(tree.leaf_nodes())} leaf nodes.")

print(f"Loading mutations from: {mutations_file}")
# Assuming the TSV has two columns: Node_ID and Mutations (comma-separated string)
# Adjust 'sep' if your file is comma-separated or uses a different delimiter.
try:
    mutations_df = pd.read_csv(mutations_file, sep='\t', header=None, names=['Node_ID', 'Mutations'], engine='python')
except pd.errors.EmptyDataError:
    print(f"Warning: {mutations_file} is empty or not formatted as expected. Proceeding with an empty mutation map.")
    mutations_df = pd.DataFrame(columns=['Node_ID', 'Mutations'])
except FileNotFoundError:
    print(f"Error: The mutations file '{mutations_file}' was not found. Please ensure it exists in the correct path.")
    exit()


# Convert the DataFrame to a dictionary for faster lookup
# Node_ID will be the key, and a list of mutation strings will be the value
node_to_mutations = defaultdict(list)
for index, row in mutations_df.iterrows():
    node_id = str(row['Node_ID']) # Ensure ID is string for consistent lookup
    mutations_str = row['Mutations']
    if pd.notna(mutations_str) and str(mutations_str).strip(): # Check for NaN or empty/whitespace-only strings
        node_to_mutations[node_id] = [m.strip() for m in str(mutations_str).split(',')]

print(f"Loaded mutations for {len(node_to_mutations)} nodes/branches.")

# Dictionary to store mutations mapped by the node_id from the TSV file
final_branch_mutations_map = {}

# Iterate through all nodes in the tree
for node in tree.postorder_node_iter():
    # The Node_ID in usher_branch_mutations.tsv is the label of the node where the branch *ends*.
    # For leaf nodes, this is usually the sample name (taxon.label).
    # For internal nodes, UShER might assign unique IDs, or they might be unlabeled in the Newick.
    # This script uses node.taxon.label for leaves and node.label for internal nodes if available.
    # If internal nodes in your Newick are unlabeled, node.oid (Dendropy's internal object ID) can be used,
    # but mapping this to UShER's dump output might require more sophisticated matching if UShER's
    # internal node IDs in the dump file are not simple OIDs.

    current_node_id_str = None
    if node.is_leaf() and node.taxon:
        current_node_id_str = str(node.taxon.label)
    elif node.label: # Check if an internal node has an explicit label in the Newick
        current_node_id_str = str(node.label)
    else: # Fallback for unlabeled internal nodes (less reliable for matching dump output)
        current_node_id_str = str(node.oid)


    mutations_on_this_branch = node_to_mutations.get(current_node_id_str, [])
    node.annotations.add_new('branch_mutations', mutations_on_this_branch)
    final_branch_mutations_map[current_node_id_str] = mutations_on_this_branch

print("Mutations associated with tree nodes.")

# 4. Output and Further Analysis
output_file_name = "annotated_branch_mutations.txt"
print(f"\nWriting branch mutations to {output_file_name}...")

with open(output_file_name, "w") as f:
    for node_id, mutations in final_branch_mutations_map.items():
        if mutations: # Only write branches that actually have mutations
            f.write(f"Branch leading to Node: {node_id}\n")
            for mut in mutations:
                f.write(f"  - {mut}\n")
            f.write("\n")

print("Finished writing results.")

# --- Optional: Further Dendropy operations ---
print("\n--- Example of annotated nodes (first 5 with mutations) ---")
count = 0
for node in tree.postorder_node_iter(): # Iterate again to show annotations
    node_muts = node.annotations.get_value('branch_mutations')
    if node_muts:
        node_display_id = node.taxon.label if node.taxon else node.label if node.label else node.oid
        print(f"Node Label (or OID): {node_display_id}")
        print(f"  Mutations: {node_muts}")
        count += 1
    if count >= 5:
        break

# Example: If you want to export the tree with annotations
# (check dendropy documentation for specific format support like NeXML)
# annotated_tree_file = "usher_tree_with_annotations.nex"
# try:
#   tree.write(path=annotated_tree_file, schema="nexus", suppress_annotations=False, unquoted_underscores=True)
#   print(f"Annotated tree saved to {annotated_tree_file}")
# except Exception as e:
#   print(f"Could not save annotated tree as NEXUS: {e}")
#   print("Consider 'newick' schema if NEXUS fails, though Newick doesn't typically support rich annotations.")
