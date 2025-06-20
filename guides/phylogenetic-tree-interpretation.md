# Guide: Interpreting Phylogenetic Trees with Python

Phylogenetic trees are graphical representations of evolutionary relationships. This guide provides a conceptual overview and a practical Python script for common tree interpretation tasks, such as loading, re-rooting, and extracting subtrees using the **Biopython** library.

### Common Tree Formats
*   **Newick:** A simple, plain-text format using nested parentheses. e.g., `(A:0.1,B:0.2,(C:0.3,D:0.4));`
*   **Nexus:** A more structured format that can contain trees, data, and analysis commands.
*   **PhyloXML:** An XML-based format designed for rich annotation and data sharing.

---

## How to Use the Python Script

### Prerequisites
1.  Install Python 3.
2.  Install the Biopython library:
    ```bash
    pip install biopython
    ```

### Quickstart Guide
1.  Save the code from the [complete script section](#the-complete-python-script) below as a Python file (e.g., `tree_tools.py`).
2.  Run the script from your terminal:
    ```bash
    python tree_tools.py
    ```
3.  The script will run a default demonstration using a sample Newick string. It will show how to load, re-root, and extract a subtree, saving the results to new files (`subtree_CD.nwk`, `midpoint_rooted_tree.nwk`).
4.  To use your own data, modify the `if __name__ == "__main__":` block at the bottom of the script to load your tree file instead of the sample string.

---

## Key Tree Operations Explained

The provided script includes functions for several fundamental tree manipulations.

### 1. Loading and Displaying a Tree
*   **`load_tree()`**: Loads a tree from a file (e.g., `.nwk`, `.nexus`) or directly from a Newick string.
*   **`print_tree_info()`**: Calculates and displays basic statistics like the number of tips, internal nodes, and total branch length.
*   **`display_ascii_tree()`**: Renders a simple text-based tree in the console, which is great for quick structural checks.

### 2. Re-rooting a Tree
Re-rooting is often necessary to correctly represent evolutionary relationships, especially when starting with an unrooted tree. The `reroot_tree()` function supports two common methods:
> **By Outgroup:** This is the preferred biological method. You provide one or more "outgroup" taxa that are known to be less related to the "ingroup." The tree is rooted on the branch leading to the common ancestor of the outgroup.
>
> **At Midpoint:** This is an algorithmic method used when a clear outgroup is unknown. It places the root at the midpoint of the longest path between any two tips in the tree. This requires meaningful branch lengths.

### 3. Extracting a Subtree
*   **`extract_subtree()`**: Creates a new, smaller tree containing only a specified subset of tips and their common ancestors. This is extremely useful for focusing on a specific clade or group of interest within a larger phylogeny. The new subtree is rooted at the Most Recent Common Ancestor (MRCA) of the specified tips.

### 4. Saving a Tree
*   **`save_tree()`**: Saves your tree object (whether original or modified) to a file in a specified format (Newick, PhyloXML, etc.).

---

## The Complete Python Script

```python
#!/usr/bin/env python3

import io
from Bio import Phylo
import copy

# --- Helper Functions ---
def print_tree_info(tree, tree_name="Tree"):
    """Prints basic information about the tree."""
    if not tree:
        print(f"[{tree_name}] No tree data to display.")
        return
    terminals = tree.get_terminals()
    total_branch_length = sum(c.branch_length for c in tree.find_clades() if c.branch_length is not None)
    
    print(f"\n--- {tree_name} Information ---")
    print(f"  - Tips (Terminals): {len(terminals)}")
    print(f"  - Internal Nodes: {len(tree.get_nonterminals())}")
    print(f"  - Tree is Rooted: {tree.rooted}")
    print(f"  - Total Branch Length: {total_branch_length:.4f}")

def display_ascii_tree(tree, tree_name="Tree"):
    """Displays a simple ASCII representation of the tree."""
    if not tree:
        return
    print(f"\n--- {tree_name} ASCII Representation ---")
    Phylo.draw_ascii(tree)

# --- Core Tree Operations ---
def load_tree(filepath_or_string, file_format="newick"):
    """Loads a tree from a file or a string."""
    try: # Try loading from a file first
        tree = Phylo.read(filepath_or_string, file_format)
        print(f"[Loader] Successfully loaded tree from file: {filepath_or_string}")
        return tree
    except FileNotFoundError: # If that fails, try loading from a string
        try:
            tree = Phylo.read(io.StringIO(filepath_or_string), file_format)
            print("[Loader] Successfully loaded tree from string input.")
            return tree
        except Exception as e_str:
            print(f"[Loader] Error: Could not parse input as a tree string. {e_str}")
            return None
    except Exception as e_file:
        print(f"[Loader] Error: Could not load tree from file. {e_file}")
        return None

def reroot_tree(tree, outgroup_names=None, at_midpoint=False):
    """Re-roots the tree by outgroup or midpoint."""
    if not tree: return None
    
    # Always work on a copy to avoid modifying the original object in place
    tree_copy = copy.deepcopy(tree)

    if at_midpoint:
        print("[Reroot] Re-rooting tree at midpoint...")
        tree_copy.root_at_midpoint()
        return tree_copy
    elif outgroup_names:
        print(f"[Reroot] Re-rooting with outgroup: {outgroup_names}...")
        try:
            tree_copy.root_with_outgroup(*outgroup_names)
            return tree_copy
        except ValueError as e:
            print(f"[Reroot] Error: Outgroup tip not found. {e}")
            return None
    else:
        print("[Reroot] No valid rooting criteria specified (outgroup or midpoint).")
        return tree

def extract_subtree(tree, tip_names):
    """Extracts a subtree containing the specified tip names."""
    if not tree or not tip_names: return None
    
    print(f"[Subtree] Extracting subtree for tips: {tip_names}...")
    try:
        # Find the Most Recent Common Ancestor (MRCA) of the target tips
        mrca = tree.common_ancestor(*tip_names)
        # Create a new tree rooted at this MRCA
        subtree = Phylo.BaseTree.Tree(root=mrca, rooted=True)
        return subtree
    except ValueError as e:
        print(f"[Subtree] Error: One or more tips not found. {e}")
        return None

def save_tree(tree, output_filepath, file_format="newick"):
    """Saves a tree to a file."""
    if not tree: return
    print(f"[Save] Saving tree to '{output_filepath}' (format: {file_format})...")
    Phylo.write(tree, output_filepath, file_format)

# --- Main Execution Example ---
if __name__ == "__main__":
    # A sample Newick string representing a simple tree
    sample_newick = "((TipA:0.1,TipB:0.2):0.05,(TipC:0.3,TipD:0.4):0.08);"
    
    print("--- 1. Loading and Inspecting Original Tree ---")
    original_tree = load_tree(sample_newick, "newick")
    if original_tree:
        print_tree_info(original_tree, "Original Tree")
        display_ascii_tree(original_tree, "Original Tree")

        # --- 2. Rerooting Examples ---
        rooted_on_a = reroot_tree(original_tree, outgroup_names=["TipA"])
        if rooted_on_a:
            print_tree_info(rooted_on_a, "Tree Rooted on TipA")
            display_ascii_tree(rooted_on_a, "Tree Rooted on TipA")

        midpoint_rooted_tree = reroot_tree(original_tree, at_midpoint=True)
        if midpoint_rooted_tree:
            print_tree_info(midpoint_rooted_tree, "Tree Rooted at Midpoint")
            display_ascii_tree(midpoint_rooted_tree, "Tree Rooted at Midpoint")
            save_tree(midpoint_rooted_tree, "midpoint_rooted_tree.nwk")

        # --- 3. Subtree Extraction Example ---
        subtree_cd = extract_subtree(original_tree, tip_names=["TipC", "TipD"])
        if subtree_cd:
            print_tree_info(subtree_cd, "Subtree C-D")
            display_ascii_tree(subtree_cd, "Subtree C-D")
            save_tree(subtree_cd, "subtree_CD.nwk")

    print("\nScript finished.")
```

---

## Beyond the Script: Further Analysis

This script provides a solid foundation. Once you have a tree loaded with `Bio.Phylo`, you can perform many other advanced operations:

*   **Traverse the tree:** Iterate over clades to access parent/child relationships.
*   **Access node attributes:** Read and write data like `clade.name`, `clade.branch_length`, and `clade.confidence` (bootstrap values).
*   **Compare trees:** Calculate topological distances between different trees.
*   **Advanced Visualization:** Export your final tree to Newick or PhyloXML format and use dedicated visualization software like [FigTree](http://tree.bio.ed.ac.uk/software/figtree/), [iTOL](https://itol.embl.de/), or the R package `ggtree` for publication-quality figures.
