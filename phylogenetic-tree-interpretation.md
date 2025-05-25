# Interpreting Phylogenetic Trees with Python

Phylogenetic trees are graphical representations of the evolutionary relationships among a group of organisms, genes, or other entities. Interpreting these trees is fundamental to understanding evolutionary history, relatedness, and patterns of diversification. This guide provides a Python script to help with common tree interpretation tasks, such as loading trees, re-rooting, and extracting subtrees.

**Common Phylogenetic Tree Formats:**

*   **Newick (or New Hampshire format):** A widely used format using nested parentheses to represent tree topology, often including branch lengths and node labels. Example: `(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;`
*   **Nexus:** A more structured format that can contain trees, data matrices, and commands for analysis. Trees within Nexus files often use a Newick-like representation.
*   **PhyloXML:** An XML-based format designed to share and store phylogenetic trees along with associated data.

Our Python script will primarily leverage the `Bio.Phylo` module from Biopython, which can handle these formats and more.

## Python Script for Phylogenetic Tree Interpretation

```python
#!/usr/bin/env python3

import io # For handling string IO for Newick examples
from Bio import Phylo

# --- Helper Functions ---
def print_tree_info(tree, tree_name="Tree"):
    """Prints basic information about the tree."""
    if not tree:
        print(f"[{tree_name}] No tree data to display.")
        return
    terminals = tree.get_terminals()
    non_terminals = tree.get_nonterminals()
    print(f"\n--- {tree_name} Information ---")
    print(f"Number of tips (terminals): {len(terminals)}")
    print(f"Number of internal nodes: {len(non_terminals)}")
    if tree.rooted:
        print("Tree is rooted.")
        if tree.root.branch_length is not None:
             print(f"Root branch length: {tree.root.branch_length}")
    else:
        print("Tree is unrooted.")
    
    total_branch_length = sum(clade.branch_length for clade in tree.find_clades() if clade.branch_length is not None)
    print(f"Total branch length (if available): {total_branch_length:.4f}")
    # print("Terminals:", [t.name for t in terminals]) # Uncomment for list of tip names

def display_ascii_tree(tree, tree_name="Tree"):
    """Displays a simple ASCII representation of the tree."""
    if not tree:
        print(f"[{tree_name}] No tree data to display.")
        return
    print(f"\n--- {tree_name} ASCII Representation ---")
    Phylo.draw_ascii(tree)

# --- Core Tree Operations ---
def load_tree(filepath_or_string, file_format="newick"):
    """
    Loads a phylogenetic tree from a file or a string.

    Args:
        filepath_or_string (str): Path to the tree file or a string containing the tree.
        file_format (str): Format of the tree (e.g., "newick", "nexus", "phyloxml").

    Returns:
        Bio.Phylo.BaseTree.Tree: A Biopython tree object, or None if loading fails.
    """
    try:
        # Try to read as a file first
        with open(filepath_or_string, 'r') as f:
            tree = Phylo.read(f, file_format)
        print(f"[Loader] Successfully loaded tree from file: {filepath_or_string}")
        return tree
    except FileNotFoundError:
        # If file not found, try to interpret as a string
        try:
            tree_string_io = io.StringIO(filepath_or_string)
            tree = Phylo.read(tree_string_io, file_format)
            print(f"[Loader] Successfully loaded tree from string input.")
            return tree
        except Exception as e_str:
            print(f"[Loader] Error loading tree from string: {e_str}")
            return None
    except Exception as e_file:
        print(f"[Loader] Error loading tree from file '{filepath_or_string}': {e_file}")
        return None

def reroot_tree(tree, outgroup_tip_names=None, at_midpoint=False):
    """
    Re-roots the tree.

    Args:
        tree (Bio.Phylo.BaseTree.Tree): The tree to re-root.
        outgroup_tip_names (list of str, optional): List of tip names forming the outgroup.
                                                    The tree will be rooted on the MRCA of these tips.
        at_midpoint (bool, optional): If True, roots the tree at its midpoint.
                                      Overrides outgroup_tip_names if both are specified.
    Returns:
        Bio.Phylo.BaseTree.Tree: The re-rooted tree (modifies the original tree object too).
    """
    if not tree:
        print("[Reroot] No tree provided.")
        return None

    if at_midpoint:
        try:
            tree.root_at_midpoint()
            print("[Reroot] Tree re-rooted at midpoint.")
            return tree
        except Exception as e:
            print(f"[Reroot] Error rooting at midpoint: {e}")
            return tree # Return original tree if error
    elif outgroup_tip_names and isinstance(outgroup_tip_names, list) and len(outgroup_tip_names) > 0:
        try:
            # Find the outgroup clade(s)
            outgroup_clades = []
            for name in outgroup_tip_names:
                clade = next(tree.find_clades(name), None)
                if clade:
                    outgroup_clades.append(clade)
                else:
                    print(f"[Reroot] Warning: Outgroup tip '{name}' not found in tree.")
            
            if not outgroup_clades:
                print("[Reroot] No valid outgroup tips found. Cannot reroot.")
                return tree

            if len(outgroup_clades) == 1:
                tree.root_with_outgroup(outgroup_clades[0])
                print(f"[Reroot] Tree re-rooted with outgroup: {outgroup_tip_names[0]}")
            else: # Multiple tips for outgroup, root on their MRCA
                mrca_outgroup = tree.common_ancestor(*outgroup_clades)
                tree.root_with_outgroup(mrca_outgroup)
                print(f"[Reroot] Tree re-rooted with MRCA of outgroup tips: {outgroup_tip_names}")
            return tree
        except Exception as e:
            print(f"[Reroot] Error re-rooting with outgroup '{outgroup_tip_names}': {e}")
            return tree # Return original tree if error
    else:
        print("[Reroot] No valid rooting criteria specified (outgroup or midpoint).")
        return tree

def extract_subtree(tree, tip_names_for_subtree):
    """
    Extracts a subtree containing the specified tip names (and their MRCA).

    Args:
        tree (Bio.Phylo.BaseTree.Tree): The original tree.
        tip_names_for_subtree (list of str): A list of tip names to include in the subtree.
                                             The subtree will be rooted at their MRCA.
    Returns:
        Bio.Phylo.BaseTree.Tree: A new tree object representing the subtree, or None.
    """
    if not tree:
        print("[Subtree] No tree provided.")
        return None
    if not tip_names_for_subtree or not isinstance(tip_names_for_subtree, list) or len(tip_names_for_subtree) < 1:
        print("[Subtree] No tip names provided for subtree extraction.")
        return None

    try:
        target_clades = []
        for name in tip_names_for_subtree:
            clade = next(tree.find_clades(name), None)
            if clade:
                target_clades.append(clade)
            else:
                print(f"[Subtree] Warning: Tip '{name}' not found in tree for subtree extraction.")
        
        if not target_clades:
            print("[Subtree] None of the specified tips found. Cannot extract subtree.")
            return None
        
        if len(target_clades) == 1 and target_clades[0].is_terminal(): # Subtree of a single tip is just the tip
             # Create a new minimal tree for a single tip
            subtree_root = Phylo.BaseTree.Clade(name=target_clades[0].name)
            subtree = Phylo.BaseTree.Tree(root=subtree_root, rooted=True) # Treat as rooted
            print(f"[Subtree] Extracted subtree for single tip: {target_clades[0].name}")
            return subtree

        # Find the Most Recent Common Ancestor (MRCA) of the target clades
        mrca = tree.common_ancestor(*target_clades)
        
        # Create a new tree from this MRCA clade
        # Phylo.BaseTree.Tree() constructor can take a Clade object as root
        subtree = Phylo.BaseTree.Tree(root=mrca, rooted=True) # The extracted part is rooted at the MRCA
        print(f"[Subtree] Extracted subtree for tips: {tip_names_for_subtree}")

        # Prune branches not leading to the desired tips (important if MRCA has other descendants)
        # This is implicitly handled by creating a new Tree from the MRCA clade if that's what's desired.
        # For a more explicit prune:
        all_subtree_tips = {c.name for c in subtree.get_terminals()}
        original_target_tips = set(tip_names_for_subtree)
        
        tips_to_prune_from_subtree = all_subtree_tips - original_target_tips
        
        for tip_name_to_prune in tips_to_prune_from_subtree:
            try:
                subtree.prune(tip_name_to_prune)
            except Exception as e_prune:
                # This can happen if the tip is part of the direct path to another target tip
                # or if the pruning logic in BioPython has specific constraints.
                # For simple MRCA-based subtrees, this step might not always be needed
                # if the goal is just the clade descending from the MRCA of the targets.
                print(f"[Subtree] Note: Could not prune '{tip_name_to_prune}' from subtree or it's not necessary: {e_prune}")

        return subtree
    except Exception as e:
        print(f"[Subtree] Error extracting subtree for tips '{tip_names_for_subtree}': {e}")
        return None

def save_tree(tree, output_filepath, file_format="newick"):
    """
    Saves a tree to a file.

    Args:
        tree (Bio.Phylo.BaseTree.Tree): The tree to save.
        output_filepath (str): Path to the output file.
        file_format (str): Format to save the tree in.
    """
    if not tree:
        print("[Save] No tree provided to save.")
        return
    try:
        Phylo.write(tree, output_filepath, file_format)
        print(f"[Save] Tree successfully saved to: {output_filepath} (format: {file_format})")
    except Exception as e:
        print(f"[Save] Error saving tree to '{output_filepath}': {e}")

# --- Main Execution Example ---
if __name__ == "__main__":
    # Example Newick string (unrooted, F is an arbitrary root for this string representation)
    # ((A,B),(C,D)); is unrooted
    # (A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F; is a rooted Newick string where F is the root.
    # Let's use a clearly unrooted conceptual tree for rooting examples:
    #    --- A
    #   |
    # --X -- B
    #   |
    #    --- Y -- C
    #         |
    #          --- D
    # This can be represented in Newick as (A,B,(C,D));
    # Or with an arbitrary root for string representation: ((A,B)Internal1,(C,D)Internal2)Root;
    
    sample_newick_unrooted = "((TipA:0.1,TipB:0.2):0.05,(TipC:0.3,TipD:0.4):0.08)RootNode;"
    # This Newick string implies RootNode is the root. If we want to treat it as unrooted and then root it:
    # An unrooted tree has a trifurcating root or is drawn without a clear single ancestor.
    # BioPython might interpret Newick like " ( (A,B), (C,D) ); " as rooted at the basal split.
    # Forcing an unrooted interpretation can be tricky just from Newick string,
    # as Newick implies a root. We will assume it's loaded and then we decide the root.

    print("--- Loading Example Tree ---")
    my_tree = load_tree(sample_newick_unrooted, "newick")

    if my_tree:
        print_tree_info(my_tree, "Original Tree")
        display_ascii_tree(my_tree, "Original Tree")

        # --- Rerooting Example ---
        print("\n--- Rerooting Examples ---")
        # 1. Reroot with 'TipA' as outgroup
        #    Note: BioPython's root_with_outgroup modifies the tree in place.
        #    To keep original, you might need to deepcopy: import copy; tree_copy = copy.deepcopy(my_tree)
        reroot_tree(my_tree, outgroup_tip_names=["TipA"]) 
        print_tree_info(my_tree, "Tree Rooted on TipA")
        display_ascii_tree(my_tree, "Tree Rooted on TipA")
        
        # Reload original for next rooting example
        my_tree = load_tree(sample_newick_unrooted, "newick") 
        # 2. Reroot at midpoint (if branch lengths are present and meaningful)
        reroot_tree(my_tree, at_midpoint=True)
        print_tree_info(my_tree, "Tree Rooted at Midpoint")
        display_ascii_tree(my_tree, "Tree Rooted at Midpoint")

        # --- Subtree Extraction Example ---
        print("\n--- Subtree Extraction Example ---")
        # Reload original tree for subtree extraction from a known state
        my_tree = load_tree(sample_newick_unrooted, "newick") 
        
        # Extract subtree for TipC and TipD
        # Their MRCA is the node ancestral to C and D.
        subtree_cd = extract_subtree(my_tree, tip_names_for_subtree=["TipC", "TipD"])
        if subtree_cd:
            print_tree_info(subtree_cd, "Subtree C-D")
            display_ascii_tree(subtree_cd, "Subtree C-D")
            save_tree(subtree_cd, "subtree_CD.nwk", "newick")

        # Extract subtree for TipA, TipB, TipC (MRCA would be the root of the original sample tree)
        subtree_abc = extract_subtree(my_tree, tip_names_for_subtree=["TipA", "TipB", "TipC"])
        if subtree_abc:
            print_tree_info(subtree_abc, "Subtree A-B-C")
            display_ascii_tree(subtree_abc, "Subtree A-B-C")
            save_tree(subtree_abc, "subtree_ABC.nwk", "newick")

        # --- Saving the (last modified) main tree ---
        # Current state of my_tree is midpoint rooted.
        # save_tree(my_tree, "midpoint_rooted_tree.phyloxml", "phyloxml")
        save_tree(my_tree, "midpoint_rooted_tree.nwk", "newick")

    print("\nScript finished.")
    print("Remember to replace sample tree data with your actual tree file/string.")
```

## Explanation and How to Use the Script

1.  **Prerequisites:**
    *   Install Python 3.
    *   Install the Biopython library:
        ```bash
        pip install biopython
        ```

2.  **Loading a Tree (`load_tree`):**
    *   This function can load a tree from a file path or directly from a string (e.g., a Newick string).
    *   **Usage:** `my_tree = load_tree("your_tree_file.nwk", "newick")` or `my_tree = load_tree("(A,B,(C,D));", "newick")`
    *   Specify the `file_format` (e.g., `"newick"`, `"nexus"`, `"phyloxml"`).

3.  **Getting Basic Tree Information (`print_tree_info`):**
    *   Prints the number of tips (terminal nodes), internal nodes, whether the tree is rooted, and total branch length.
    *   **Usage:** `print_tree_info(my_tree, "My Awesome Tree")`

4.  **Displaying an ASCII Tree (`display_ascii_tree`):**
    *   Provides a simple text-based visualization of the tree structure in your console. Useful for quick checks.
    *   **Usage:** `display_ascii_tree(my_tree)`

5.  **Re-rooting a Tree (`reroot_tree`):**
    *   This function allows you to change the root of the tree. This is a common step, especially when starting with an unrooted tree.
    *   **Methods for re-rooting:**
        *   **By Outgroup:** Provide a list of `outgroup_tip_names`. The tree will be rooted on the branch leading to the Most Recent Common Ancestor (MRCA) of these specified outgroup tips. This is biologically the most common way to root a tree if a known outgroup is available.
            *   Usage: `reroot_tree(my_tree, outgroup_tip_names=["OutgroupTaxon1", "OutgroupTaxon2"])`
        *   **At Midpoint:** Set `at_midpoint=True`. This method finds the longest path between any two tips and places the root at the midpoint of that path. It's a common approach when a clear outgroup is not known, but requires meaningful branch lengths.
            *   Usage: `reroot_tree(my_tree, at_midpoint=True)`
    *   **Note:** Re-rooting operations typically modify the tree object in place. If you need to preserve the original tree, make a deep copy first (`import copy; original_tree = copy.deepcopy(my_tree)`).

6.  **Extracting a Subtree (`extract_subtree`):**
    *   This function creates a new tree object containing only the specified `tip_names_for_subtree` and their common ancestors, up to their MRCA.
    *   The new subtree will be rooted at the MRCA of the specified tips.
    *   **Usage:** `my_subtree = extract_subtree(my_tree, tip_names_for_subtree=["SpeciesA", "SpeciesB", "SpeciesC"])`
    *   This is useful for focusing on a specific clade or group of interest within a larger phylogeny.

7.  **Saving a Tree (`save_tree`):**
    *   Saves your (potentially modified) tree object to a new file in the specified format.
    *   **Usage:** `save_tree(my_tree, "modified_tree.nwk", "newick")` or `save_tree(my_subtree, "my_clade.phyloxml", "phyloxml")`

8.  **Main Execution Block (`if __name__ == "__main__":`)**
    *   The script includes an example section that demonstrates loading a sample Newick tree string, printing info, re-rooting it in two ways (outgroup and midpoint), extracting a couple of subtrees, and saving results. You can adapt this section to work with your own tree files.

## Further Possibilities (Beyond this script)

Once you have a tree loaded and manipulated with `Bio.Phylo`, you can perform many other operations:

*   **Traversing the tree:** Iterate over clades, access parent/child relationships (`clade.clades`, `clade.root.get_path(clade)`).
*   **Accessing node attributes:** Get `clade.name`, `clade.branch_length`, `clade.confidence` (bootstrap values).
*   **Comparing trees:** Calculate distances between trees (e.g., Robinson-Foulds distance).
*   **Annotating trees:** Add custom information to clades.
*   **Advanced visualization:** While `draw_ascii` is basic, for publication-quality figures, you'd typically export the tree (e.g., Newick) and use tools like FigTree, iTOL, ETE Toolkit (which also has Python bindings), or R packages like `ggtree`.
*   **Ancestral state reconstruction:** Infer characteristics of ancestral nodes.
*   **Testing phylogenetic hypotheses:** Use statistical tests related to tree topology or branch lengths.

This script provides a foundational toolkit for programmatically interacting with phylogenetic trees in Python. Remember to consult the [BioPython Phylo documentation](https://biopython.org/wiki/Phylo) for more advanced features.

---