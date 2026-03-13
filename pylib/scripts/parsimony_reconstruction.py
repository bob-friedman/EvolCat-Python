"""
parsimony_reconstruction.py — Ancestral Sequence Reconstruction by Fitch Parsimony
====================================================================================
Parses a Newick-formatted tree file with mutations annotated on branches and a
reference sequence at the root node.  Reconstructs biological sequences at all
internal nodes and leaves using a two-pass algorithm.

Algorithm (two passes)
----------------------
Pass 1 — root → tips (calculate_initial_sequences):
    Propagates the root/reference sequence down every branch, applying the
    mutations annotated in the Newick string at each node.  This produces a
    first estimate of every node's sequence from the known substitution history
    (the format produced by UShER / matUtils).

Pass 2 — tips → root (reconstruct_ancestral_sequences_parsimony):
    Applies Fitch's (1971) parsimony algorithm bottom-up.  At each internal
    node and for each sequence site:
      • If the intersection of children's IUPAC base sets is non-empty, the
        intersection is assigned (parsimony score unchanged).
      • If the intersection is empty, the union is assigned and encoded as the
        appropriate IUPAC ambiguity code (parsimony score +1).
    The two-pass design is the correct approach for UShER-derived trees, where
    branch mutations are already recorded: Pass 1 provides informed leaf
    sequences; Pass 2 refines the internal nodes without re-inventing the
    mutation history.

Pipeline context
----------------
This script implements the Ancestral State Reconstruction step of the
sequence-to-sequence transformer training pipeline described in:

    EvolCat-Python/guides/PIPELINE.md

In that pipeline, ancestral sequences reconstructed here serve as the source
("ancestor") sequences in (ancestor, descendant) training pairs for the seq2seq
transformer model.  IQ-TREE maximum-likelihood ASR is the higher-accuracy
alternative for the same role; this script provides a dependency-free baseline
that runs on any annotated Newick tree without external tools.

Use --fasta to write all node sequences to a FASTA file suitable for direct
input to the transformer training pipeline.

See also:
    EvolCat-Python/guides/PIPELINE_2.md  — Spike-focused interpretable model
    data-sarscov2-genomes/scripts/       — Origin repository and SARS-CoV-2
                                           application context

Input
-----
A two-line text file:
    Line 1: Newick tree string with mutations annotated on branches.
            Example: ((Leaf1:[1C>A,3T>A],Leaf2:[1C>T])Anc1)Root;
            Mutations are 0-indexed: position, reference base (optional),
            '>', new base.  E.g., "1C>A" = position 1, C mutates to A.
    Line 2: Reference / root sequence string (A/C/G/T + IUPAC codes).
            Example: CGTGA

Output
------
Default (stdout): human-readable tree with node names and inferred sequences.
--fasta <file>  : FASTA file — one record per node, suitable for pipeline use.
--pairs         : emit only (internal-node, child) FASTA pairs for transformer
                  training — writes ancestor.fas and descendant.fas.

Usage
-----
    python parsimony_reconstruction.py <input_file>
    python parsimony_reconstruction.py <input_file> --fasta all_nodes.fas
    python parsimony_reconstruction.py <input_file> --pairs

References
----------
Fitch WM (1971) Toward defining the course of evolution: minimum change for
    a specific tree topology.  Syst Zool 20:406-416.
Turakhia Y et al. (2021) Ultrafast Sample placement on Existing tRees (UShER).
    Nat Genet 53:809-816.

Authors
-------
    Bob Friedman
    Claude (Anthropic) — migration, --fasta/--pairs output modes, validation
"""

import re
import sys
import argparse
from typing import List, Dict, Set, Optional

# ---------------------------------------------------------------------------
# IUPAC data structures
# ---------------------------------------------------------------------------

IUPAC_EXPANSION: Dict[str, Set[str]] = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'}, 'R': {'A', 'G'},
    'Y': {'C', 'T'}, 'M': {'A', 'C'}, 'K': {'G', 'T'}, 'W': {'A', 'T'},
    'S': {'C', 'G'}, 'H': {'A', 'C', 'T'}, 'D': {'A', 'G', 'T'},
    'V': {'A', 'C', 'G'}, 'B': {'C', 'G', 'T'}, 'N': {'A', 'C', 'G', 'T'}
}

IUPAC_MAP: Dict[str, str] = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'AG': 'R', 'CT': 'Y',
    'AC': 'M', 'GT': 'K', 'AT': 'W', 'CG': 'S', 'ACT': 'H',
    'AGT': 'D', 'ACG': 'V', 'CGT': 'B', 'ACGT': 'N'
}

# ---------------------------------------------------------------------------
# Tree node
# ---------------------------------------------------------------------------

class PhyloNode:
    """Represents a node in a phylogenetic tree."""
    def __init__(self, name: str, branch_info: Optional[str] = None):
        self.name: str = name
        self.branch_info: Optional[str] = branch_info
        self.children: List['PhyloNode'] = []
        self.inferred_sequence: List[str] = []

    def __repr__(self) -> str:
        return f"PhyloNode(name='{self.name}')"

    def is_leaf(self) -> bool:
        return len(self.children) == 0

# ---------------------------------------------------------------------------
# Newick parser
# ---------------------------------------------------------------------------

def _split_top_level(s: str, delimiter: str) -> List[str]:
    """Split by delimiter only at depth 0 (not inside parentheses or brackets)."""
    parts: List[str] = []
    buf: List[str] = []
    p_depth = b_depth = 0
    for ch in s:
        if   ch == '(': p_depth += 1
        elif ch == ')': p_depth -= 1
        elif ch == '[': b_depth += 1
        elif ch == ']': b_depth -= 1
        if ch == delimiter and p_depth == 0 and b_depth == 0:
            parts.append(''.join(buf)); buf = []
        else:
            buf.append(ch)
    parts.append(''.join(buf))
    return [p.strip() for p in parts if p.strip()]


def parse_newick(s: str) -> PhyloNode:
    """
    Parse a Newick string (with optional branch mutation annotations) into a
    PhyloNode tree.

    Example:
        ((Leaf1:[1C>A],Leaf2:[3T>G])Anc1)Root;
    """
    s = s.strip().rstrip(';')

    def _recursive_parse(substring: str) -> PhyloNode:
        substring = substring.strip()
        if not substring.startswith('('):
            parts = substring.split(':', 1)
            return PhyloNode(
                name=parts[0],
                branch_info=parts[1] if len(parts) > 1 else None
            )
        depth = 0
        match_idx = -1
        for i, ch in enumerate(substring):
            if   ch == '(': depth += 1
            elif ch == ')':
                depth -= 1
                if depth == 0: match_idx = i; break
        if match_idx == -1:
            raise ValueError(f"Unmatched parentheses in Newick substring: '{substring}'")
        children_str    = substring[1:match_idx]
        parent_info_str = substring[match_idx + 1:]
        parent_parts    = parent_info_str.split(':', 1)
        node = PhyloNode(
            name=parent_parts[0],
            branch_info=parent_parts[1] if len(parent_parts) > 1 else None
        )
        for cs in _split_top_level(children_str, ','):
            node.children.append(_recursive_parse(cs))
        return node

    return _recursive_parse(s)

# ---------------------------------------------------------------------------
# Pass 1 — root → tips: apply branch mutations
# ---------------------------------------------------------------------------

def apply_mutations(base_seq: List[str], branch_info: Optional[str]) -> List[str]:
    """
    Apply mutations encoded in a branch annotation string to a copy of base_seq.

    Mutation format: [POS(REF)>NEW, ...]
        POS  — 0-indexed position (integer)
        REF  — reference nucleotide (optional, for human readability)
        NEW  — new nucleotide (single IUPAC character)

    Examples: "1C>A"  "0>G"  "[2T>C,10A>G]"
    """
    new_seq = list(base_seq)
    if not branch_info:
        return new_seq
    for pos_str, new_nt in re.compile(r'(\d+)[A-Z]?>([A-Z])').findall(branch_info):
        pos = int(pos_str)
        if 0 <= pos < len(new_seq):
            new_seq[pos] = new_nt
        else:
            print(
                f"Warning: mutation at position {pos} is out of bounds for sequence "
                f"of length {len(new_seq)} — skipped.",
                file=sys.stderr
            )
    return new_seq


def calculate_initial_sequences(node: PhyloNode, parent_seq: List[str]) -> None:
    """
    Pass 1 (root → tips): propagate the root sequence down the tree, applying
    annotated mutations at each branch.
    """
    node.inferred_sequence = apply_mutations(parent_seq, node.branch_info)
    for child in node.children:
        calculate_initial_sequences(child, node.inferred_sequence)

# ---------------------------------------------------------------------------
# Pass 2 — tips → root: Fitch parsimony refinement
# ---------------------------------------------------------------------------

def _get_expanded_bases(base_char: str, node_name: str, site_idx: int) -> Set[str]:
    """Expand a (possibly IUPAC-ambiguous) base to its set of nucleotides."""
    if base_char not in IUPAC_EXPANSION:
        print(
            f"Warning: unknown base '{base_char}' at site {site_idx} of node "
            f"'{node_name}' — treated as 'N'.",
            file=sys.stderr
        )
        return set(IUPAC_EXPANSION['N'])
    return set(IUPAC_EXPANSION[base_char])


def reconstruct_ancestral_sequences_parsimony(node: PhyloNode) -> None:
    """
    Pass 2 (tips → root): Fitch (1971) parsimony algorithm.

    For each site at each internal node:
      - Intersection of children's base sets non-empty → assign intersection.
      - Intersection empty → assign union (parsimony cost +1), encode as IUPAC.
    """
    if node.is_leaf():
        return
    for child in node.children:
        reconstruct_ancestral_sequences_parsimony(child)

    site_count = len(node.children[0].inferred_sequence)
    if not all(len(c.inferred_sequence) == site_count for c in node.children):
        raise ValueError(
            f"Children of node '{node.name}' have inconsistent sequence lengths."
        )
    if site_count == 0:
        node.inferred_sequence = []
        return

    new_seq: List[str] = []
    for i in range(site_count):
        # Start with the first child's base set
        current = _get_expanded_bases(
            node.children[0].inferred_sequence[i],
            node.children[0].name, i
        )
        # Intersect with each remaining child
        for child in node.children[1:]:
            current.intersection_update(
                _get_expanded_bases(child.inferred_sequence[i], child.name, i)
            )
        # If intersection is empty, take union (Fitch rule)
        if not current:
            current = set()
            for child in node.children:
                current.update(
                    _get_expanded_bases(child.inferred_sequence[i], child.name, i)
                )
        new_seq.append(IUPAC_MAP.get(''.join(sorted(current)), 'N'))
    node.inferred_sequence = new_seq

# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def collect_all_nodes(node: PhyloNode, result: Optional[List[PhyloNode]] = None
                      ) -> List[PhyloNode]:
    """Pre-order traversal — returns all nodes in the tree."""
    if result is None:
        result = []
    result.append(node)
    for child in node.children:
        collect_all_nodes(child, result)
    return result


def print_tree_sequences(node: PhyloNode, depth: int = 0) -> None:
    """Recursively print each node's name and inferred sequence."""
    indent    = "  " * depth
    node_type = "Internal" if node.children else "Leaf"
    seq_str   = ''.join(node.inferred_sequence) if node.inferred_sequence else "N/A"
    print(f"{indent}- {node_type} Node: {node.name}")
    print(f"{indent}  - Final Sequence: {seq_str}")
    for child in node.children:
        print_tree_sequences(child, depth + 1)


def write_fasta(nodes: List[PhyloNode], path: str) -> None:
    """Write all node sequences to a FASTA file."""
    with open(path, 'w') as fh:
        for node in nodes:
            seq = ''.join(node.inferred_sequence)
            if seq:
                fh.write(f">{node.name}\n{seq}\n")
    print(f"[INFO] FASTA written: {path} ({len(nodes)} sequences)", file=sys.stderr)


def write_pairs(root: PhyloNode, anc_path: str, desc_path: str) -> None:
    """
    Write (internal-node, child) sequence pairs to two FASTA files — the
    format expected by the seq2seq transformer training pipeline in
    EvolCat-Python/guides/PIPELINE.md.

    ancestor.fas    — one record per parent–child edge, labelled ParentName
    descendant.fas  — matching record per edge, labelled ChildName
    """
    anc_records:  List[str] = []
    desc_records: List[str] = []

    def _collect_pairs(node: PhyloNode) -> None:
        parent_seq = ''.join(node.inferred_sequence)
        for child in node.children:
            child_seq = ''.join(child.inferred_sequence)
            if parent_seq and child_seq:
                anc_records.append(f">{node.name}\n{parent_seq}")
                desc_records.append(f">{child.name}\n{child_seq}")
            _collect_pairs(child)

    _collect_pairs(root)

    with open(anc_path, 'w') as fa:
        fa.write('\n'.join(anc_records) + '\n')
    with open(desc_path, 'w') as fd:
        fd.write('\n'.join(desc_records) + '\n')

    n_pairs = len(anc_records)
    print(
        f"[INFO] {n_pairs} training pair(s) written:\n"
        f"       ancestors   → {anc_path}\n"
        f"       descendants → {desc_path}",
        file=sys.stderr
    )

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Reconstruct ancestral sequences using a two-pass parsimony method "
            "(Fitch 1971) on a Newick tree with branch-annotated mutations."
        )
    )
    parser.add_argument(
        "input_file_path",
        type=str,
        help=(
            "Path to input file.  Line 1: Newick string with branch mutations "
            "(e.g. ((Leaf1:[1C>A],Leaf2)Anc1)Root;).  "
            "Line 2: root/reference sequence."
        )
    )
    parser.add_argument(
        "--fasta",
        metavar="FILE",
        default=None,
        help="Write all node sequences to a FASTA file instead of printing."
    )
    parser.add_argument(
        "--pairs",
        action="store_true",
        help=(
            "Write (ancestor, descendant) FASTA pairs for transformer training. "
            "Produces ancestor.fas and descendant.fas in the current directory."
        )
    )
    args = parser.parse_args()

    # --- Read input ---
    try:
        with open(args.input_file_path, 'r') as fh:
            newick_string   = fh.readline().strip()
            root_seq_str    = fh.readline().strip().upper()
        if not newick_string:
            raise ValueError("Line 1 (Newick string) is missing or empty.")
        if not root_seq_str:
            raise ValueError("Line 2 (root sequence) is missing or empty.")
        if not re.match(r'^[ACGTNRYMKSWBDHV]+$', root_seq_str):
            raise ValueError(
                f"Invalid characters in root sequence: '{root_seq_str}'. "
                "Only standard IUPAC nucleotide codes are permitted."
            )
    except FileNotFoundError:
        print(f"[ERROR] File not found: '{args.input_file_path}'", file=sys.stderr)
        sys.exit(1)
    except ValueError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        sys.exit(1)

    print(f"Input file   : {args.input_file_path}", file=sys.stderr)
    print(f"Newick       : {newick_string}",        file=sys.stderr)
    print(f"Root sequence: {root_seq_str}\n",       file=sys.stderr)

    # --- Parse and reconstruct ---
    try:
        tree_root = parse_newick(newick_string)
    except ValueError as exc:
        print(f"[ERROR] Newick parse failed: {exc}", file=sys.stderr)
        sys.exit(1)

    # Pass 1: root → tips, applying annotated mutations
    calculate_initial_sequences(tree_root, list(root_seq_str))

    # Pass 2: tips → root, Fitch parsimony refinement
    reconstruct_ancestral_sequences_parsimony(tree_root)

    # --- Output ---
    all_nodes = collect_all_nodes(tree_root)

    if args.pairs:
        write_pairs(tree_root, "ancestor.fas", "descendant.fas")
    elif args.fasta:
        write_fasta(all_nodes, args.fasta)
    else:
        print("--- Final Inferred Sequences for All Nodes (with Parsimony) ---")
        print_tree_sequences(tree_root)


if __name__ == '__main__':
    main()
