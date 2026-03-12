#!/usr/bin/env python3
"""
pal2nal_enforce.py

Enforces an amino acid sequence alignment (FASTA) onto the corresponding
unaligned nucleotide sequences (FASTA), inserting '---' for every gap ('-')
in the amino acid alignment. Output is a codon-aware nucleotide alignment.

Usage:
    python pal2nal_enforce.py -a <aligned_aa.fasta> -n <unaligned_nt.fasta> -o <output_nt_aligned.fasta>

Arguments:
    -a / --aa_aligned     Aligned amino acid sequences in FASTA format
    -n / --nt_unaligned   Unaligned nucleotide sequences in FASTA format
    -o / --output         Output file for aligned nucleotide sequences (default: stdout)
    --gap_char            Gap character to use in output (default: -)
"""

import argparse
import sys
from pathlib import Path


# ----------------------------------------------------------------------
# FASTA I/O helpers
# ----------------------------------------------------------------------

def read_fasta(filepath):
    """
    Parse a FASTA file and return an ordered list of (header, sequence) tuples.
    The header is stored without the leading '>'.
    """
    records = []
    header = None
    seq_parts = []

    with open(filepath, "r") as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header = line[1:]          # strip '>'
                seq_parts = []
            else:
                seq_parts.append(line)

    if header is not None:
        records.append((header, "".join(seq_parts)))

    return records


def write_fasta(records, filehandle, line_width=60):
    """
    Write a list of (header, sequence) tuples to a file handle in FASTA format.
    Sequences are wrapped at line_width characters.
    """
    for header, seq in records:
        filehandle.write(f">{header}\n")
        for i in range(0, len(seq), line_width):
            filehandle.write(seq[i : i + line_width] + "\n")


# ----------------------------------------------------------------------
# Core logic
# ----------------------------------------------------------------------

def extract_seq_id(header):
    """
    Return the first whitespace-delimited token of a FASTA header.
    This is used to match amino acid and nucleotide records by ID.
    """
    return header.split()[0]


def enforce_aa_alignment_on_nt(aa_aligned_seq, nt_seq, gap_char="-"):
    """
    Given an aligned amino acid sequence and its corresponding unaligned
    nucleotide sequence, return the nucleotide sequence with gaps inserted
    to match the amino acid alignment.

    Each amino acid position maps to exactly 3 nucleotide positions (a codon).
    Each gap '-' in the amino acid alignment becomes '---' in the nucleotide output.

    Raises ValueError if the ungapped amino acid length * 3 does not match
    the nucleotide sequence length.
    """
    ungapped_aa_len = len(aa_aligned_seq.replace(gap_char, ""))
    expected_nt_len = ungapped_aa_len * 3

    if len(nt_seq) != expected_nt_len:
        raise ValueError(
            f"Length mismatch: ungapped AA length ({ungapped_aa_len}) × 3 = "
            f"{expected_nt_len}, but nucleotide sequence has {len(nt_seq)} bases."
        )

    nt_aligned = []
    nt_pos = 0

    for aa_char in aa_aligned_seq:
        if aa_char == gap_char:
            nt_aligned.append("---")
        else:
            codon = nt_seq[nt_pos : nt_pos + 3]
            nt_aligned.append(codon)
            nt_pos += 3

    return "".join(nt_aligned)


def build_nt_alignment(aa_records, nt_records, gap_char="-"):
    """
    Match amino acid and nucleotide records by sequence ID and enforce
    the amino acid alignment on each nucleotide sequence.

    Returns an ordered list of (header, aligned_nt_sequence) tuples,
    preserving the order of the amino acid alignment file.
    """
    # Index nucleotide records by their sequence ID for fast lookup
    nt_index = {extract_seq_id(hdr): (hdr, seq) for hdr, seq in nt_records}

    aligned_nt_records = []
    warnings = []

    for aa_header, aa_seq in aa_records:
        seq_id = extract_seq_id(aa_header)

        if seq_id not in nt_index:
            warnings.append(
                f"  WARNING: '{seq_id}' found in AA alignment but not in NT file — skipping."
            )
            continue

        nt_header, nt_seq = nt_index[seq_id]

        try:
            aligned_nt_seq = enforce_aa_alignment_on_nt(aa_seq, nt_seq.upper(), gap_char)
        except ValueError as err:
            warnings.append(f"  ERROR for '{seq_id}': {err} — skipping.")
            continue

        # Use the nucleotide header to preserve any original annotation
        aligned_nt_records.append((nt_header, aligned_nt_seq))

    return aligned_nt_records, warnings


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Enforce an amino acid alignment onto unaligned nucleotide sequences. "
            "Gaps in the AA alignment are expanded to '---' (one codon gap) in the output."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "-a", "--aa_aligned",
        required=True,
        metavar="FILE",
        help="Aligned amino acid sequences in FASTA format.",
    )
    parser.add_argument(
        "-n", "--nt_unaligned",
        required=True,
        metavar="FILE",
        help="Unaligned nucleotide sequences in FASTA format.",
    )
    parser.add_argument(
        "-o", "--output",
        metavar="FILE",
        default=None,
        help="Output file for the aligned nucleotide sequences (default: stdout).",
    )
    parser.add_argument(
        "--gap_char",
        default="-",
        metavar="CHAR",
        help="Gap character used in the amino acid alignment (default: -).",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # --- Validate inputs ---
    for label, path in [("AA alignment", args.aa_aligned), ("NT sequences", args.nt_unaligned)]:
        if not Path(path).is_file():
            sys.exit(f"Error: {label} file not found: '{path}'")

    # --- Read inputs ---
    print(f"Reading AA alignment:      {args.aa_aligned}", file=sys.stderr)
    aa_records = read_fasta(args.aa_aligned)

    print(f"Reading NT sequences:      {args.nt_unaligned}", file=sys.stderr)
    nt_records = read_fasta(args.nt_unaligned)

    print(f"AA sequences in alignment: {len(aa_records)}", file=sys.stderr)
    print(f"NT sequences provided:     {len(nt_records)}", file=sys.stderr)

    # --- Build alignment ---
    aligned_nt_records, warnings = build_nt_alignment(aa_records, nt_records, args.gap_char)

    # --- Report any issues ---
    if warnings:
        print("\nIssues encountered:", file=sys.stderr)
        for w in warnings:
            print(w, file=sys.stderr)

    print(f"\nSuccessfully aligned:      {len(aligned_nt_records)} sequences", file=sys.stderr)

    # --- Write output ---
    if args.output:
        with open(args.output, "w") as out_fh:
            write_fasta(aligned_nt_records, out_fh)
        print(f"Output written to:         {args.output}", file=sys.stderr)
    else:
        write_fasta(aligned_nt_records, sys.stdout)


if __name__ == "__main__":
    main()
