#!/usr/bin/env python3

"""
Analyzes a Multiple Sequence Alignment (MSA) file to:
- Calculate a consensus sequence.
- Compute basic statistics.
- Convert to a different MSA format.
"""

import argparse
import sys
import os
from Bio import Align # For Align.parse and Align.write
from Bio.Seq import Seq # Potentially for alphabet
from Bio.Align import MultipleSeqAlignment
from collections import Counter

# --- Helper Function for Alphabet Inference ---
def infer_alphabet_from_msa(alignment, sample_size=5, col_sample_size=20):
    """
    Infers the alphabet type (DNA, RNA, Protein) from a sample of the MSA.
    Returns 'DNA', 'RNA', 'Protein', or 'Unknown'.
    """
    if not alignment or len(alignment) == 0:
        return 'Unknown'

    sequences_to_sample = alignment[:min(sample_size, len(alignment))]
    
    all_chars = Counter()
    for record in sequences_to_sample:
        seq_str = str(record.seq).upper()
        # Sample columns to avoid very long sequences slowing this down
        cols_to_sample = min(len(seq_str), col_sample_size)
        sampled_seq_chars = seq_str[:cols_to_sample] + seq_str[-cols_to_sample:] # Check beginning and end
        
        for char in sampled_seq_chars:
            if char not in ['-', '.']: # Ignore gaps
                all_chars[char] += 1
    
    if not all_chars:
        return 'Unknown' # All gaps or empty

    # Heuristic thresholds
    dna_chars = set(['A', 'C', 'G', 'T'])
    rna_chars = set(['A', 'C', 'G', 'U'])
    
    total_alpha_chars = sum(all_chars.values())
    
    # Check for RNA (presence of U, absence of T)
    if 'U' in all_chars and 'T' not in all_chars:
        non_rna_chars = sum(count for char, count in all_chars.items() if char not in rna_chars)
        if (non_rna_chars / total_alpha_chars) < 0.1: # Allow for some non-standard chars
            return 'RNA'
            
    # Check for DNA (presence of T, absence of U)
    if 'T' in all_chars and 'U' not in all_chars:
        non_dna_chars = sum(count for char, count in all_chars.items() if char not in dna_chars)
        if (non_dna_chars / total_alpha_chars) < 0.1:
            return 'DNA'

    # If ambiguous (e.g., contains both T and U, or mostly other characters)
    # A more sophisticated check might look at amino acid frequencies vs nucleotide frequencies
    # For now, if it's not clearly DNA or RNA, assume Protein or Unknown
    # Count typical protein characters (e.g., L, S, E, K etc.) vs. A,C,G,T,U
    protein_specific_chars = set(all_chars.keys()) - dna_chars - rna_chars - set(['N', 'X']) # N,X can be ambiguous
    if protein_specific_chars: # If any character is clearly not DNA/RNA/N/X
        return 'Protein'

    # If only A,C,G,N,X etc. it's harder to distinguish DNA from protein without more context
    # However, if we only have DNA/RNA chars but couldn't decide above (e.g. only ACG)
    if all(c in dna_chars or c in ['N','X'] for c in all_chars.keys()):
        return 'DNA' # Default to DNA if only DNA-like chars are present

    return 'Protein' # Default assumption if not clearly DNA/RNA


# --- Main Consensus Logic ---
def calculate_consensus(alignment, threshold, ambiguous_char, require_multiple, inferred_alphabet):
    """Calculates and returns the consensus sequence."""
    consensus_sequence = []
    alignment_len = alignment.get_alignment_length()

    # Adjust ambiguous character for DNA/RNA if default 'X' is used
    if inferred_alphabet in ['DNA', 'RNA'] and ambiguous_char == 'X':
        effective_ambiguous_char = 'N'
    else:
        effective_ambiguous_char = ambiguous_char

    for i in range(alignment_len):
        column_str = alignment[:, i]
        counts = Counter(char.upper() for char in column_str)
        
        valid_chars_counts = Counter()
        for char, count in counts.items():
            if char not in ['-', '.']:
                valid_chars_counts[char] += count
        
        if not valid_chars_counts:
            consensus_sequence.append('-') # Or effective_ambiguous_char based on preference
            continue

        total_valid = sum(valid_chars_counts.values())
        
        best_char = None
        best_freq = 0.0
        
        # Find the character with the highest frequency
        sorted_chars = sorted(valid_chars_counts.items(), key=lambda item: item[1], reverse=True)
        
        best_char_candidate, best_char_count = sorted_chars[0]
        best_freq = best_char_count / total_valid

        if best_freq >= threshold:
            # Check the require_multiple condition
            # If the best char meets threshold, but there are other chars present (len > 1),
            # and require_multiple is, for example, 2, meaning we need at least 2 distinct chars
            # for ambiguity, then if only one char type meets threshold, it's the consensus.
            # If require_multiple is 1, this condition is always met for ambiguity if threshold is not met by a single char.
            if len(valid_chars_counts) < require_multiple and require_multiple > 1:
                 # This logic needs refinement based on exact interpretation of require_multiple.
                 # Let's assume: if best_freq >= threshold, it's the consensus_char,
                 # UNLESS len(valid_chars_counts) >= require_multiple, in which case it becomes ambiguous.
                 # This seems counter-intuitive.
                 # A more standard interpretation:
                 # If best_freq >= threshold, use best_char_candidate.
                 # Else (best_freq < threshold), use ambiguous_char.
                 # The require_multiple is typically to force ambiguous if, say, A=60%, T=40% (distinct_chars=2)
                 # and threshold is 50%, but we want ambiguity unless one char is > threshold AND others are very few.
                 # Let's stick to a simpler model first: if threshold met, use char, else ambiguous.
                 # The provided example logic: "if best_freq >= threshold and len(valid_chars) < require_multiple_threshold:"
                 # This implies if there are *few* character types, we are *more* likely to pick the best_char.
                 # If there are *many* character types (>= require_multiple), we are *more* likely to be ambiguous.

                # Simpler interpretation: if best frequency passes threshold, it's the consensus.
                # The 'require_multiple' seems to be intended to force ambiguity if diversity is high,
                # even if one char passes the threshold.
                # Let's re-evaluate the example:
                # if best_freq >= threshold and len(valid_chars) < require_multiple: consensus_char = best_char
                # else: consensus_char = ambiguous_char
                # This means if best_freq hits threshold AND character diversity is LOW, pick best_char.
                # Otherwise (either threshold not met OR diversity is HIGH), pick ambiguous.

                if len(valid_chars_counts) < require_multiple:
                    consensus_sequence.append(best_char_candidate)
                else: # Diversity is high (>= require_multiple), so even if threshold is met by one, call it ambiguous
                    consensus_sequence.append(effective_ambiguous_char)
            else: # Standard case: threshold met, diversity doesn't force ambiguity by itself under this rule
                 consensus_sequence.append(best_char_candidate)

        else: # Threshold not met by any single character
            consensus_sequence.append(effective_ambiguous_char)
            
    return "".join(consensus_sequence)

# --- Main Statistics Logic ---
def calculate_stats(alignment, inferred_alphabet):
    """Calculates and returns a dictionary of MSA statistics."""
    stats = {}
    num_sequences = len(alignment)
    alignment_length = alignment.get_alignment_length()
    stats["Number of sequences"] = num_sequences
    stats["Alignment length"] = alignment_length

    if num_sequences == 0 or alignment_length == 0:
        stats["GC content (%)"] = "N/A"
        stats["Columns entirely gaps (%)"] = "N/A"
        stats["Columns with any gap (%)"] = "N/A"
        return stats

    # GC Content
    if inferred_alphabet in ['DNA', 'RNA']:
        gc_count = 0
        total_bases = 0
        for i in range(alignment_length):
            col_str = alignment[:, i].upper()
            for char in col_str:
                if char in ['G', 'C']:
                    gc_count += 1
                if char in ['A', 'T', 'G', 'C', 'U']: # Include U for RNA
                    total_bases +=1
        stats["GC content (%)"] = f"{(gc_count / total_bases * 100):.2f}" if total_bases > 0 else "N/A"
    else:
        stats["GC content (%)"] = "N/A"

    # Gap Statistics
    cols_all_gaps = 0
    cols_any_gaps = 0
    for i in range(alignment_length):
        col_str = alignment[:, i]
        is_all_gaps = all(char in ['-', '.'] for char in col_str)
        has_any_gaps = any(char in ['-', '.'] for char in col_str)
        if is_all_gaps:
            cols_all_gaps += 1
        if has_any_gaps:
            cols_any_gaps += 1
    
    stats["Columns entirely gaps (%)"] = f"{(cols_all_gaps / alignment_length * 100):.2f}"
    stats["Columns with any gap (%)"] = f"{(cols_any_gaps / alignment_length * 100):.2f}"
    
    return stats

# --- Main Function ---
def main():
    parser = argparse.ArgumentParser(
        description="Analyze a Multiple Sequence Alignment (MSA) file.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "msa_file",
        help="Path to the input MSA file."
    )
    parser.add_argument(
        "--informat",
        required=True,
        help="Format of the input MSA (e.g., 'fasta', 'clustal', 'phylip', 'nexus', 'stockholm')."
    )
    
    # Task flags
    parser.add_argument(
        "--get_consensus",
        action="store_true",
        help="Calculate and output the consensus sequence."
    )
    parser.add_argument(
        "--consensus_threshold",
        type=float,
        default=0.5,
        help="Threshold for a character to be part of the consensus (default: 0.5)."
    )
    parser.add_argument(
        "--consensus_ambiguous_char",
        type=str,
        default='X',
        help="Ambiguous character for consensus if threshold not met (default: 'X')."
    )
    parser.add_argument(
        "--consensus_require_multiple",
        type=int,
        default=1, # Default 1 means any diversity below threshold makes it ambiguous.
                   # If set to e.g. 2, it means if only 1 char type exists (even if > threshold), it's consensus.
                   # If 2+ char types exist AND best is > threshold, it becomes ambiguous. This is the spec's example.
        help="Minimum number of distinct characters in a column to use the ambiguous character, even if one character meets the threshold. (default: 1, meaning this rule is less restrictive by default. Set >1 to make it more likely to use ambiguous_char if diversity is high)."
    )
    
    parser.add_argument(
        "--get_stats",
        action="store_true",
        help="Calculate and output basic MSA statistics."
    )
    
    parser.add_argument(
        "--convert_to",
        type=str,
        metavar="OUTPUT_FORMAT",
        help="Convert MSA to this output format (e.g., 'fasta', 'clustal')."
    )
    parser.add_argument(
        "--outfile",
        type=str,
        metavar="CONVERTED_MSA_FILEPATH",
        help="Output file for converted MSA (required if --convert_to is used)."
    )
    
    parser.add_argument(
        "--output_report",
        type=str,
        metavar="REPORT_FILEPATH",
        help="File to write text reports (consensus, stats) instead of stdout."
    )

    args = parser.parse_args()

    # --- Check if any task is specified ---
    if not (args.get_consensus or args.get_stats or args.convert_to):
        parser.print_help(sys.stderr)
        sys.stderr.write("\nError: No task specified. Please choose at least one action "
                         "(--get_consensus, --get_stats, --convert_to).\n")
        sys.exit(1)

    # --- Read MSA ---
    try:
        alignment_iter = Align.parse(args.msa_file, args.informat)
        alignment = next(alignment_iter)
        # Check if there are more alignments in the file, warn if so
        try:
            next(alignment_iter)
            sys.stderr.write(f"Warning: Multiple alignments found in '{args.msa_file}'. Only processing the first one.\n")
        except StopIteration:
            pass # Expected if only one alignment
            
    except FileNotFoundError:
        sys.stderr.write(f"Error: MSA file '{args.msa_file}' not found.\n")
        sys.exit(1)
    except ValueError as e: # Biopython often raises ValueError for format issues
        sys.stderr.write(f"Error parsing MSA file '{args.msa_file}' with format '{args.informat}'. Details: {e}\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"An unexpected error occurred while reading MSA file '{args.msa_file}': {e}\n")
        sys.exit(1)

    if not alignment: # Should be caught by next(alignment_iter) if file is empty/no alignment
        sys.stderr.write(f"Error: No alignment found in '{args.msa_file}'.\n")
        sys.exit(1)

    # --- Setup Report Stream ---
    report_stream = sys.stdout
    if args.output_report:
        try:
            report_stream = open(args.output_report, 'w')
        except IOError as e:
            sys.stderr.write(f"Error: Could not open report file '{args.output_report}' for writing: {e}\n")
            sys.exit(1)
    
    # --- Infer Alphabet ---
    inferred_alphabet = infer_alphabet_from_msa(alignment)
    if args.get_consensus or args.get_stats: # Only print if relevant
        print(f"Inferred alphabet type: {inferred_alphabet}", file=sys.stderr if report_stream != sys.stdout else report_stream)


    # --- Perform Tasks ---
    if args.get_consensus:
        if report_stream != sys.stdout: sys.stderr.write(f"Calculating consensus sequence...\n")
        consensus_seq = calculate_consensus(
            alignment, 
            args.consensus_threshold, 
            args.consensus_ambiguous_char, 
            args.consensus_require_multiple,
            inferred_alphabet
        )
        report_stream.write(f">Consensus_from_{os.path.basename(args.msa_file)}\n")
        report_stream.write(f"{consensus_seq}\n")
        if report_stream != sys.stdout: sys.stderr.write(f"Consensus sequence written to '{args.output_report}'.\n")

    if args.get_stats:
        if report_stream != sys.stdout: sys.stderr.write(f"Calculating MSA statistics...\n")
        stats = calculate_stats(alignment, inferred_alphabet)
        report_stream.write("\nMSA Statistics:\n")
        for key, value in stats.items():
            report_stream.write(f"  {key}: {value}\n")
        if report_stream != sys.stdout: sys.stderr.write(f"MSA statistics written to '{args.output_report}'.\n")

    if args.convert_to:
        if not args.outfile:
            sys.stderr.write("Error: --outfile is required when --convert_to is specified.\n")
            if report_stream != sys.stdout: report_stream.close()
            sys.exit(1)
        
        sys.stderr.write(f"Converting MSA to '{args.convert_to}' format, saving to '{args.outfile}'...\n")
        try:
            # Bio.Align.write expects an iterable of alignments
            Align.write([alignment], args.outfile, args.convert_to)
            sys.stderr.write(f"Successfully converted MSA to '{args.outfile}'.\n")
        except Exception as e:
            sys.stderr.write(f"Error during MSA conversion to format '{args.convert_to}': {e}\n")
            if report_stream != sys.stdout: report_stream.close()
            sys.exit(1)

    # --- Close Report File if Opened ---
    if report_stream != sys.stdout:
        report_stream.close()

if __name__ == '__main__':
    main()
```
