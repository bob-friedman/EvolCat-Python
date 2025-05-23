#!/usr/bin/env python3

import argparse
import sys
import math
import itertools
from pylib.utils import seq_parser

def main():
    """
    Calculates Kimura 2-Parameter (K2P) distances, Ts/Tv ratios, and standard error
    for all unique pairs of sequences in a FASTA file.
    """
    parser = argparse.ArgumentParser(
        description="Calculate K2P distances between all unique pairs of sequences in a FASTA file."
    )
    parser.add_argument(
        "input_fasta_file",
        help="Path to the input FASTA file (aligned sequences)."
    )
    args = parser.parse_args()

    # Print header
    print("Gene1\tGene2\tK2P_Distance\tSE_K2P\tTs_Tv_Ratio")

    records = []
    try:
        for record in seq_parser.parse_fasta_file(args.input_fasta_file):
            records.append(record)
    except FileNotFoundError:
        print(f"Error: Input file not found at {args.input_fasta_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading FASTA file: {e}", file=sys.stderr)
        sys.exit(1)

    # Validate Input
    if len(records) < 2:
        print("Error: At least two sequences are required for K2P calculation.", file=sys.stderr)
        sys.exit(1)

    alignment_length = len(records[0].seq)
    for i, record in enumerate(records):
        if len(record.seq) != alignment_length:
            print(f"Error: Sequences must be aligned and of equal length. "
                  f"Sequence '{records[0].id}' has length {alignment_length}, "
                  f"but sequence '{record.id}' (index {i}) has length {len(record.seq)}.",
                  file=sys.stderr)
            sys.exit(1)
    
    if alignment_length == 0:
        print("Error: Sequences are empty. Cannot calculate K2P distance.", file=sys.stderr)
        sys.exit(1)

    # Define purines and pyrimidines
    purines = {'A', 'G'}
    pyrimidines = {'C', 'T'}
    valid_bases = purines.union(pyrimidines)

    # Iterate through all unique pairs of sequences
    for seq1_record, seq2_record in itertools.combinations(records, 2):
        s1 = str(seq1_record.seq).upper()
        s2 = str(seq2_record.seq).upper()

        transitions = 0 # P (total transitions)
        transversions = 0 # Q (total transversions)
        valid_sites = 0

        for i in range(alignment_length):
            char1 = s1[i]
            char2 = s2[i]

            if char1 not in valid_bases or char2 not in valid_bases:
                continue # Skip site if either base is not valid DNA (A, T, C, G)
            
            valid_sites += 1

            if char1 != char2:
                is_char1_purine = char1 in purines
                is_char2_purine = char2 in purines
                
                if (is_char1_purine and is_char2_purine) or \
                   (not is_char1_purine and not is_char2_purine): # Both pyrimidines
                    transitions += 1
                else:
                    transversions += 1
        
        K_str, SE_K_str, R_str = "N/A", "N/A", "N/A"

        if valid_sites == 0:
            print(f"Warning: No valid comparable sites found between {seq1_record.id} and {seq2_record.id}.", file=sys.stderr)
        else:
            P_prop = transitions / valid_sites
            Q_prop = transversions / valid_sites

            val_for_log1 = 1 - 2*P_prop - Q_prop
            val_for_log2 = 1 - 2*Q_prop

            if val_for_log1 > 0 and val_for_log2 > 0:
                try:
                    K = -0.5 * math.log(val_for_log1) - 0.25 * math.log(val_for_log2)
                    K_str = f"{K:.5f}"

                    # Ts/Tv Ratio (R)
                    # The Perl script calculates s_val and v_val which are components of K for R.
                    # s_val = -0.5 * log(1-2P-Q)
                    # v_val = -0.25 * log(1-2Q)
                    # R = s_val / (2*v_val) effectively if P & Q are proportions.
                    # Or, more directly, R = P/Q if Q_prop > 0.
                    # The Perl script's R calculation:
                    # S = 0.5 * log(w1) - 0.25 * log(w2) -> -0.5*log(1-2P-Q) - (-0.25*log(1-2Q))
                    # V = 0.5 * log(w2) -> -0.5*log(1-2Q)
                    # R = S/V. This is (Ts rate) / (Tv rate part)
                    # For K2P, a common Ts/Tv ratio is simply transitions/transversions counts.
                    # Let's use the ratio of proportions P_prop / Q_prop
                    K = -0.5 * math.log(val_for_log1) - 0.25 * math.log(val_for_log2)
                    K_str = f"{K:.5f}"

                    # Ts/Tv Ratio (R) - based on Perl script's S/V components
                    # S_perl = -0.5*math.log(1-2P-Q) + 0.25*math.log(1-2Q)
                    # V_perl = -0.5*math.log(1-2Q)
                    s_val_perl = -0.5 * math.log(val_for_log1) + 0.25 * math.log(val_for_log2)
                    v_val_perl = -0.5 * math.log(val_for_log2)
                    
                    if v_val_perl != 0:
                        R = s_val_perl / v_val_perl
                        R_str = f"{R:.5f}"
                    else:
                        R_str = "N/A" # Or "Inf" if s_val_perl is non-zero

                    # Standard Error (SE_K) - using formula from prompt, derived from Perl
                    w1_perl = 1 / val_for_log1 
                    w2_perl = 1 / val_for_log2
                    w3_perl = 0.5 * (w1_perl + w2_perl)
                    
                    term_P_variance_like = (w1_perl**2 * P_prop)
                    term_Q_variance_like = (w3_perl**2 * Q_prop)
                    covariance_like_term = ( (w1_perl*P_prop) + (w3_perl*Q_prop) )**2
                    
                    # Ensure valid_sites is not zero before division
                    if valid_sites > 0:
                        variance_K = (term_P_variance_like + term_Q_variance_like - covariance_like_term) / valid_sites
                        if variance_K >= 0:
                            SE_K = math.sqrt(variance_K)
                            SE_K_str = f"{SE_K:.5f}"
                        else:
                            SE_K_str = "N/A" # Variance is negative
                    else: # Should not happen if we already checked valid_sites > 0 for K calculation
                        SE_K_str = "N/A"
                        
                except (ValueError, ZeroDivisionError, OverflowError) as e:
                    # Catch math errors if logs or divisions fail
                    print(f"Math error during K2P calculation for {seq1_record.id} and {seq2_record.id}: {e}", file=sys.stderr)
                    # K_str, SE_K_str, R_str remain "N/A"
        
        print(f"{seq1_record.id}\t{seq2_record.id}\t{K_str}\t{SE_K_str}\t{R_str}")

if __name__ == "__main__":
    main()
