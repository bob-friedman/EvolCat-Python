#!/usr/bin/env python3

import argparse
import sys
import math
import itertools
from pylib.utils import seq_parser

# --- Helper functions for Distance Calculations ---

def calculate_jc69(P_prop, Q_prop, valid_sites):
    """Calculates Jukes-Cantor distance and its SE."""
    D_prop = P_prop + Q_prop
    jc_dist_str, jc_se_str = "N/A", "N/A"

    if D_prop >= 0.75 or (1.0 - (4.0/3.0) * D_prop) <= 1e-9: # Added tolerance for float comparison
        # Condition where log argument would be <= 0 or distance is undefined
        pass
    else:
        try:
            jc_dist = -0.75 * math.log(1.0 - (4.0/3.0) * D_prop)
            jc_dist_str = f"{jc_dist:.5f}"
            
            # SE calculation
            variance_denom = (1.0 - (4.0/3.0) * D_prop)**2 * valid_sites
            if variance_denom > 0:
                jc_se_variance_val = (D_prop * (1.0 - D_prop)) / variance_denom
                if jc_se_variance_val >=0: # ensure variance is not negative due to float issues
                     jc_se = math.sqrt(jc_se_variance_val)
                     jc_se_str = f"{jc_se:.5f}"
        except (ValueError, ZeroDivisionError):
            pass # Keep as N/A
    return jc_dist_str, jc_se_str

def calculate_k2p(P_prop, Q_prop, valid_sites):
    """Calculates Kimura 2-Parameter distance, SE, and Ts/Tv ratio."""
    k2p_dist_str, k2p_se_str, ts_tv_ratio_str = "N/A", "N/A", "N/A"

    val_for_log1_k2p = 1.0 - 2*P_prop - Q_prop
    val_for_log2_k2p = 1.0 - 2*Q_prop

    if val_for_log1_k2p > 1e-9 and val_for_log2_k2p > 1e-9: # Added tolerance
        try:
            k2p_dist = -0.5 * math.log(val_for_log1_k2p) - 0.25 * math.log(val_for_log2_k2p)
            k2p_dist_str = f"{k2p_dist:.5f}"

            # Ts/Tv Ratio (R) - based on Perl script's S/V components
            s_val_perl = -0.5 * math.log(val_for_log1_k2p) + 0.25 * math.log(val_for_log2_k2p)
            v_val_perl = -0.5 * math.log(val_for_log2_k2p)
            
            if abs(v_val_perl) > 1e-9: # Check for non-zero denominator with tolerance
                R = s_val_perl / v_val_perl
                ts_tv_ratio_str = f"{R:.5f}"
            elif P_prop == 0 and Q_prop == 0: # Both zero, Ts/Tv is effectively 0
                 ts_tv_ratio_str = "0.00000"
            elif P_prop > 0 and Q_prop == 0: # Transitions but no transversions
                 ts_tv_ratio_str = "Inf"
            else: # v_val_perl is near zero, Q_prop might be non-zero but leads to undefined R
                 ts_tv_ratio_str = "N/A"


            # Standard Error (SE_K) - using formula from prompt, derived from Perl
            w1_perl = 1 / val_for_log1_k2p
            w3_perl = 0.5 * (w1_perl + (1/val_for_log2_k2p))
            
            term_P_variance_like = (w1_perl**2 * P_prop)
            term_Q_variance_like = (w3_perl**2 * Q_prop)
            covariance_like_term = ( (w1_perl*P_prop) + (w3_perl*Q_prop) )**2
            
            if valid_sites > 0:
                variance_K = (term_P_variance_like + term_Q_variance_like - covariance_like_term) / valid_sites
                if variance_K >= 0:
                    SE_K = math.sqrt(variance_K)
                    k2p_se_str = f"{SE_K:.5f}"
        except (ValueError, ZeroDivisionError, OverflowError):
            pass # Keep as N/A
    return k2p_dist_str, k2p_se_str, ts_tv_ratio_str


def calculate_tajima_nei(pair_counts, valid_sites, P_prop, Q_prop):
    """Calculates Tajima-Nei distance and its SE."""
    D_prop = P_prop + Q_prop
    dist_str, se_str = "N/A", "N/A"

    if valid_sites == 0: return dist_str, se_str

    # Calculate base frequencies from pair counts
    f = {k: v / valid_sites for k, v in pair_counts.items()}
    
    qA = f['AA'] + (f['AT'] + f['AG'] + f['AC']) / 2.0
    qT = f['TT'] + (f['AT'] + f['TG'] + f['TC']) / 2.0
    qG = f['GG'] + (f['AG'] + f['TG'] + f['GC']) / 2.0
    qC = f['CC'] + (f['AC'] + f['TC'] + f['GC']) / 2.0
    
    # Ensure sum of freqs is close to 1 (can have minor float issues)
    # Normalization might be needed if sums deviate significantly
    
    qI2 = qA**2 + qT**2 + qG**2 + qC**2
    
    sum_of_fXY_sq_terms_over_2qXqY = 0
    # Terms like (fAT**2)/(2*qA*qT)
    if qA > 1e-9 and qT > 1e-9: sum_of_fXY_sq_terms_over_2qXqY += (f['AT']**2) / (2*qA*qT)
    if qA > 1e-9 and qG > 1e-9: sum_of_fXY_sq_terms_over_2qXqY += (f['AG']**2) / (2*qA*qG)
    if qA > 1e-9 and qC > 1e-9: sum_of_fXY_sq_terms_over_2qXqY += (f['AC']**2) / (2*qA*qC)
    if qT > 1e-9 and qG > 1e-9: sum_of_fXY_sq_terms_over_2qXqY += (f['TG']**2) / (2*qT*qG)
    if qT > 1e-9 and qC > 1e-9: sum_of_fXY_sq_terms_over_2qXqY += (f['TC']**2) / (2*qT*qC)
    if qG > 1e-9 and qC > 1e-9: sum_of_fXY_sq_terms_over_2qXqY += (f['GC']**2) / (2*qG*qC)

    c_tajima = sum_of_fXY_sq_terms_over_2qXqY

    if c_tajima > 1e-9: # c_tajima must be positive
        try:
            b_tajima = 0.5 * (1.0 - qI2 + (D_prop**2 / c_tajima))
            log_arg = 1.0 - (D_prop / b_tajima) if abs(b_tajima) > 1e-9 else -1.0 # Avoid division by zero

            if log_arg > 1e-9:
                dist = -b_tajima * math.log(log_arg)
                dist_str = f"{dist:.5f}"
                
                se_denom = ((b_tajima - D_prop)**2 * valid_sites)
                if se_denom > 1e-9:
                    var_val = (b_tajima**2 * D_prop * (1.0 - D_prop)) / se_denom
                    if var_val >=0:
                        se_str = f"{math.sqrt(var_val):.5f}"
        except (ValueError, ZeroDivisionError, OverflowError):
            pass # Keep as N/A
    return dist_str, se_str


def calculate_tamura_nei(p1_transitions_pur, p2_transitions_pyr, q_transversions, base_counts_diff, base_counts_ident, valid_sites):
    """Calculates Tamura-Nei distance and its SE."""
    dist_str, se_str = "N/A", "N/A"
    if valid_sites == 0: return dist_str, se_str

    # Base frequencies (e.g., freq_A_trn = (0.5 * nA_tn_diff + nAA_tn_ident) / valid_sites)
    # base_counts_diff: nA_tn_diff, nT_tn_diff, nC_tn_diff, nG_tn_diff
    # base_counts_ident: nAA_tn_ident, nTT_tn_ident, nCC_tn_ident, nGG_tn_ident
    freq_A = (0.5 * base_counts_diff['A'] + base_counts_ident['A']) / valid_sites
    freq_T = (0.5 * base_counts_diff['T'] + base_counts_ident['T']) / valid_sites
    freq_C = (0.5 * base_counts_diff['C'] + base_counts_ident['C']) / valid_sites
    freq_G = (0.5 * base_counts_diff['G'] + base_counts_ident['G']) / valid_sites

    ag_trn = freq_A + freq_G  # Sum of purine frequencies
    tc_trn = freq_T + freq_C  # Sum of pyrimidine frequencies

    P1_prop = p1_transitions_pur / valid_sites
    P2_prop = p2_transitions_pyr / valid_sites
    Q_prop = q_transversions / valid_sites

    k1_trn = (2 * freq_A * freq_G / ag_trn) if ag_trn > 1e-9 else 0
    k2_trn = (2 * freq_T * freq_C / tc_trn) if tc_trn > 1e-9 else 0
    # k3 for the distance formula is 2 * ag_trn * tc_trn
    k3_dist_coeff = 2 * ag_trn * tc_trn

    w1_trn_log_arg = 1.0 - P1_prop / k1_trn - Q_prop / (2 * ag_trn) if k1_trn > 1e-9 and ag_trn > 1e-9 else -1.0
    w2_trn_log_arg = 1.0 - P2_prop / k2_trn - Q_prop / (2 * tc_trn) if k2_trn > 1e-9 and tc_trn > 1e-9 else -1.0
    w3_trn_log_arg = 1.0 - Q_prop / (2 * ag_trn * tc_trn) if ag_trn > 1e-9 and tc_trn > 1e-9 else -1.0
    
    if w1_trn_log_arg > 1e-9 and w2_trn_log_arg > 1e-9 and w3_trn_log_arg > 1e-9:
        try:
            term1 = -k1_trn * math.log(w1_trn_log_arg) if k1_trn > 1e-9 else 0 # if k1_trn is 0, P1_prop must be 0 for w1_arg > 0
            term2 = -k2_trn * math.log(w2_trn_log_arg) if k2_trn > 1e-9 else 0
            term3 = -k3_dist_coeff * math.log(w3_trn_log_arg) if k3_dist_coeff > 1e-9 else 0
            
            dist = term1 + term2 + term3
            dist_str = f"{dist:.5f}"
            
            # SE for Tamura-Nei is very complex, setting to N/A as per fallback plan.
            se_str = "N/A" 
        except (ValueError, ZeroDivisionError, OverflowError):
            pass # Keep as N/A
    return dist_str, se_str


def main():
    parser = argparse.ArgumentParser(
        description="Calculate various DNA distances between all unique pairs of sequences in a FASTA file."
    )
    parser.add_argument(
        "input_fasta_file",
        help="Path to the input FASTA file (aligned sequences)."
    )
    args = parser.parse_args()

    header = (
        "Gene1\tGene2\tTs/Tv_Ratio\tJC_Distance\tSE_JC\t"
        "K2P_Distance\tSE_K2P\tTamuraNei_Distance\tSE_TamuraNei\t"
        "TajimaNei_Distance\tSE_TajimaNei"
    )
    print(header)

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

    if len(records) < 2:
        print("Error: At least two sequences are required.", file=sys.stderr)
        sys.exit(1)

    alignment_length = len(records[0].seq)
    for i, record in enumerate(records):
        if len(record.seq) != alignment_length:
            print(f"Error: Sequences must be aligned and of equal length. "
                  f"'{records[0].id}' (len {alignment_length}) vs '{record.id}' (len {len(record.seq)}).",
                  file=sys.stderr)
            sys.exit(1)
    
    if alignment_length == 0:
        print("Error: Sequences are empty.", file=sys.stderr)
        sys.exit(1)

    purines = {'A', 'G'}
    pyrimidines = {'C', 'T'}
    valid_bases_set = purines.union(pyrimidines)

    for seq1_rec, seq2_rec in itertools.combinations(records, 2):
        s1 = str(seq1_rec.seq).upper()
        s2 = str(seq2_rec.seq).upper()

        p1_transitions_pur = 0
        p2_transitions_pyr = 0
        q_transversions = 0
        valid_sites = 0
        
        # For Tajima-Nei pair counts
        pair_counts_tnj = {
            'AA':0, 'AT':0, 'AG':0, 'AC':0, 'TA':0, 'TT':0, 'TG':0, 'TC':0,
            'GA':0, 'GT':0, 'GG':0, 'GC':0, 'CA':0, 'CT':0, 'CG':0, 'CC':0
        }
        
        # For Tamura-Nei base counts (from differing and identical sites)
        base_counts_diff_trn = {'A':0, 'T':0, 'C':0, 'G':0} # nA_tn_diff etc.
        base_counts_ident_trn = {'A':0, 'T':0, 'C':0, 'G':0} # nAA_tn_ident etc.

        for i in range(alignment_length):
            b1, b2 = s1[i], s2[i]
            if b1 not in valid_bases_set or b2 not in valid_bases_set:
                continue
            valid_sites += 1

            pair_str = b1 + b2
            pair_counts_tnj[pair_str] += 1 # For Tajima-Nei

            if b1 == b2:
                base_counts_ident_trn[b1] += 1 # For Tamura-Nei
            else: # b1 != b2
                # For Tamura-Nei differing site base counts
                base_counts_diff_trn[b1] += 1
                base_counts_diff_trn[b2] += 1

                # Transitions/Transversions
                b1_is_pur = b1 in purines
                b2_is_pur = b2 in purines
                if (b1_is_pur and b2_is_pur): # Purine-Purine transition
                    p1_transitions_pur += 1
                elif (not b1_is_pur and not b2_is_pur): # Pyrimidine-Pyrimidine transition
                    p2_transitions_pyr += 1
                else: # Transversion
                    q_transversions += 1
        
        # Normalize Tajima-Nei pair counts (e.g., AT includes TA)
        # This logic might be simpler if done inside calculate_tajima_nei directly from f[pair_str]
        # The current f['AT'] in calculate_tajima_nei will be count(AT)/valid_sites
        # For fXY, it needs to be (count(XY) + count(YX)) / (2*valid_sites) for X!=Y
        # Or, if pair_counts_tnj is passed directly:
        # fAA = pair_counts_tnj['AA']/valid_sites
        # fAT = (pair_counts_tnj['AT'] + pair_counts_tnj['TA'])/(2*valid_sites)
        # This correction needs to be applied when calculating qA, qT etc.
        # For now, the helper `calculate_tajima_nei` assumes `pair_counts` are raw counts of ordered pairs.
        # Let's simplify by making pair_counts symmetric before passing
        symmetric_pair_counts_tnj = {
            'AA': pair_counts_tnj['AA'], 'TT': pair_counts_tnj['TT'], 
            'GG': pair_counts_tnj['GG'], 'CC': pair_counts_tnj['CC'],
            'AT': pair_counts_tnj['AT'] + pair_counts_tnj['TA'],
            'AG': pair_counts_tnj['AG'] + pair_counts_tnj['GA'],
            'AC': pair_counts_tnj['AC'] + pair_counts_tnj['CA'],
            'TG': pair_counts_tnj['TG'] + pair_counts_tnj['GT'],
            'TC': pair_counts_tnj['TC'] + pair_counts_tnj['CT'],
            'GC': pair_counts_tnj['GC'] + pair_counts_tnj['CG']
        }


        jc_dist, jc_se = "N/A", "N/A"
        k2p_dist, k2p_se, ts_tv_ratio = "N/A", "N/A", "N/A"
        tam_nei_dist, tam_nei_se = "N/A", "N/A"
        taj_nei_dist, taj_nei_se = "N/A", "N/A"

        if valid_sites == 0:
            print(f"Warning: No valid comparable sites between {seq1_rec.id} and {seq2_rec.id}.", file=sys.stderr)
        else:
            P_prop = (p1_transitions_pur + p2_transitions_pyr) / valid_sites
            Q_prop = q_transversions / valid_sites
            
            jc_dist, jc_se = calculate_jc69(P_prop, Q_prop, valid_sites)
            k2p_dist, k2p_se, ts_tv_ratio = calculate_k2p(P_prop, Q_prop, valid_sites)
            # Pass symmetric_pair_counts_tnj, which holds combined counts like 'AT' = count(AT)+count(TA)
            # The helper function then divides by valid_sites (for fAA like) or 2*valid_sites (for fAT like)
            taj_nei_dist, taj_nei_se = calculate_tajima_nei(symmetric_pair_counts_tnj, valid_sites, P_prop, Q_prop)
            tam_nei_dist, tam_nei_se = calculate_tamura_nei(
                p1_transitions_pur, p2_transitions_pyr, q_transversions,
                base_counts_diff_trn, base_counts_ident_trn, valid_sites
            )
            
        print(f"{seq1_rec.id}\t{seq2_rec.id}\t{ts_tv_ratio}\t"
              f"{jc_dist}\t{jc_se}\t{k2p_dist}\t{k2p_se}\t"
              f"{tam_nei_dist}\t{tam_nei_se}\t{taj_nei_dist}\t{taj_nei_se}")

if __name__ == "__main__":
    main()
