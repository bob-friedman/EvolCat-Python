import argparse
import math
import sys # For stderr
from collections import defaultdict
import itertools # For permutations in count_substitutions_for_pair

# Define genetic codes
STANDARD_GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': 'Z', 'TAG': 'Z',
    'TGT': 'C', 'TGC': 'C', 'TGA': 'Z', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

MITO_GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': 'Z', 'TAG': 'Z',
    'TGT': 'C', 'TGC': 'C', 'TGA': 'W', 'TGG': 'W', # TGA is W in mito
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'M', 'ATG': 'M', # ATA is M in mito
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'Z', 'AGG': 'Z', # AGA/AGG are Z in mito
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Base to number mapping for codon conversion
BASE_TO_NUM = {'T': 0, 'C': 1, 'A': 2, 'G': 3}

# Codon sync sites
# These are precomputed values from the Perl scripts.
# The structure is: codon_number: [val1, val2, val3, val4]
# It's not directly used in the provided Python structure, 
# but the logic to calculate these values is within calculate_syn_site.
# For now, these are kept as reference or if a direct lookup approach is chosen later.
CODON_SYN_SITES_STANDARD = {
    # Copied from dsdn_counts_dist.pl @codon_syn_sites
    # Example: 0 => [0.0, 0.0, 0.6666666666666666, 0.0] (for TTT -> F)
    # This will be a large dictionary if fully transcribed.
    # Instead of transcribing, the logic from syn_site will be used.
}

CODON_SYN_SITES_MITO = {
    # Copied from dsdn_counts_dist_mito.pl @codon_syn_sites_mito
    # Similar to above, will rely on calculate_syn_site logic.
}

# Bases for iteration
BASES = ['T', 'C', 'A', 'G']

def parse_fasta(filename):
    """Parses a FASTA file and returns a dictionary of sequences."""
    sequences = {}
    current_seq_name = None
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_seq_name = line[1:]
                sequences[current_seq_name] = []
            elif current_seq_name:
                sequences[current_seq_name].append(line.upper())
    
    for name, seq_parts in sequences.items():
        sequences[name] = "".join(seq_parts)
    return sequences


def validate_sequences(sequences_dict):
    """
    Validates sequences for allowed characters (ACGTN-) and length (multiple of 3).
    Removes gaps ('-') and 'N's from sequences.
    Returns a new dictionary with valid, cleaned sequences, or raises ValueError.
    """
    validated_sequences = {}
    allowed_chars = set("ACGTN-")
    
    first_seq_len = None

    for name, seq_str in sequences_dict.items():
        if not all(char in allowed_chars for char in seq_str):
            raise ValueError(f"Sequence {name} contains invalid characters. Allowed: ACGTN-")

        # Remove gaps and Ns
        cleaned_seq = "".join(c for c in seq_str.upper() if c in "ACGT")

        if not cleaned_seq:
            print(f"Warning: Sequence {name} is empty after removing non-ACGT characters. Skipping.", file=sys.stderr)
            continue

        if len(cleaned_seq) % 3 != 0:
            raise ValueError(f"Sequence {name} (after cleaning) has length {len(cleaned_seq)}, which is not a multiple of 3.")

        if first_seq_len is None:
            first_seq_len = len(cleaned_seq)
        elif len(cleaned_seq) != first_seq_len:
            # This mimics the Perl script's behavior of forcing all sequences to the length of the first one.
            # Consider if this is the desired behavior or if pairs with different lengths should be skipped/error.
            # For now, let's raise an error as unaligned sequences of different lengths are problematic for this analysis.
            raise ValueError(f"Sequence {name} (length {len(cleaned_seq)}) has a different length from the first sequence (length {first_seq_len}) after cleaning. Please provide aligned sequences of equal length.")

        validated_sequences[name] = cleaned_seq
        
    if not validated_sequences:
        raise ValueError("No valid sequences found in the input file after validation.")
        
    return validated_sequences


def get_genetic_code(code_name):
    """Returns the specified genetic code."""
    if code_name == 'standard':
        return STANDARD_GENETIC_CODE
    elif code_name == 'mito':
        return MITO_GENETIC_CODE
    else:
        raise ValueError(f"Unknown genetic code: {code_name}")

def codon_to_aa(codon, genetic_code):
    """Converts a codon to an amino acid using the given genetic code."""
    return genetic_code.get(codon, 'X') # X for unknown codons

def get_codon_number(codon):
    """Converts a DNA codon to a numerical representation (0-63)."""
    num = 0
    for i, base in enumerate(codon):
        num += BASE_TO_NUM[base] * (4**(2-i))
    return num

def estimate_r_value(sequences_list):
    """
    Estimates the transition/transversion ratio (R) from the 3rd codon positions
    of the input sequences. Mirrors the Perl script's logic.
    """
    third_pos_bases = []
    for seq in sequences_list:
        if len(seq) % 3 != 0:
            # Skip sequences not multiple of 3, or handle as error
            # Perl script might implicitly truncate or ignore; being explicit is better
            print(f"Warning: Sequence of length {len(seq)} is not a multiple of 3 and will be skipped for R estimation.")
            continue
        for i in range(2, len(seq), 3): # 3rd position is index 2
            third_pos_bases.append(seq[i])

    if not third_pos_bases:
        print("Warning: No 3rd codon position bases found for R estimation.")
        return None

    counts = {'T': 0, 'C': 0, 'A': 0, 'G': 0}
    total_bases = 0
    for base in third_pos_bases:
        if base in counts:
            counts[base] += 1
            total_bases += 1
    
    if total_bases == 0:
        print("Warning: No valid bases found at 3rd codon positions for R estimation.")
        return None

    # Frequencies
    freq_t = counts['T'] / total_bases
    freq_c = counts['C'] / total_bases
    freq_a = counts['A'] / total_bases
    freq_g = counts['G'] / total_bases

    # Parameters from Perl script (P1, Q, P2)
    # These represent observed proportions of changes, but the Perl script calculates them from frequencies
    # This part of the Perl script is a bit confusing as it seems to be based on expected changes
    # under a model, rather than direct observation of differences between sequences.
    # Let's follow the direct calculation of P, Q, D from frequencies as in Kimura's 2-parameter model
    # The Perl script's P1, Q, P2 seems to be calculating expected proportions of transitions/transversions
    # if all sites were variable.
    
    # Re-interpreting Perl's logic for P, Q based on Kimura (1980) for nucleotide substitution
    # P is proportion of sites differing by transition
    # Q is proportion of sites differing by transversion
    # The Perl script's R estimation is specific and seems to assume equilibrium frequencies.
    
    # Let's try to match the Perl script's variable names and logic as closely as possible.
    # p1 = freq_t * freq_c * 2 (T<->C transitions)
    # q  = (freq_t + freq_c) * (freq_a + freq_g) (T/C <-> A/G transversions)
    # p2 = freq_a * freq_g * 2 (A<->G transitions)

    # The Perl script calculates g1, g2, g3, g4 (frequencies)
    g1, g2, g3, g4 = freq_t, freq_c, freq_a, freq_g # T, C, A, G
    
    # PA = transitions, PB = transversions (these are not P, Q directly from Kimura)
    # This part is from Nei's book or similar, for estimating transition/transversion rates
    # from nucleotide frequencies at equilibrium.
    
    # P1 = Transitions between T and C
    # P2 = Transitions between A and G
    # Q = Transversions between (T,C) and (A,G)
    
    # The Perl script calculates:
    # P1 = $g1*$g2/($g1+$g2); (This is not a proportion, it's a term used in rate estimation)
    # Q = ($g1+$g2)*($g3+$g4);
    # P2 = $g3*$g4/($g3+$g4);
    
    # Avoid division by zero if some base groups are absent
    denom_tc = g1 + g2
    denom_ag = g3 + g4

    if denom_tc == 0 and denom_ag == 0: # All bases are same, or only one type, no changes
        return None # Cannot estimate R
    
    # If only pyrimidines (T,C) or only purines (A,G) exist, Q will be 0.
    # If only T and A exist (for example), P1 and P2 will be 0.

    p1_term = (g1 * g2 / denom_tc) if denom_tc > 0 else 0.0
    q_term = denom_tc * denom_ag # This is Q in Nei & Gojobori's context for divergence
    p2_term = (g3 * g4 / denom_ag) if denom_ag > 0 else 0.0
    
    # These are terms used in calculating A, B, C for Kimura's 3ST model or similar.
    # The Perl script then calculates P, Q for divergence (D)
    # P = ($p[0]*$p[1]+$p[2]*$p[3]) / ( ($p[0]+$p[1])*($p[2]+$p[3]) + ($p[0]+$p[1])**2 + ($p[2]+$p[3])**2 );
    # Q = 1 - P;
    # This P and Q are different. It seems to be proportion of sites that are transitions vs transversions.

    # Let's assume the goal is to get proportions of transitional and transversional differences
    # for Kimura's two-parameter model (K2P).
    # For K2P, P is proportion of sites with transitional diffs, Q for transversional.
    # This requires pairwise comparison of sequences, not just frequencies.
    
    # The Perl script's R estimation is quite specific. It uses this:
    # w1 = $g1+$g2; # freq of pyrimidines
    # w2 = $g3+$g4; # freq of purines
    # s = ($p1_term + $p2_term) / ( ($w1**2 * $w2**2) / $q_term ) if $q_term > 0 else some_value
    # v = 0.5 / ( ($w1**2 * $w2**2) / $q_term ) if $q_term > 0 else some_value
    # ratio = s/v
    
    # This seems to be from Kimura (1980) "A simple method for estimating evolutionary rates..."
    # Equation 8: alpha = (1/(2*T)) * ln(X/(1-2P)) and beta = (1/(2*T)) * ln(1/(1-2Q))
    # R = alpha / beta.
    # P and Q there are proportions of sites that differ by transition/transversion.
    # The Perl script's `estimate_R` seems to be calculating something else or using a shortcut.

    # Let's simplify and use the direct formulas for P and Q (proportions of differences)
    # This means we need to compare pairs of 3rd position bases.
    # The Perl script's `estimate_R` does *not* do pairwise comparison for P and Q.
    # It calculates P and Q based on *frequencies* g1,g2,g3,g4.

    # Trying to exactly replicate Perl's estimate_R:
    g = [freq_t, freq_c, freq_a, freq_g] # T, C, A, G as in Perl's @p
    
    # Denominators for P1, P2 terms, check for zero
    d1 = g[0] + g[1] # T+C
    d2 = g[2] + g[3] # A+G

    if d1 == 0 or d2 == 0: # Only purines or only pyrimidines, so no transversions possible in this model
        # This situation would make Q_kimura (transversion rate) zero, R undefined or infinite.
        # The Perl script would likely lead to division by zero here for $Q_kimura
        print("Warning: Only purines or only pyrimidines found at 3rd pos. R estimation may be unreliable/undef.")
        # If d1=0, P_kimura for T-C transitions is undefined. If d2=0, P_kimura for A-G is undefined.
        # The script calculates P (overall transitions) and Q (overall transversions)
        # P_overall = g[0]*g[1]/d1 + g[2]*g[3]/d2; (This is not a probability but a sum of terms)
        # Q_overall = (g[0]+g[1])*(g[2]+g[3]); 
        # The Perl script's P and Q for divergence calculation (D) is:
        # P_divergence = (g[0]*g[1] + g[2]*g[3]) / ((g[0]+g[1])*(g[2]+g[3]) + (g[0]+g[1])**2 + (g[2]+g[3])**2)
        # This is very specific. Let's assume it's from a particular model.
        # If d1 or d2 is zero, it implies Q_overall is zero.
        # If Q_overall is zero, then w1=0 or w2=0.
        # Then D calculation will have division by zero.
        return None # Cannot reliably estimate R

    # P_kimura and Q_kimura are the proportions of sites that differ by transition/transversion.
    # The Perl script uses these:
    # P_trans_TC = 2 * g[0] * g[1] / d1  (Expected transitions T-C if sites were T or C)
    # P_trans_AG = 2 * g[2] * g[3] / d2  (Expected transitions A-G if sites were A or G)
    # Q_transver = 2 * d1 * d2          (Expected transversions)

    # These are used to calculate expected P (transition rate) and Q (transversion rate)
    # from Kimura's model using equilibrium frequencies.
    # P_rate_numerator = (g[0]*g[1]/d1 + g[2]*g[3]/d2)
    # P_rate_denominator = ( (g[0]+g[1])**2 * (g[2]+g[3])**2 ) / ( (g[0]+g[1])*(g[2]+g[3]) ) # Simplified later
    # Q_rate_numerator = 0.5
    # Q_rate_denominator = ( (g[0]+g[1])**2 * (g[2]+g[3])**2 ) / ( (g[0]+g[1])*(g[2]+g[3]) )
    
    # Let P_exp and Q_exp be expected proportions of transitional/transversional differences
    # P_exp = 2 * (g[0]*g[1]/d1 * d1**2 + g[2]*g[3]/d2 * d2**2) # This is not right
    # P_exp = g[0]*g[1] + g[2]*g[3] # Sum of ways to have transition
    # Q_exp = (g[0]+g[1])*(g[2]+g[3]) # Sum of ways to have transversion

    # The formulas in the Perl script for s, v, w1, w2, D are from Nei and Gojobori (1986)
    # equations A5, A6 for estimating transition (s) and transversion (v) rates per site.
    # Here, P, Q are actual counts of transitional/transversional differences, not frequencies.
    # The Perl script's estimate_R seems to be a mix or a specific variant.

    # Sticking to the Perl code's `estimate_R` direct translation:
    # $p[0] = T, $p[1] = C, $p[2] = A, $p[3] = G frequencies
    if d1 == 0: term1 = 0.0
    else: term1 = g[0]*g[1]/d1
    
    if d2 == 0: term2 = 0.0
    else: term2 = g[2]*g[3]/d2
        
    P_val = term1 + term2 # Corresponds to $P in Perl's estimate_R
    Q_val = d1 * d2       # Corresponds to $Q in Perl's estimate_R

    if Q_val == 0: # No transversions possible based on frequencies (e.g. only T and C present)
        return None # R would be undefined or infinite

    try:
        # D = -0.5 * log(1 - P_val/(d1*d2) - Q_val/((d1+d2)**2 - d1**2 - d2**2) ); # This is not what Perl does.
        # Perl D = -0.5*log(1 - $P_val/$Q_val - $Q_val / ( ( ($p[0]+$p[1])**2*($p[2]+$p[3])**2 )/$Q_val ) );
        # This D is not Kimura's D. It's part of a specific R estimation.
        
        # w1 = d1 (freq pyrimidines), w2 = d2 (freq purines)
        # s_numerator = P_val
        # v_numerator = Q_val / 2.0 # (Perl: $Q/2)
        
        # Denominator for both s and v rates in that model:
        # denom_sv = 1.0 - (g[0]**2/d1 + g[1]**2/d1) - (g[2]**2/d2 + g[3]**2/d2) # This is not in Perl.
        # The Perl script calculation for D (divergence) is:
        # $D = -0.5*log(1 - $P/$Q - $Q / ( (($p[0]+$p[1])**2*($p[2]+$p[3])**2)/$Q ) );
        # This D is not used for s/v directly.
        # It calculates s and v as:
        # s = (1/$D) * 0.5 * log(1/(1 - 2*$P_val_prime - $Q_val_prime));
        # v = (1/$D) * 0.25 * log(1/(1 - 2*$Q_val_prime));
        # Where P_val_prime and Q_val_prime are proportions of sites with transi/transver diffs.
        # This means the Perl script's estimate_R is NOT using the P and Q (freq-based terms)
        # directly to calculate s and v as R = (P/d1 + P/d2) / (Q/2).
        
        # It seems the P and Q in estimate_R are intermediate terms for something else.
        # The final $ratio = $s/$v uses $s and $v calculated based on P,Q from count_diffs_for_R.
        # And count_diffs_for_R *is* doing pairwise comparisons.
        
        # So, estimate_R in Perl is only for setting up genetic code and calling count_diffs_for_R.
        # The actual R calculation (s/v) uses P and Q from pairwise comparisons.
        
        # The task is to mirror Perl's R estimation.
        # If 'R' is passed, the Perl script calls estimate_R, which then calls count_diffs_for_R.
        # count_diffs_for_R:
        #   - takes all sequences
        #   - iterates all pairs of sequences
        #   - for each pair, counts actual P (transitions) and Q (transversions) at 3rd sites
        #   - averages these P and Q over all pairs
        #   - then calculates s, v, and ratio = s/v using these averaged P, Q.

        # This means `estimate_r_value` needs to perform pairwise comparisons on 3rd sites.
        
        num_seqs = len(sequences_list)
        if num_seqs < 2:
            print("Warning: R estimation requires at least 2 sequences.")
            return None

        total_transitions = 0
        total_transversions = 0
        total_comparisons = 0 # Total number of 3rd sites compared

        for i in range(num_seqs):
            for j in range(i + 1, num_seqs):
                seq1 = sequences_list[i]
                seq2 = sequences_list[j]

                # Ensure sequences are valid for comparison
                if len(seq1) % 3 != 0 or len(seq2) % 3 != 0 or len(seq1) != len(seq2):
                    # Perl script seems to use first sequence length as reference.
                    # For simplicity here, require equal, valid lengths for pairs in R estimation.
                    continue 

                codons1 = [seq1[k:k+3] for k in range(0, len(seq1), 3)]
                codons2 = [seq2[k:k+3] for k in range(0, len(seq2), 3)]

                for c_idx in range(len(codons1)):
                    base1 = codons1[c_idx][2] # 3rd base
                    base2 = codons2[c_idx][2] # 3rd base

                    if base1 not in BASE_TO_NUM or base2 not in BASE_TO_NUM:
                        continue # Skip if invalid base

                    if base1 == base2:
                        continue
                    
                    total_comparisons += 1
                    
                    # Check for transition/transversion
                    b1_is_purine = base1 in ['A', 'G']
                    b2_is_purine = base2 in ['A', 'G']

                    if (b1_is_purine and b2_is_purine) or (not b1_is_purine and not b2_is_purine):
                        total_transitions += 1
                    else:
                        total_transversions += 1
        
        if total_comparisons == 0:
            print("Warning: No differing 3rd codon sites found for R estimation.")
            return None # Avoid division by zero later

        P_prime = total_transitions / total_comparisons # Proportion of transitional differences
        Q_prime = total_transversions / total_comparisons # Proportion of transversional differences

        # Now use Kimura's formulas for s and v (alpha and beta in his paper)
        # s = alpha rate, v = beta rate
        # alpha = -0.5 * log(1 - 2*P_prime - Q_prime)  (This is for divergence, not rate directly)
        # beta  = -0.5 * log(1 - 2*Q_prime)
        
        # The Perl script uses specific formulas from Nei & Gojobori (1986) or similar:
        # w1 = 1 - 2*P_prime - Q_prime
        # w2 = 1 - 2*Q_prime
        # If w1 <= 0 or w2 <= 0, log will fail.
        
        val_for_s_log = 1.0 - 2.0 * P_prime - Q_prime
        val_for_v_log = 1.0 - 2.0 * Q_prime

        if val_for_s_log <= 0 or val_for_v_log <= 0:
            print("Warning: Cannot compute log for s or v in R estimation (P', Q' too high).")
            return None

        s_rate = 0.5 * math.log(1.0 / val_for_s_log) # Renamed from 's' to 's_rate' to avoid conflict
        v_rate = 0.25 * math.log(1.0 / val_for_v_log) # Renamed from 'v' to 'v_rate'
        
        # Perl: $s = (1/$D) * ... ; $v = (1/$D) * ...
        # But $D (divergence time) cancels out in $s/$v.
        # So, ratio = (0.5 * log(1/(1-2P-Q))) / (0.25 * log(1/(1-2Q)))
        # Which is 2 * log(1/(1-2P-Q)) / log(1/(1-2Q))

        if v_rate == 0:
            # If Q_prime is 0, v_rate is 0. R is infinite or undefined.
            # Perl script would output 'undef'.
            print("Warning: Estimated v_rate is 0 in R estimation.")
            return None 
            
        estimated_r = s_rate / v_rate
        
        if estimated_r < 0: # Should not happen if logs are valid
            print("Warning: Estimated R is negative.")
            return None
            
        return estimated_r

    except Exception as e:
        print(f"Error during R estimation calculation: {e}")
        return None


def calculate_syn_site(codon, genetic_code, r_value):
    """
    Calculates potential synonymous sites for a given codon.
    This replicates the logic of the syn_site subroutine in the Perl scripts.
    """
    if len(codon) != 3 or not all(b in BASE_TO_NUM for b in codon):
        return 0.0 # Invalid codon

    original_aa = codon_to_aa(codon, genetic_code)
    if original_aa == 'Z': # Stop codons have no synonymous sites
        return 0.0

    syn_sites = 0.0
    
    # Iterate over each position in the codon
    for i in range(3):
        original_base = codon[i]
        # Iterate over possible mutations (T, C, A, G)
        for new_base_char in BASES:
            if new_base_char == original_base:
                continue

            # Create the new codon
            mutated_codon_list = list(codon)
            mutated_codon_list[i] = new_base_char
            mutated_codon = "".join(mutated_codon_list)

            mutated_aa = codon_to_aa(mutated_codon, genetic_code)

            if original_aa == mutated_aa:
                # Check if transition or transversion
                # Transitions: A <-> G, C <-> T
                is_transition = False
                if (original_base in ['A', 'G'] and new_base_char in ['A', 'G']) or \
                   (original_base in ['C', 'T'] and new_base_char in ['C', 'T']):
                    is_transition = True
                
                if is_transition:
                    syn_sites += r_value
                else:
                    syn_sites += 1.0 
                    
    # The Perl script divides by sum of (R * num_transitions + num_transversions)
    # This normalization factor depends on R and the codon's possible changes.
    # Let's calculate that normalization factor, which is effectively the denominator in Perl's $codon_syn_sites[$codon_num][1]
    
    norm_factor = 0.0
    for i in range(3):
        original_base = codon[i]
        for new_base_char in BASES:
            if new_base_char == original_base:
                continue
            is_transition = False
            if (original_base in ['A', 'G'] and new_base_char in ['A', 'G']) or \
               (original_base in ['C', 'T'] and new_base_char in ['C', 'T']):
                is_transition = True
            
            if is_transition:
                norm_factor += r_value
            else:
                norm_factor += 1.0
                
    if norm_factor == 0: # Should not happen for valid codons
        return 0.0

    return syn_sites / norm_factor if norm_factor > 0 else 0.0


def get_diff_positions(codon1, codon2):
    """Returns a list of positions where two codons differ."""
    diffs = []
    for i in range(3):
        if codon1[i] != codon2[i]:
            diffs.append(i)
    return diffs

def count_substitutions_for_pair(seq1_name, seq1, seq2_name, seq2, genetic_code, r_value):
    """
    Counts synonymous and non-synonymous substitutions for a pair of sequences.
    Mirrors the logic of countsubstitutions in the Perl scripts.
    """
    syn_subs = 0.0  # Use float for syn/nonsyn counts as they can be averaged
    nonsyn_subs = 0.0
    
    num_codons = len(seq1) // 3 # Assumes seq1 and seq2 have same length, multiple of 3

    for i in range(num_codons):
        codon1_str = seq1[i*3 : (i+1)*3]
        codon2_str = seq2[i*3 : (i+1)*3]

        if not all(b in BASE_TO_NUM for b in codon1_str) or \
           not all(b in BASE_TO_NUM for b in codon2_str):
            # Skip if codon contains invalid characters (e.g., gaps after initial validation)
            # The main validation should catch this earlier for whole sequences.
            continue

        if codon1_str == codon2_str:
            continue

        aa1 = codon_to_aa(codon1_str, genetic_code)
        aa2 = codon_to_aa(codon2_str, genetic_code)

        if aa1 == 'Z' or aa2 == 'Z': # If either original codon is a stop codon, skip (as per Perl)
            continue

        diff_positions = get_diff_positions(codon1_str, codon2_str)
        num_diffs = len(diff_positions)

        if num_diffs == 1:
            pos = diff_positions[0]
            # Determine if transition or transversion for weighting
            is_transition = False
            original_base = codon1_str[pos]
            mutated_base = codon2_str[pos]
            if (original_base in ['A', 'G'] and mutated_base in ['A', 'G']) or \
               (original_base in ['C', 'T'] and mutated_base in ['C', 'T']):
                is_transition = True
            
            weight = r_value if is_transition else 1.0
            
            if aa1 == aa2: # Synonymous
                syn_subs += weight
            else: # Non-synonymous
                nonsyn_subs += weight
        
        elif num_diffs == 2:
            # Average over two possible pathways
            path_syn_sum = 0.0
            path_nonsyn_sum = 0.0
            num_paths = 0

            for k in range(2): # Two intermediate codons
                intermediate_codon_list = list(codon1_str)
                # Path 1: Change first differing base
                intermediate_codon_list[diff_positions[k]] = codon2_str[diff_positions[k]]
                intermediate_codon1_str = "".join(intermediate_codon_list)
                
                # Path 2: Change second differing base (relative to codon1)
                # This is effectively creating the other intermediate:
                # codon1 -> intermediate_codon2_str -> codon2
                intermediate_codon_list_alt = list(codon1_str)
                intermediate_codon_list_alt[diff_positions[1-k]] = codon2_str[diff_positions[1-k]]
                intermediate_codon2_str = "".join(intermediate_codon_list_alt)


                # Path: codon1 -> intermediate_codon1_str -> codon2
                aa_intermediate1 = codon_to_aa(intermediate_codon1_str, genetic_code)
                
                if aa_intermediate1 == 'Z': # Path via stop codon is not counted
                    continue 
                
                num_paths +=1 # Count this valid path

                # Step 1: codon1 -> intermediate_codon1_str
                original_base_s1 = codon1_str[diff_positions[k]]
                mutated_base_s1 = intermediate_codon1_str[diff_positions[k]]
                is_transition_s1 = (original_base_s1 in ['A','G'] and mutated_base_s1 in ['A','G']) or \
                                   (original_base_s1 in ['C','T'] and mutated_base_s1 in ['C','T'])
                weight1 = r_value if is_transition_s1 else 1.0

                if aa1 == aa_intermediate1:
                    path_syn_sum += weight1
                else:
                    path_nonsyn_sum += weight1

                # Step 2: intermediate_codon1_str -> codon2_str
                # The differing base for this step is at diff_positions[1-k]
                original_base_s2 = intermediate_codon1_str[diff_positions[1-k]]
                mutated_base_s2 = codon2_str[diff_positions[1-k]]
                is_transition_s2 = (original_base_s2 in ['A','G'] and mutated_base_s2 in ['A','G']) or \
                                   (original_base_s2 in ['C','T'] and mutated_base_s2 in ['C','T'])
                weight2 = r_value if is_transition_s2 else 1.0

                if aa_intermediate1 == aa2:
                    path_syn_sum += weight2
                else:
                    path_nonsyn_sum += weight2
            
            if num_paths > 0:
                syn_subs += path_syn_sum / num_paths
                nonsyn_subs += path_nonsyn_sum / num_paths

        elif num_diffs == 3:
            # Average over six possible pathways (3! = 6)
            # Path: c1 -> i1 -> i2 -> c2
            path_syn_sum = 0.0
            path_nonsyn_sum = 0.0
            num_valid_paths = 0
            
            # Permutations of the differing positions (0, 1, 2 for the diff_positions list)
            import itertools
            # diff_positions already contains the indices that are different, e.g., [0, 1, 2] if all differ
            
            for p in itertools.permutations(diff_positions):
                # p is a permutation of the indices where codons differ.
                # e.g., if codon1=TTT, codon2=CCC, diff_positions=[0,1,2]. p could be (0,1,2)
                # Path step 1: codon1 -> intermediate1 (change at p[0])
                intermediate1_list = list(codon1_str)
                intermediate1_list[p[0]] = codon2_str[p[0]]
                intermediate1_str = "".join(intermediate1_list)
                aa_intermediate1 = codon_to_aa(intermediate1_str, genetic_code)
                if aa_intermediate1 == 'Z': continue

                original_base_s1 = codon1_str[p[0]]
                mutated_base_s1 = intermediate1_str[p[0]]
                is_transition_s1 = (original_base_s1 in ['A','G'] and mutated_base_s1 in ['A','G']) or \
                                   (original_base_s1 in ['C','T'] and mutated_base_s1 in ['C','T'])
                weight1 = r_value if is_transition_s1 else 1.0
                
                current_path_syn = 0.0
                current_path_nonsyn = 0.0

                if aa1 == aa_intermediate1: current_path_syn += weight1
                else: current_path_nonsyn += weight1

                # Path step 2: intermediate1 -> intermediate2 (change at p[1])
                intermediate2_list = list(intermediate1_str)
                intermediate2_list[p[1]] = codon2_str[p[1]]
                intermediate2_str = "".join(intermediate2_list)
                aa_intermediate2 = codon_to_aa(intermediate2_str, genetic_code)
                if aa_intermediate2 == 'Z': continue
                
                original_base_s2 = intermediate1_str[p[1]]
                mutated_base_s2 = intermediate2_str[p[1]]
                is_transition_s2 = (original_base_s2 in ['A','G'] and mutated_base_s2 in ['A','G']) or \
                                   (original_base_s2 in ['C','T'] and mutated_base_s2 in ['C','T'])
                weight2 = r_value if is_transition_s2 else 1.0

                if aa_intermediate1 == aa_intermediate2: current_path_syn += weight2
                else: current_path_nonsyn += weight2
                
                # Path step 3: intermediate2 -> codon2 (change at p[2])
                original_base_s3 = intermediate2_str[p[2]] # Base from intermediate2
                mutated_base_s3 = codon2_str[p[2]]    # Base from final codon2
                is_transition_s3 = (original_base_s3 in ['A','G'] and mutated_base_s3 in ['A','G']) or \
                                   (original_base_s3 in ['C','T'] and mutated_base_s3 in ['C','T'])
                weight3 = r_value if is_transition_s3 else 1.0

                if aa_intermediate2 == aa2: current_path_syn += weight3
                else: current_path_nonsyn += weight3
                
                num_valid_paths += 1
                path_syn_sum += current_path_syn
                path_nonsyn_sum += current_path_nonsyn

            if num_valid_paths > 0:
                syn_subs += path_syn_sum / num_valid_paths
                nonsyn_subs += path_nonsyn_sum / num_valid_paths
                
    # Calculate potential synonymous and non-synonymous sites (SA_Nei, SB_Nei)
    # SA_Nei and SB_Nei are sums of f_i values for each codon in seq1 and seq2 respectively
    # where f_i is the fraction of synonymous sites for codon i.
    SA_Nei = sum(calculate_syn_site(seq1[i*3:(i+1)*3], genetic_code, r_value) 
                 for i in range(num_codons) if codon_to_aa(seq1[i*3:(i+1)*3], genetic_code) != 'Z')
    SB_Nei = sum(calculate_syn_site(seq2[i*3:(i+1)*3], genetic_code, r_value)
                 for i in range(num_codons) if codon_to_aa(seq2[i*3:(i+1)*3], genetic_code) != 'Z')
    # Number of codons considered for potential sites (excluding those that are stop codons in original)
    # This count might differ slightly from num_codons if original sequences had stop codons.
    # However, the Perl script uses $num_codons (total codons in alignment) for potential_nonsyn.
    # Let's stick to num_codons as per Perl for `potential_nonsyn` calculation.
    # The SA_Nei/SB_Nei sums already exclude contributions from initial stop codons.
    
    if num_codons > 0:
        # The Perl output divides SA_Nei and SB_Nei by 3.
        # This is because their syn_site function returns sum for 3 positions, 
        # while calculate_syn_site here returns a fraction (0 to 1) for the whole codon.
        # So SA_Nei/3 in Perl is equivalent to sum(f_i)/num_codons if f_i was per site.
        # Here, SA_Nei = sum(f_i_codon).
        # The formula (SA_Nei/3 + SB_Nei/3)/2 in Perl becomes (sum(f_i_c1) + sum(f_i_c2))/2 in Python,
        # which is what Nei & Gojobori's S = (S1+S2)/2 is, where S1 = sum(f_i) for seq1.
        # The division by 3 in Perl output: $SA_Nei = sprintf("%1.2f", $SA_Nei/3);
        # This means the values SA_Nei, SB_Nei printed by Perl are 1/3 of what's used in potential_syn.
        # The Python script should output SA_Nei and SB_Nei directly as sum of f_i for each sequence.
        # But for `potential_syn` calculation, we need to match the Perl's $S_sites = ($SA_Nei + $SB_Nei) / 2;
        # where $SA_Nei and $SB_Nei were already divided by 3.
        # So, potential_syn = ( (sum_fi_seq1/3) + (sum_fi_seq2/3) ) / 2
        # This is equivalent to (sum_fi_seq1 + sum_fi_seq2) / 6
        # Let's re-check Nei & Gojobori (1986).
        # S = sum over codons (L_i) / (number of codons). L_i is number of syn sites at codon i (0 to 3).
        # The syn_site in Perl calculates L_i effectively. So SA_Nei in Perl is sum(L_i).
        # Then potential_syn in Perl is ( (sum(L_i_s1)/3) + (sum(L_i_s2)/3) ) / 2.
        # This is ( sum(L_i_s1) + sum(L_i_s2) ) / 6.
        # My calculate_syn_site returns f_i (fraction, L_i / (total paths weighted)).
        # The Perl script's $codon_syn_sites[$codon_num][0] is sum of (paths * R_or_1)
        # $codon_syn_sites[$codon_num][1] is sum of (R_or_1 for all changes)
        # syn_site returns [0]/[1] which is the f_i.
        # So SA_Nei in Perl is sum(f_i).
        # Then potential_syn = (sum(f_i_s1) + sum(f_i_s2)) / 2. This matches my current code.
        # The division by 3 for printing SA_Nei, SB_Nei in Perl is separate.
        # The formula for potential_nonsyn = (num_codons * 3 - potential_syn_sites) is correct
        # if potential_syn_sites is the sum of L_i values (0 to 3 scale).
        # If potential_syn is sum of f_i (0 to 1 scale), then potential_nonsyn needs adjustment.
        # Let's assume calculate_syn_site gives f_i (fraction of changes that are synonymous from this codon).
        # Then SA_Nei = sum(f_i for seq1), SB_Nei = sum(f_i for seq2).
        # S_sites (total number of synonymous sites) = (SA_Nei + SB_Nei) / 2. This S_sites is on the scale of codons (0 to num_codons).
        # N_sites (total number of non-synonymous sites) = (num_codons * 3) - S_sites_scaled_to_nucleotides
        
        # Re-reading Nei & Gojobori (1986) p.419:
        # L_i = number of synonymous sites at the i-th codon (e.g. 1/3, 2/3, 1, 4/3, ... up to 3)
        # s_ij = L_i for codon j in sequence i.
        # S = (1/m) * sum over j=1 to m ( (s_1j + s_2j)/2 ) -- this is average number of syn sites per codon
        # Total number of synonymous sites in the two sequences = m * S = sum_j ( (s_1j + s_2j)/2 )
        # This is `potential_syn` if s_ij are L_i values.
        # The Perl `syn_site` returns $s/$N which is f_i (fraction of synonymous changes from a codon, weighted by R).
        # This f_i is not L_i directly. L_i is the count of sites.
        # The $SA_Nei, $SB_Nei in Perl are sum(f_i).
        # The output $SA_Nei/3 and $SB_Nei/3.
        # $S_sites = ($SA_Nei + $SB_Nei)/2. This is sum of f_i values, averaged.
        # $N_sites = $num_codons*3 - $S_sites. This is the formula used in Perl.
        # This implies $S_sites is treated as the total number of synonymous sites.
        # This means f_i values are directly used as if they are L_i values.
        # This has been a common interpretation/simplification. Let's stick to it.
        
        potential_syn_sites = (SA_Nei + SB_Nei) / 2.0
        # Ensure potential_syn_sites is not more than total sites (num_codons * 3)
        # This should not happen if calculate_syn_site returns a fraction <=1, as SA_Nei <= num_codons.
        # Then potential_syn_sites <= num_codons.
        
        # If SA_Nei is sum(f_i) where f_i is fraction of changes from codon i that are synonymous,
        # then potential_syn_sites is average number of "synonymous codons" if we consider fractions.
        # The number of non-synonymous sites is total_sites - synonymous_sites.
        # Total sites = num_codons * 3 (at nucleotide level).
        # The S_sites in Perl is sum of f_i for codons. Max value is num_codons.
        # N_sites = num_codons*3 - S_sites. This is what Perl does.
        
        potential_nonsyn_sites = (num_codons * 3.0) - potential_syn_sites
        
        # Defensive check:
        if potential_syn_sites < 0: potential_syn_sites = 0
        if potential_nonsyn_sites < 0: potential_nonsyn_sites = 0 # Should not happen if potential_syn_sites <= num_codons*3
        
    else:
        potential_syn_sites = 0
        potential_nonsyn_sites = 0

    return syn_subs, nonsyn_subs, potential_syn_sites, potential_nonsyn_sites


def main():
    parser = argparse.ArgumentParser(description="Calculate dS/dN ratios for sequences in a FASTA file.")
    parser.add_argument("input_file", help="Path to the input FASTA file.")
    parser.add_argument("ratio", nargs='?', help="Transition/transversion ratio (R). Use 'R' to estimate from data or provide a float.", default="0.5")
    parser.add_argument("--genetic_code", choices=['standard', 'mito'], default='standard', help="Genetic code to use.")
    
    args = parser.parse_args()

    raw_sequences = {}
    try:
        raw_sequences = parse_fasta(args.input_file)
    except FileNotFoundError:
        print(f"Error: Input file '{args.input_file}' not found.", file=sys.stderr)
        return
    except Exception as e:
        print(f"Error parsing FASTA file: {e}", file=sys.stderr)
        return

    sequences = {}
    try:
        sequences = validate_sequences(raw_sequences)
    except ValueError as e:
        print(f"Error validating sequences: {e}", file=sys.stderr)
        return
        
    if len(sequences) < 2:
        print("Error: At least two valid sequences are required for comparison.", file=sys.stderr)
        return

    selected_genetic_code = get_genetic_code(args.genetic_code)
    
    r_value_to_use = 0.0
    r_value_to_use = None # Initialize to handle 'undef' case
    r_estimation_failed_globally = False

    if isinstance(args.ratio, str) and args.ratio.upper() == 'R':
        # Check if there are enough sequences for R estimation
        if len(sequences) < 2: # Should have been caught by earlier check, but good to be safe
            print("Error: At least two sequences are required to estimate R.", file=sys.stderr)
            return
        
        estimated_r = estimate_r_value(list(sequences.values()))
        if estimated_r is None or estimated_r <= 0:
            print("Warning: R estimation resulted in an undefined or non-positive value. 'undef' will be used for the ratio in output.", file=sys.stderr)
            r_estimation_failed_globally = True 
            # No specific r_value_to_use is set yet, will be handled per pair.
        else:
            r_value_to_use = estimated_r
    else:
        try:
            r_value_to_use = float(args.ratio)
            if r_value_to_use <= 0:
                print("Error: Ratio must be a positive float or 'R'.", file=sys.stderr)
                return
        except ValueError:
            print(f"Error: Invalid ratio value '{args.ratio}'. Must be a positive float or 'R'.", file=sys.stderr)
            return
    
    seq_names = list(sequences.keys())


    for i in range(len(seq_names)):
        for j in range(i + 1, len(seq_names)):
            seq1_name = seq_names[i]
            seq2_name = seq_names[j]
            seq1 = sequences[seq1_name]
            seq2 = sequences[seq2_name]

            # Basic length check - needs refinement based on alignment strategy
            if len(seq1) != len(seq2) or len(seq1) % 3 != 0:
                # This check needs to be more robust, considering gaps and alignment.
                # The Perl script seems to assume aligned sequences of equal length.
                # For now, skip pairs with different lengths or non-codon lengths.
                print(f"Skipping pair {seq1_name}-{seq2_name} due to length mismatch or non-multiple of 3 length.")
                continue
            
            current_r_for_pair = r_value_to_use
            formatted_ratio_str = "undef"

            if r_estimation_failed_globally:
                # If global R estimation failed, formatted_ratio_str remains "undef"
                # and current_r_for_pair is None. Calculations depending on R will need to handle None.
                # However, the problem asks to print "undef" for the ratio and continue.
                # The count_substitutions_for_pair function expects a float R.
                # The Perl script, if R estimation fails, prints undef for ratio and seems to stop for that pair or all.
                # "if ($ratio eq "undef") {print "$seq1\t$seq2\tundef\n"; next;}"
                # This implies if R is "undef", we don't proceed with syn/nonsyn counts for the pair.
                print(f"{seq1_name}\t{seq2_name}\tundef")
                continue
            
            # If R was successfully provided or estimated globally
            if current_r_for_pair is not None:
                 formatted_ratio_str = f"{current_r_for_pair:.2f}"
            else:
                # This case should ideally not be reached if r_estimation_failed_globally is true
                # Or if a fixed ratio was invalid (caught earlier).
                # If it's reached, it means a specific R for a pair failed, which is not current logic.
                print(f"Error: R value is unexpectedly None for pair {seq1_name}-{seq2_name} when not globally failed.", file=sys.stderr)
                print(f"{seq1_name}\t{seq2_name}\tundef")
                continue


            # The core calculation logic will go here
            # Ensure current_r_for_pair is a float before passing
            if not isinstance(current_r_for_pair, float):
                 print(f"Error: R value '{current_r_for_pair}' is not a float for pair {seq1_name}-{seq2_name}.", file=sys.stderr)
                 print(f"{seq1_name}\t{seq2_name}\tundef")
                 continue

            syn_count, nonsyn_count, pot_syn_sites, pot_nonsyn_sites = count_substitutions_for_pair(
                seq1_name, seq1, seq2_name, seq2, 
                selected_genetic_code, current_r_for_pair
            )

            # Handle division by zero for potential_syn or potential_nonsyn
            # The Perl script exits with "denom 0". Python script should print specific message for the pair.
            # It seems the Perl script condition is `if ($S_sites == 0 || $N_sites == 0)`
            if pot_syn_sites == 0 or pot_nonsyn_sites == 0:
                # Output matches Perl's "denom 0" like message, but includes values.
                # Perl has: print "$seq1\t$seq2\t$ratio\tdenom 0\n"; # And then it typically exits or skips the pair.
                # Here, we print for the pair and continue.
                print(f"{seq1_name}\t{seq2_name}\t{formatted_ratio_str}\tdenom_0_syn={pot_syn_sites:.2f}_nonsyn={pot_nonsyn_sites:.2f}")
                continue

            # Output format: gene1_name gene2_name formatted_ratio syn_codons nonsyn_codons potential_syn potential_nonsyn
            # syn_codons and nonsyn_codons are $S_changes, $N_changes in Perl (actual counts/sums of weights)
            # potential_syn and potential_nonsyn are $S_sites, $N_sites in Perl
            print(f"{seq1_name}\t{seq2_name}\t{formatted_ratio_str}\t{syn_count:.2f}\t{nonsyn_count:.2f}\t{pot_syn_sites:.2f}\t{pot_nonsyn_sites:.2f}")


if __name__ == "__main__":
    main()
