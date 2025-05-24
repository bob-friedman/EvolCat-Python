import argparse
import math
import sys
from collections import defaultdict
import itertools

# --- Constants and Genetic Code (similar to calculate_ds_dn.py) ---
STANDARD_GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': 'Z', 'TAG': 'Z', # Z for Stop
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

BASE_TO_NUM = {'T': 0, 'C': 1, 'A': 2, 'G': 3}
BASES = ['T', 'C', 'A', 'G']

# --- Helper Functions (some can be reused or adapted from calculate_ds_dn.py) ---

def parse_fasta(filename):
    """Parses a FASTA file and returns a dictionary of sequences."""
    sequences = {}
    current_seq_name = None
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    current_seq_name = line[1:]
                    sequences[current_seq_name] = []
                elif current_seq_name:
                    sequences[current_seq_name].append(line.upper())
    except FileNotFoundError:
        raise
    except Exception as e:
        raise Exception(f"Error reading FASTA file {filename}: {e}")
    
    for name, seq_parts in sequences.items():
        sequences[name] = "".join(seq_parts)
    return sequences

def validate_sequences(sequences_dict):
    """
    Validates sequences for allowed characters (ACGTN-) and length (multiple of 3).
    Removes gaps ('-') and 'N's. Ensures all sequences have the same length.
    """
    validated_sequences = {}
    allowed_chars = set("ACGTN-")
    first_seq_len = None

    for name, seq_str in sequences_dict.items():
        if not all(char in allowed_chars for char in seq_str):
            raise ValueError(f"Sequence {name} contains invalid characters. Allowed: ACGTN-")
        
        cleaned_seq = "".join(c for c in seq_str.upper() if c in "ACGT")

        if not cleaned_seq:
            print(f"Warning: Sequence {name} is empty after removing non-ACGT characters. Skipping.", file=sys.stderr)
            continue

        if len(cleaned_seq) % 3 != 0:
            raise ValueError(f"Sequence {name} (cleaned length {len(cleaned_seq)}) is not a multiple of 3.")

        if first_seq_len is None:
            first_seq_len = len(cleaned_seq)
        elif len(cleaned_seq) != first_seq_len:
            raise ValueError(f"Cleaned sequence {name} (length {len(cleaned_seq)}) has a different length from the first sequence (length {first_seq_len}). Please provide aligned sequences of equal length.")
        
        validated_sequences[name] = cleaned_seq
    
    if not validated_sequences:
        raise ValueError("No valid sequences found after validation.")
    if len(validated_sequences) < 2 and len(sequences_dict) >=2 : # Check if reduction to <2 was due to validation
        raise ValueError("Less than two sequences remained after validation. Need at least two for comparison.")
    return validated_sequences

def parse_property_file(filename):
    """Parses an amino acid property file."""
    properties = {}
    header = ""
    try:
        with open(filename, 'r') as f:
            header = f.readline().strip()
            for line in f:
                parts = line.strip().split()
                if len(parts) == 2:
                    aa, value = parts[0], float(parts[1]) # Assuming property value is float
                    properties[aa] = value
                elif line.strip(): # Non-empty line that isn't two parts
                    raise ValueError(f"Malformed line in property file: {line.strip()}")
    except FileNotFoundError:
        raise
    except ValueError as e:
        raise ValueError(f"Error parsing property file {filename}: {e}")
    except Exception as e:
        raise Exception(f"Error reading property file {filename}: {e}")
    
    if not header:
        raise ValueError("Property file is empty or header is missing.")
    if not properties:
        raise ValueError("No properties found in property file.")
    return header, properties

def codon_to_aa(codon, genetic_code):
    """Converts a codon to an amino acid."""
    return genetic_code.get(codon, 'X') # X for unknown, Z for stop

def estimate_r_value(sequences_list):
    """Estimates transition/transversion ratio (R) from 3rd codon positions."""
    # This function can be identical to the one in calculate_ds_dn.py
    # For brevity, assuming it will be copied or imported if this were a real library.
    # Placeholder for now.
    num_seqs = len(sequences_list)
    if num_seqs < 2:
        print("Warning: R estimation requires at least 2 sequences.", file=sys.stderr)
        return None

    total_transitions = 0
    total_transversions = 0
    total_comparisons = 0 

    for i in range(num_seqs):
        for j in range(i + 1, num_seqs):
            seq1 = sequences_list[i]
            seq2 = sequences_list[j]

            if len(seq1) != len(seq2) or len(seq1) % 3 != 0 : # Should be caught by validation
                continue

            for k_codon in range(0, len(seq1), 3):
                base1 = seq1[k_codon + 2] # 3rd base
                base2 = seq2[k_codon + 2] 

                if base1 not in BASE_TO_NUM or base2 not in BASE_TO_NUM:
                    continue 
                if base1 == base2:
                    continue
                
                total_comparisons += 1
                b1_is_purine = base1 in ['A', 'G']
                b2_is_purine = base2 in ['A', 'G']

                if (b1_is_purine and b2_is_purine) or \
                   (not b1_is_purine and not b2_is_purine):
                    total_transitions += 1
                else:
                    total_transversions += 1
    
    if total_comparisons == 0:
        print("Warning: No differing 3rd codon sites found for R estimation.", file=sys.stderr)
        return None

    P_prime = total_transitions / total_comparisons
    Q_prime = total_transversions / total_comparisons

    val_for_s_log = 1.0 - 2.0 * P_prime - Q_prime
    val_for_v_log = 1.0 - 2.0 * Q_prime

    if val_for_s_log <= 1e-9 or val_for_v_log <= 1e-9: # Check for log(0) or log(negative)
        print("Warning: Cannot compute log for s or v in R estimation (P', Q' too high). R set to None.", file=sys.stderr)
        return None

    try:
        s_rate = 0.5 * math.log(1.0 / val_for_s_log)
        v_rate = 0.25 * math.log(1.0 / val_for_v_log)
    except ValueError: # Should be caught by check above, but as safeguard
         print("Warning: Math domain error during R estimation log calculation. R set to None.", file=sys.stderr)
         return None


    if abs(v_rate) < 1e-9: # Avoid division by zero if v_rate is effectively zero
        print("Warning: Estimated v_rate is ~0 in R estimation. R set to None.", file=sys.stderr)
        return None
            
    estimated_r = s_rate / v_rate
    
    if estimated_r < 0:
        print("Warning: Estimated R is negative. R set to None.", file=sys.stderr)
        return None
            
    return estimated_r

# --- Core Logic Placeholder Functions ---
def calculate_syn_site_fraction(codon, genetic_code, r_value):
    """
    Calculates f_i, the fraction of synonymous changes from a codon, weighted by R.
    Similar to calculate_syn_site in calculate_ds_dn.py.
    """
    original_aa = codon_to_aa(codon, genetic_code)
    if original_aa == 'Z': return 0.0 # Stop codons

    syn_weighted_paths = 0.0
    total_weighted_paths = 0.0

    for i in range(3): # Iterate over codon positions
        original_base = codon[i]
        for new_base in BASES:
            if new_base == original_base:
                continue
            
            mutated_codon_list = list(codon)
            mutated_codon_list[i] = new_base
            mutated_codon = "".join(mutated_codon_list)
            
            mutated_aa = codon_to_aa(mutated_codon, genetic_code)
            
            is_transition = False
            if (original_base in ['A', 'G'] and new_base in ['A', 'G']) or \
               (original_base in ['C', 'T'] and new_base in ['C', 'T']):
                is_transition = True
            
            weight = r_value if is_transition else 1.0
            
            if mutated_aa != 'Z': # Only consider paths not leading to a stop codon
                total_weighted_paths += weight
                if original_aa == mutated_aa:
                    syn_weighted_paths += weight
    
    return syn_weighted_paths / total_weighted_paths if total_weighted_paths > 0 else 0.0


def calculate_potential_conservative_sites(codon_str, genetic_code, aa_properties, r_value):
    """
    Calculates potential conservative sites for a single codon.
    This is the `conserv_site` logic from Perl.
    Returns the sum of weighted fractions of conservative changes from this codon.
    """
    original_aa_char = codon_to_aa(codon_str, genetic_code)
    if original_aa_char == 'Z':  # Stop codons cannot have conservative changes
        return 0.0
    
    original_property = aa_properties.get(original_aa_char)
    if original_property is None:
        # This AA from the codon is not in the property file.
        # Perl script would crash or give warnings. Here, we assume it can't be conservative.
        # Or, one could raise an error or print to stderr. For now, treat as non-conservative.
        # print(f"Warning: Amino acid {original_aa_char} from codon {codon_str} not in property file.", file=sys.stderr)
        return 0.0

    total_weighted_conservative_fraction = 0.0

    for i in range(3):  # Iterate over codon positions (0, 1, 2)
        original_base = codon_str[i]
        
        # Denominator for this site's weighting, initially R+2 (1 transition, 2 transversions)
        # These are relative weights for changes at THIS specific position i.
        # The Perl script's $flag/$flag1 logic dynamically adjusts this denominator.
        site_denominator = r_value + 2.0 
        
        temp_site_conservative_numerator = 0.0

        # Check paths to stop codons from this position to adjust site_denominator
        adjusted_site_denominator = r_value + 2.0 # Start with full denominator for this site
        for alt_base_char_for_stop_check in BASES:
            if alt_base_char_for_stop_check == original_base:
                continue
            
            temp_codon_list = list(codon_str)
            temp_codon_list[i] = alt_base_char_for_stop_check
            temp_codon_str = "".join(temp_codon_list)
            temp_aa_char = codon_to_aa(temp_codon_str, genetic_code)

            if temp_aa_char == 'Z': # If this change leads to a stop codon
                is_transition_to_stop = False
                if (original_base in ['A', 'G'] and alt_base_char_for_stop_check in ['A', 'G']) or \
                   (original_base in ['C', 'T'] and alt_base_char_for_stop_check in ['C', 'T']):
                    is_transition_to_stop = True
                
                if is_transition_to_stop:
                    adjusted_site_denominator -= r_value
                else:
                    adjusted_site_denominator -= 1.0
        
        if adjusted_site_denominator <= 1e-9: # Avoid division by zero if all paths from this site lead to stop
            continue # No valid paths from this site, so it contributes 0 to conservative sites.

        # Now calculate conservative contributions from this site
        for alt_base_char in BASES:
            if alt_base_char == original_base:
                continue

            mutated_codon_list = list(codon_str)
            mutated_codon_list[i] = alt_base_char
            mutated_codon_str = "".join(mutated_codon_list)

            mutated_aa_char = codon_to_aa(mutated_codon_str, genetic_code)

            if mutated_aa_char == 'Z':  # Paths to stop codons are not conservative
                continue
            
            if mutated_aa_char == original_aa_char: # Synonymous, not conservative non-synonymous
                continue

            mutated_property = aa_properties.get(mutated_aa_char)
            if mutated_property is None:
                # print(f"Warning: Amino acid {mutated_aa_char} from mutated codon {mutated_codon_str} not in property file.", file=sys.stderr)
                continue # Cannot be conservative if property unknown

            if mutated_property == original_property:
                # This is a non-synonymous but conservative change
                is_transition = False
                if (original_base in ['A', 'G'] and alt_base_char in ['A', 'G']) or \
                   (original_base in ['C', 'T'] and alt_base_char in ['C', 'T']):
                    is_transition = True
                
                weight = r_value if is_transition else 1.0
                temp_site_conservative_numerator += weight
        
        # The fraction for this site is its conservative numerator over its adjusted denominator
        total_weighted_conservative_fraction += (temp_site_conservative_numerator / adjusted_site_denominator)

    # The Perl script sums these fractions from each of the 3 sites.
    # Each site contributes between 0 and 1 ( (sum of weights of cons. paths at site) / (sum of weights of valid paths at site) )
    # So the total sum can be between 0 and 3.
    return total_weighted_conservative_fraction

def get_codon_diff_positions(codon1_str, codon2_str):
    """Returns a list of positions (0,1,2) where two codons differ."""
    diffs = []
    for k in range(3):
        if codon1_str[k] != codon2_str[k]:
            diffs.append(k)
    return diffs

def count_substitutions_pnc_pnr(seq1_str, seq2_str, genetic_code, aa_properties, r_value):
    """
    Main calculation function for a pair of sequences.
    Mirrors `countsubstitutions` from pnc_pnr_dist.pl.
    """
    num_codons = len(seq1_str) // 3
    
    actual_syn_subs = 0.0
    actual_nonsyn_subs = 0.0
    actual_conserv_subs = 0.0 # Non-synonymous changes that are conservative

    # Calculate SA_Nei, SB_Nei (sum of f_i for each sequence)
    # f_i = fraction of synonymous changes from codon i, weighted by R
    SA_Nei = sum(calculate_syn_site_fraction(seq1_str[k*3:(k+1)*3], genetic_code, r_value) 
                 for k in range(num_codons) if codon_to_aa(seq1_str[k*3:(k+1)*3], genetic_code) != 'Z')
    SB_Nei = sum(calculate_syn_site_fraction(seq2_str[k*3:(k+1)*3], genetic_code, r_value)
                 for k in range(num_codons) if codon_to_aa(seq2_str[k*3:(k+1)*3], genetic_code) != 'Z')

    # Calculate potential conservative sites for each sequence
    # This is sum over codons of (sum over 3 sites of (weighted conservative changes / weighted valid changes from site))
    potential_conserv_sites_A = sum(calculate_potential_conservative_sites(seq1_str[k*3:(k+1)*3], genetic_code, aa_properties, r_value)
                                    for k in range(num_codons) if codon_to_aa(seq1_str[k*3:(k+1)*3], genetic_code) != 'Z')
    potential_conserv_sites_B = sum(calculate_potential_conservative_sites(seq2_str[k*3:(k+1)*3], genetic_code, aa_properties, r_value)
                                    for k in range(num_codons) if codon_to_aa(seq2_str[k*3:(k+1)*3], genetic_code) != 'Z')
    
    # Averaged potential conservative sites for the pair
    avg_potential_total_conservative_sites = (potential_conserv_sites_A + potential_conserv_sites_B) / 2.0

    # Iterate through codon pairs to count actual substitutions
    for i in range(num_codons):
        codon1 = seq1_str[i*3:(i+1)*3]
        codon2 = seq2_str[i*3:(i+1)*3]

        if codon1 == codon2:
            continue

        aa1_char = codon_to_aa(codon1, genetic_code)
        aa2_char = codon_to_aa(codon2, genetic_code)
        
        # Skip comparison if either codon is a stop codon or cannot be mapped in property file
        if aa1_char == 'Z' or aa2_char == 'Z':
            continue
        
        prop1 = aa_properties.get(aa1_char)
        prop2 = aa_properties.get(aa2_char)
        if prop1 is None or prop2 is None: # Should not happen if AAs are standard
            # print(f"Warning: Codon pair {codon1}({aa1_char}) - {codon2}({aa2_char}) involves AA not in property file.", file=sys.stderr)
            continue


        diff_positions = get_codon_diff_positions(codon1, codon2)
        num_diffs = len(diff_positions)

        if num_diffs == 1:
            pos = diff_positions[0]
            is_transition = (codon1[pos] in ['A','G'] and codon2[pos] in ['A','G']) or \
                            (codon1[pos] in ['C','T'] and codon2[pos] in ['C','T'])
            weight = r_value if is_transition else 1.0

            if aa1_char == aa2_char: # Synonymous
                actual_syn_subs += weight
            else: # Non-synonymous
                actual_nonsyn_subs += weight
                if prop1 == prop2: # Conservative non-synonymous
                    actual_conserv_subs += weight
        
        elif num_diffs == 2:
            path_syn_sum = 0.0
            path_nonsyn_sum = 0.0
            path_conserv_sum = 0.0
            num_valid_paths = 0

            # Two intermediate codons: int_A (change at diff_pos[0]), int_B (change at diff_pos[1])
            # Path 1: codon1 -> int_A -> codon2
            # Path 2: codon1 -> int_B -> codon2
            for k_path in range(2): # k_path=0 means first intermediate by changing diff_positions[0] of codon1
                                    # k_path=1 means first intermediate by changing diff_positions[1] of codon1
                
                # Step 1: codon1 -> intermediate_codon
                intermediate_codon_list = list(codon1)
                change_pos_s1 = diff_positions[k_path]
                intermediate_codon_list[change_pos_s1] = codon2[change_pos_s1]
                intermediate_codon = "".join(intermediate_codon_list)
                
                aa_intermediate_char = codon_to_aa(intermediate_codon, genetic_code)
                if aa_intermediate_char == 'Z': continue # Path via stop codon is invalid
                
                prop_intermediate = aa_properties.get(aa_intermediate_char)
                if prop_intermediate is None: continue


                is_transition_s1 = (codon1[change_pos_s1] in ['A','G'] and intermediate_codon[change_pos_s1] in ['A','G']) or \
                                   (codon1[change_pos_s1] in ['C','T'] and intermediate_codon[change_pos_s1] in ['C','T'])
                weight1 = r_value if is_transition_s1 else 1.0
                
                current_path_syn = 0.0
                current_path_nonsyn = 0.0
                current_path_conserv = 0.0

                if aa1_char == aa_intermediate_char:
                    current_path_syn += weight1
                else:
                    current_path_nonsyn += weight1
                    if prop1 == prop_intermediate:
                        current_path_conserv += weight1
                
                # Step 2: intermediate_codon -> codon2
                change_pos_s2 = diff_positions[1-k_path] # The other differing position
                is_transition_s2 = (intermediate_codon[change_pos_s2] in ['A','G'] and codon2[change_pos_s2] in ['A','G']) or \
                                   (intermediate_codon[change_pos_s2] in ['C','T'] and codon2[change_pos_s2] in ['C','T'])
                weight2 = r_value if is_transition_s2 else 1.0

                if aa_intermediate_char == aa2_char:
                    current_path_syn += weight2
                else:
                    current_path_nonsyn += weight2
                    if prop_intermediate == prop2:
                        current_path_conserv += weight2
                
                num_valid_paths += 1
                path_syn_sum += current_path_syn
                path_nonsyn_sum += current_path_nonsyn
                path_conserv_sum += current_path_conserv
            
            if num_valid_paths > 0:
                actual_syn_subs += path_syn_sum / num_valid_paths
                actual_nonsyn_subs += path_nonsyn_sum / num_valid_paths
                actual_conserv_subs += path_conserv_sum / num_valid_paths

        elif num_diffs == 3:
            path_syn_sum = 0.0
            path_nonsyn_sum = 0.0
            path_conserv_sum = 0.0
            num_valid_paths = 0

            for p_order in itertools.permutations(diff_positions): # p_order is like (pos0, pos1, pos2)
                # Path: codon1 -> int1 (change at p_order[0]) -> int2 (change at p_order[1]) -> codon2 (change at p_order[2])
                
                # Step 1: codon1 -> int1
                int1_list = list(codon1)
                change_pos_s1 = p_order[0]
                int1_list[change_pos_s1] = codon2[change_pos_s1]
                int1_str = "".join(int1_list)
                aa_int1_char = codon_to_aa(int1_str, genetic_code)
                if aa_int1_char == 'Z': continue
                prop_int1 = aa_properties.get(aa_int1_char)
                if prop_int1 is None: continue
                
                is_trans_s1 = (codon1[change_pos_s1] in ['A','G'] and int1_str[change_pos_s1] in ['A','G']) or \
                              (codon1[change_pos_s1] in ['C','T'] and int1_str[change_pos_s1] in ['C','T'])
                w1 = r_value if is_trans_s1 else 1.0
                
                current_path_syn = 0.0
                current_path_nonsyn = 0.0
                current_path_conserv = 0.0

                if aa1_char == aa_int1_char: current_path_syn += w1
                else:
                    current_path_nonsyn += w1
                    if prop1 == prop_int1: current_path_conserv += w1

                # Step 2: int1 -> int2
                int2_list = list(int1_str)
                change_pos_s2 = p_order[1]
                int2_list[change_pos_s2] = codon2[change_pos_s2]
                int2_str = "".join(int2_list)
                aa_int2_char = codon_to_aa(int2_str, genetic_code)
                if aa_int2_char == 'Z': continue
                prop_int2 = aa_properties.get(aa_int2_char)
                if prop_int2 is None: continue

                is_trans_s2 = (int1_str[change_pos_s2] in ['A','G'] and int2_str[change_pos_s2] in ['A','G']) or \
                              (int1_str[change_pos_s2] in ['C','T'] and int2_str[change_pos_s2] in ['C','T'])
                w2 = r_value if is_trans_s2 else 1.0

                if aa_int1_char == aa_int2_char: current_path_syn += w2
                else:
                    current_path_nonsyn += w2
                    if prop_int1 == prop_int2: current_path_conserv += w2
                
                # Step 3: int2 -> codon2
                change_pos_s3 = p_order[2]
                is_trans_s3 = (int2_str[change_pos_s3] in ['A','G'] and codon2[change_pos_s3] in ['A','G']) or \
                              (int2_str[change_pos_s3] in ['C','T'] and codon2[change_pos_s3] in ['C','T'])
                w3 = r_value if is_trans_s3 else 1.0

                if aa_int2_char == aa2_char: current_path_syn += w3
                else:
                    current_path_nonsyn += w3
                    if prop_int2 == prop2: current_path_conserv += w3
                
                num_valid_paths +=1
                path_syn_sum += current_path_syn
                path_nonsyn_sum += current_path_nonsyn
                path_conserv_sum += current_path_conserv

            if num_valid_paths > 0:
                actual_syn_subs += path_syn_sum / num_valid_paths
                actual_nonsyn_subs += path_nonsyn_sum / num_valid_paths
                actual_conserv_subs += path_conserv_sum / num_valid_paths

    # Final calculations for pNc, pNr and their SEs
    potential_syn_total = (SA_Nei + SB_Nei) / 2.0
    
    # num_codons_for_nonsyn is total codons excluding pairs where either is stop.
    # The SA_Nei/SB_Nei already exclude contributions from initial stop codons.
    # The iteration for actual subs also skips pairs with stop codons.
    # So num_codons is a good proxy for total comparable codons.
    potential_nonsyn_total = (num_codons * 3.0) - potential_syn_total
    if potential_nonsyn_total < 0: potential_nonsyn_total = 0 # Should not happen

    # avg_potential_total_conservative_sites is already calculated.
    # This is $ConSite in Perl. $NonSite = $N_sites - $ConSite.
    potential_rad_sites = potential_nonsyn_total - avg_potential_total_conservative_sites
    if potential_rad_sites < 0: potential_rad_sites = 0

    actual_rad_subs = actual_nonsyn_subs - actual_conserv_subs
    if actual_rad_subs < 0: actual_rad_subs = 0 # Ensure non-negative due to float precision

    pnc, SEpnc, pnr, SEpnr = "NA", "NA", "NA", "NA"

    if avg_potential_total_conservative_sites > 1e-9:
        val_pnc = actual_conserv_subs / avg_potential_total_conservative_sites
        if 0 <= val_pnc <= 1.0: # pNc must be a proportion
             pnc = val_pnc
             SEpnc = math.sqrt(pnc * (1.0 - pnc) / avg_potential_total_conservative_sites) if pnc * (1.0 - pnc) >= -1e-9 else "NA" # allow for small neg due to precision
        elif val_pnc > 1.0 and val_pnc < 1.0 + 1e-6: # Allow for small float inaccuracies if slightly above 1
            pnc = 1.0
            SEpnc = 0.0 # p(1-p) = 0
        else: # val_pnc is significantly > 1 or < 0
            # print(f"Warning: pNc calculation out of bounds ({val_pnc}). Setting to NA.", file=sys.stderr)
            pnc = "NA"
            SEpnc = "NA"
    
    if potential_rad_sites > 1e-9:
        val_pnr = actual_rad_subs / potential_rad_sites
        if 0 <= val_pnr <= 1.0: # pNr must be a proportion
            pnr = val_pnr
            SEpnr = math.sqrt(pnr * (1.0 - pnr) / potential_rad_sites) if pnr * (1.0 - pnr) >= -1e-9 else "NA"
        elif val_pnr > 1.0 and val_pnr < 1.0 + 1e-6:
            pnr = 1.0
            SEpnr = 0.0
        else:
            # print(f"Warning: pNr calculation out of bounds ({val_pnr}). Setting to NA.", file=sys.stderr)
            pnr = "NA"
            SEpnr = "NA"
            
    # The dS/dN calculations are not requested for the final output of this script.
    # ps = actual_syn_subs / potential_syn_total if potential_syn_total > 1e-9 else 0
    # pn = actual_nonsyn_subs / potential_nonsyn_total if potential_nonsyn_total > 1e-9 else 0
    # ds, dn calculations would follow here if needed.
    
    return pnc, SEpnc, pnr, SEpnr


# --- Main Execution ---
def main():
    parser = argparse.ArgumentParser(description="Calculate pNc/pNr ratios for sequences in a FASTA file.")
    parser.add_argument("input_file", help="Path to the input FASTA file.")
    parser.add_argument("property_file", help="Path to the amino acid property file.")
    parser.add_argument("ratio", nargs='?', default="0.5", help="Transition/transversion ratio (R). Use 'R' to estimate from data or provide a float (default: 0.5).")
    
    args = parser.parse_args()

    try:
        raw_sequences = parse_fasta(args.input_file)
        sequences = validate_sequences(raw_sequences)
        
        if len(sequences) < 2:
            print("Error: At least two valid sequences are required for comparison.", file=sys.stderr)
            sys.exit(1)
            
        property_header, aa_properties = parse_property_file(args.property_file)

    except FileNotFoundError as e:
        print(f"Error: File not found. {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e: # Catch any other unexpected error during setup
        print(f"An unexpected error occurred during setup: {e}", file=sys.stderr)
        sys.exit(1)

    r_value_to_use = None
    r_estimation_failed_globally = False

    if isinstance(args.ratio, str) and args.ratio.upper() == 'R':
        estimated_r = estimate_r_value(list(sequences.values()))
        if estimated_r is None or estimated_r <= 0:
            print("Warning: R estimation resulted in an undefined or non-positive value. 'undef' will be used for the ratio in output.", file=sys.stderr)
            r_estimation_failed_globally = True
        else:
            r_value_to_use = estimated_r
    else:
        try:
            r_value_to_use = float(args.ratio)
            if r_value_to_use <= 0:
                print("Error: Ratio must be a positive float or 'R'.", file=sys.stderr)
                sys.exit(1)
        except ValueError:
            print(f"Error: Invalid ratio value '{args.ratio}'. Must be a positive float or 'R'.", file=sys.stderr)
            sys.exit(1)

    seq_names = list(sequences.keys())
    print(f"Seq1\tSeq2\tProperty\tRatio\tpNc\tSEpNc\tpNr\tSEpNr")


    for i in range(len(seq_names)):
        for j in range(i + 1, len(seq_names)):
            seq1_name = seq_names[i]
            seq2_name = seq_names[j]
            seq1 = sequences[seq1_name]
            seq2 = sequences[seq2_name]

            formatted_ratio_str = "undef"
            current_r_for_calc = None

            if r_estimation_failed_globally:
                # Output 'undef' for ratio, other values will be NA as calculations depend on R.
                # The count_substitutions_pnc_pnr must handle None R if we want NA for pNc etc.
                # Or, we can assign a default R for calculation if global estimation fails but proceed.
                # The Perl script seems to print "undef" for ratio and then likely stops or produces NA.
                # For now, if R is undef globally, we'll print undef for ratio and subsequent NA values.
                 pass # formatted_ratio_str is already "undef"
            elif r_value_to_use is not None:
                formatted_ratio_str = f"{r_value_to_use:.2f}"
                current_r_for_calc = r_value_to_use
            else: # Should not happen if logic is correct
                print(f"Internal Error: R value ambiguous for pair {seq1_name}-{seq2_name}.", file=sys.stderr)
                # formatted_ratio_str remains "undef"

            if current_r_for_calc is None and not r_estimation_failed_globally:
                # This implies a fixed ratio was not parseable, but that's caught earlier.
                # Or 'R' was specified, but it's not globally failed yet not set.
                # This state indicates an issue, likely means we should output undef/NA and continue.
                 print(f"{seq1_name}\t{seq2_name}\t{property_header}\t{formatted_ratio_str}\tNA\tNA\tNA\tNA")
                 continue


            # If R estimation failed globally, current_r_for_calc will be None.
            # The count_substitutions_pnc_pnr and subsequent calculations need to handle this
            # to produce "NA" for pNc, SEpNc, pNr, SEpNr.
            # The Perl script's `countsubstitutions` uses $ratio directly. If it's "undef",
            # it would fail if not checked. The Perl script checks `if ($ratio eq "undef")`
            # at the top of the loop for pairs.
            if r_estimation_failed_globally:
                 print(f"{seq1_name}\t{seq2_name}\t{property_header}\t{formatted_ratio_str}\tNA\tNA\tNA\tNA")
                 continue


            pnc, SEpnc, pnr, SEpnr = count_substitutions_pnc_pnr(
                seq1, seq2, STANDARD_GENETIC_CODE, aa_properties, current_r_for_calc
            )
            
            # Format output values
            pnc_str = f"{pnc:.4f}" if isinstance(pnc, float) else "NA"
            SEpnc_str = f"{SEpnc:.4f}" if isinstance(SEpnc, float) else "NA"
            pnr_str = f"{pnr:.4f}" if isinstance(pnr, float) else "NA"
            SEpnr_str = f"{SEpnr:.4f}" if isinstance(SEpnr, float) else "NA"

            print(f"{seq1_name}\t{seq2_name}\t{property_header}\t{formatted_ratio_str}\t{pnc_str}\t{SEpnc_str}\t{pnr_str}\t{SEpnr_str}")

if __name__ == "__main__":
    main()
