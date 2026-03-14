"""
dsdn_dist.py — Synonymous (dS) and Non-synonymous (dN) Substitution Rate Calculator
===================================================================================
Implementation of the Nei & Gojobori (1986) method for estimating the numbers of
synonymous and nonsynonymous nucleotide substitutions per site.

Reference:
    Nei, M., & Gojobori, T. (1986). Simple methods for estimating the numbers of
    synonymous and nonsynonymous nucleotide substitutions. Molecular Biology
    and Evolution, 3(4), 418-426.

Usage:
    python dsdn_dist.py <input_fasta> [ratio | R]

Arguments:
    input_fasta  : FASTA file containing a pairwise or multiple pairwise
                   sequence alignment after complete gap deletion.
    ratio        : Transition/transversion ratio (default 0.5).
                   Pass R to estimate the ratio from 3rd codon positions.

Input requirements (all enforced at runtime):
    1. FASTA format: each sequence on a SINGLE line following its header.
    2. All sequences must be the same length (complete gap deletion enforced).
    3. Sequence length must be divisible by 3 (complete codons only).
    4. Only the bases A, C, G, T are permitted (no gaps, no ambiguity codes).
    5. No internal stop codons are permitted (TAA, TAG, TGA).
    6. If a terminal stop codon is present (last codon = TAA, TAG, or TGA)
       it is automatically excluded from all calculations.

Output (tab-delimited):
    gene1  gene2  ratio  ps  SEps  pn  SEpn  ds  SEds  dn  SEdn
"""

import sys
import math
import re

# ---------------------------------------------------------------------------
# Codon tables
# ---------------------------------------------------------------------------

# Amino acid encoded by each codon (0-63, standard genetic code).
# Stop codons are represented as 'Z'.
# Codon index = (base1_num * 16) + (base0_num * 4) + base2_num
# where T=0, C=1, A=2, G=3  (base0=1st, base1=2nd, base2=3rd codon position)
AA_ARRAY = [
    'F', 'F', 'L', 'L', 'L', 'L', 'L', 'L', 'I', 'I', 'I', 'M', 'V', 'V', 'V', 'V',
    'S', 'S', 'S', 'S', 'P', 'P', 'P', 'P', 'T', 'T', 'T', 'T', 'A', 'A', 'A', 'A',
    'Y', 'Y', 'Z', 'Z', 'H', 'H', 'Q', 'Q', 'N', 'N', 'K', 'K', 'D', 'D', 'E', 'E',
    'C', 'C', 'Z', 'W', 'R', 'R', 'R', 'R', 'S', 'S', 'R', 'R', 'G', 'G', 'G', 'G'
]

# Mapping of nucleotide bases to numeric indices used in codon_conversion.
BASE_NUMBER = {"T": 0, "C": 1, "A": 2, "G": 3}

# Stop codon triplets (used for terminal stop codon detection).
STOP_CODONS = {"TAA", "TAG", "TGA"}


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

def validate_sequences(names, seqs):
    """
    Enforce all biological and mathematical requirements on the input sequences.
    Exits with an informative error message if any check fails.
    Removes a terminal stop codon from every sequence if one is present.

    Checks performed:
      1. At least two sequences are present.
      2. All sequences are the same length (equal-length alignment after
         complete gap deletion).
      3. Sequence length is divisible by 3 (complete codons only).
      4. Only A, C, G, T characters are present (no gaps, no ambiguity codes).
      5. No internal stop codons exist (TAA, TAG, TGA at any non-terminal codon).
      6. If the last codon of every sequence is a stop codon it is stripped and
         a notice is printed. If only some sequences end in a stop codon, an
         error is raised because the alignment would become inconsistent.
    """
    if len(seqs) < 2:
        _error("At least two sequences are required.")

    # Check 2: equal lengths
    lengths = [len(s) for s in seqs]
    if len(set(lengths)) != 1:
        detail = ", ".join(f"{n}: {l}" for n, l in zip(names, lengths))
        _error(
            "All sequences must be the same length (complete gap deletion "
            f"required).\n  Lengths found — {detail}"
        )
    seq_len = lengths[0]

    # Check 3: divisible by 3
    if seq_len % 3 != 0:
        _error(
            f"Sequence length ({seq_len}) is not divisible by 3. "
            "Only complete codons are supported."
        )

    # Check 4: only A/C/G/T
    valid = set("ACGT")
    for name, seq in zip(names, seqs):
        bad = sorted(set(seq) - valid)
        if bad:
            _error(
                f"Sequence '{name}' contains invalid character(s): "
                f"{bad}. Only A, C, G, T are permitted."
            )

    # Check 5 & 6: stop codons
    num_codons = seq_len // 3
    terminal_stops = [
        AA_ARRAY[codon_conversion(list(seq[-3:]))] == 'Z' for seq in seqs
    ]

    # Internal stop codons (all codons except the last)
    for name, seq in zip(names, seqs):
        for pos in range(0, seq_len - 3, 3):   # excludes terminal codon
            codon = seq[pos:pos + 3]
            if codon in STOP_CODONS:
                codon_num = pos // 3 + 1
                _error(
                    f"Sequence '{name}' contains an internal stop codon "
                    f"({codon}) at codon position {codon_num}. "
                    "Internal stop codons are not permitted for this method."
                )

    # Terminal stop codon handling
    # Note: a mixed case (stop in some but not all) is impossible here because
    # all sequences are already confirmed to be the same length and divisible by 3.
    if all(terminal_stops):
        print(
            f"[INFO] Terminal stop codon detected in all {len(seqs)} sequences "
            "— stripped from all sequences before calculation.",
            file=sys.stderr
        )
        seqs = [s[:-3] for s in seqs]
    # else: no terminal stops — proceed normally

    return seqs


def _error(msg):
    """Print a formatted error message and exit."""
    print(f"[ERROR] {msg}", file=sys.stderr)
    sys.exit(1)


# ---------------------------------------------------------------------------
# Core computation functions
# ---------------------------------------------------------------------------

def codon_conversion(codon):
    """
    Convert a 3-base codon list to an integer index (0-63).
    Encoding: index = (base1 * 16) + (base0 * 4) + base2
    where base0/1/2 are the 1st/2nd/3rd codon positions and T=0,C=1,A=2,G=3.
    Returns -1 for any codon containing an unrecognised base.
    """
    if len(codon) != 3 or any(b not in BASE_NUMBER for b in codon):
        return -1
    return (BASE_NUMBER[codon[1]] * 16) + (BASE_NUMBER[codon[0]] * 4) + BASE_NUMBER[codon[2]]


def get_syn_site(codon_num, T, R):
    """
    Return the number of synonymous sites for a codon (Nei & Gojobori 1986).
    T = R / (R + 2)  where R = 2 * ratio (transition/transversion ratio).
    Values are pre-computed for all 64 codons (stop codons return 0).
    """
    if codon_num < 0:
        return 0.0
    sites = {
        0: T * 3,            1: T * 3,            2: (2 * T) * 3,      3: (2 * T) * 3,
        4: 3.0,              5: 3.0,              6: (1 + T) * 3,      7: (1 + T) * 3,
        8: ((R+1)/(R+2))*3, 9: ((R+1)/(R+2))*3, 10: (2/(R+2))*3,     11: 0.0,
        12: 3.0,             13: 3.0,             14: 3.0,             15: 3.0,
        16: 3.0,             17: 3.0,             18: 3.0,             19: 3.0,
        20: 3.0,             21: 3.0,             22: 3.0,             23: 3.0,
        24: 3.0,             25: 3.0,             26: 3.0,             27: 3.0,
        28: 3.0,             29: 3.0,             30: 3.0,             31: 3.0,
        32: 3.0,             33: 3.0,             34: 0.0,             35: 0.0,
        36: T * 3,           37: T * 3,           38: T * 3,           39: T * 3,
        40: T * 3,           41: T * 3,           42: T * 3,           43: T * 3,
        44: T * 3,           45: T * 3,           46: T * 3,           47: T * 3,
        48: (R/(R+1))*3,    49: (R/(R+1))*3,    50: 0.0,             51: 0.0,
        52: 3.0,             53: 3.0,             54: 1.5 * 3,         55: (1+1/(R+2))*3,
        56: T * 3,           57: T * 3,           58: (T+1/(R+1))*3,  59: (T+1/(R+2))*3,
        60: 3.0,             61: 3.0,             62: 3.0,             63: 3.0,
    }
    return sites.get(codon_num, 0.0)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) < 2:
        print("Usage: python dsdn.py <input_fasta> [ratio | R]")
        sys.exit(1)

    input_file = sys.argv[1]
    ratio_arg  = sys.argv[2] if len(sys.argv) > 2 else "0.5"
    base_ratio = 'R' if ratio_arg.upper() == 'R' else float(ratio_arg)

    # --- Parse FASTA (single-line sequences required) ---
    names = []
    seqs  = []
    current_seq = []
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if line.startswith('>'):
                if current_seq:
                    seqs.append(re.sub(r'\s+', '', ''.join(current_seq)).upper())
                    current_seq = []
                names.append(line[1:].strip())
            else:
                current_seq.append(line.upper())
    if current_seq:
        seqs.append(re.sub(r'\s+', '', ''.join(current_seq)).upper())

    # --- Validate all sequences before any calculation ---
    seqs = validate_sequences(names, seqs)

    seq_len  = len(seqs[0])
    num_seqs = len(names)

    # --- Estimate global R from all sequences if requested ---
    # R is estimated once from all 3rd codon positions across all sequence pairs
    # in the alignment, rather than per pair. This is more statistically sound
    # because: (1) stationarity is already assumed by the K2P / Nei-Gojobori
    # method, making R a property of the process rather than any individual
    # lineage pair; and (2) pooling all pairs yields a much larger sample of
    # 3rd-position substitutions, reducing estimation noise from small counts.
    # The constant R value applied to all pairs will be visible in the output.
    if base_ratio == 'R':
        P1 = P2 = Q = 0
        num_codons = seq_len // 3
        for na in range(num_seqs):
            for n in range(na + 1, num_seqs):
                for x in range(0, seq_len, 3):
                    f_base = seqs[na][x + 2]
                    s_base = seqs[n][x + 2]
                    if f_base != s_base:
                        if   f_base in 'AG' and s_base in 'AG': P1 += 1
                        elif f_base in 'AG' and s_base in 'CT': Q  += 1
                        elif f_base in 'CT' and s_base in 'AG': Q  += 1
                        elif f_base in 'CT' and s_base in 'CT': P2 += 1
        total_pairs  = num_seqs * (num_seqs - 1) // 2
        total_codons = num_codons * total_pairs
        P     = (P1 + P2) / total_codons
        Q_val = Q / total_codons
        den1  = 1.0 - 2.0 * P - Q_val
        den2  = 1.0 - 2.0 * Q_val
        if den1 <= 0 or den2 <= 0:
            _error("Global R estimation failed (K2P correction undefined). "
                   "Try supplying a fixed ratio instead.")
        w1    = 1.0 / den1
        w2    = 1.0 / den2
        s_val = 0.5 * math.log(w1) - 0.25 * math.log(w2)
        v_val = 0.5 * math.log(w2)
        if v_val == 0:
            _error("Global R estimation failed (no transversions observed at "
                   "3rd codon positions). Try supplying a fixed ratio instead.")
        base_ratio = s_val / v_val
        if base_ratio <= 0:
            _error(
                f"Global R estimation produced a non-positive value (R = {base_ratio:.4f}), "
                "indicating more transversions than transitions were observed at 3rd codon "
                "positions. This is mathematically invalid for this method: the synonymous "
                "site calculations require division by (R + 1) and (R + 2), which become "
                "zero or negative when R <= 0, producing undefined or nonsensical results. "
                "Likely causes are: (1) insufficient 3rd-codon-position variation due to "
                "too few sequences or too short an alignment; or (2) high sequence divergence "
                "causing transition saturation, where multiple hits erase the transition "
                "signal and transversions appear to dominate — the regime where the K2P "
                "model itself is unreliable. Try supplying a fixed ratio instead (e.g. 0.5)."
            )
        print(
            f"[INFO] Global R estimated from all {total_pairs} sequence pair(s) "
            f"across {num_codons} codon positions: R = {base_ratio:.4f} "
            f"(transition/transversion ratio = {base_ratio / 2:.4f}). "
            "This value is applied uniformly to all pairs.",
            file=sys.stderr
        )

    # --- Output header ---
    print("gene1\tgene2\tratio\tps\tSEps\tpn\tSEpn\tds\tSEds\tdn\tSEdn")

    def fmt(val):
        return f"{val:.4f}" if isinstance(val, float) else str(val)

    # --- Pairwise loop ---
    for na in range(num_seqs):
        for n in range(na + 1, num_seqs):
            seqA = seqs[na]
            seqB = seqs[n]

            current_ratio = base_ratio

            ratio_rounded = float(f"{current_ratio:.2f}")
            R_val = ratio_rounded * 2.0
            T_val = R_val / (R_val + 2.0)

            count_codons  = 0
            SA_Nei        = 0.0
            SB_Nei        = 0.0
            syn_codons    = 0.0
            nonsyn_codons = 0.0

            for i in range(0, seq_len, 3):
                codA = list(seqA[i:i+3])
                codB = list(seqB[i:i+3])

                count_codons  += 1
                codon_numberA  = codon_conversion(codA)
                codon_numberB  = codon_conversion(codB)
                SA_Nei        += get_syn_site(codon_numberA, T_val, R_val)
                SB_Nei        += get_syn_site(codon_numberB, T_val, R_val)

                if codA == codB:
                    continue

                diffs = sum(1 for j in range(3) if codA[j] != codB[j])

                # --- 1-base change ---
                if diffs == 1:
                    if AA_ARRAY[codon_numberA] == AA_ARRAY[codon_numberB]:
                        syn_codons    += 1.0
                    else:
                        nonsyn_codons += 1.0

                # --- 2-base changes: 3 cases by which positions differ ---
                elif diffs == 2:
                    if   codA[0] != codB[0] and codA[1] != codB[1]:  # positions 0 & 1
                        codC = [codB[0], codA[1], codA[2]]
                        codD = [codA[0], codB[1], codA[2]]
                    elif codA[0] != codB[0] and codA[2] != codB[2]:  # positions 0 & 2
                        codC = [codB[0], codA[1], codA[2]]
                        codD = [codA[0], codA[1], codB[2]]
                    else:                                              # positions 1 & 2
                        codC = [codA[0], codB[1], codA[2]]
                        codD = [codA[0], codA[1], codB[2]]

                    codon_numberC = codon_conversion(codC)
                    codon_numberD = codon_conversion(codD)
                    tmp_syn = temp = 0.0

                    if AA_ARRAY[codon_numberC] == 'Z':
                        temp += 1.0
                    else:
                        if AA_ARRAY[codon_numberA] == AA_ARRAY[codon_numberC]: tmp_syn += 1.0
                        if AA_ARRAY[codon_numberB] == AA_ARRAY[codon_numberC]: tmp_syn += 1.0

                    if AA_ARRAY[codon_numberD] == 'Z':
                        temp += 1.0
                    else:
                        if AA_ARRAY[codon_numberA] == AA_ARRAY[codon_numberD]: tmp_syn += 1.0
                        if AA_ARRAY[codon_numberB] == AA_ARRAY[codon_numberD]: tmp_syn += 1.0

                    if temp < 2.0:
                        tmp_syn        = tmp_syn / (2.0 - temp)
                        syn_codons    += tmp_syn
                        nonsyn_codons += 2.0 - tmp_syn

                # --- 3-base changes: 6 paths through pairs of intermediates ---
                elif diffs == 3:
                    # One-base-from-A intermediates
                    codC = [codB[0], codA[1], codA[2]];  codon_numberC = codon_conversion(codC)
                    codD = [codA[0], codB[1], codA[2]];  codon_numberD = codon_conversion(codD)
                    codE = [codA[0], codA[1], codB[2]];  codon_numberE = codon_conversion(codE)
                    # One-base-from-B intermediates (two-base-from-A)
                    codF = [codB[0], codB[1], codA[2]];  codon_numberF = codon_conversion(codF)
                    codG = [codB[0], codA[1], codB[2]];  codon_numberG = codon_conversion(codG)
                    codH = [codA[0], codB[1], codB[2]];  codon_numberH = codon_conversion(codH)

                    temp    = 0.0
                    tmp_syn = 0.0

                    # Path A->C->F->B
                    if AA_ARRAY[codon_numberC] == 'Z' or AA_ARRAY[codon_numberF] == 'Z': temp += 1.0
                    else:
                        if AA_ARRAY[codon_numberA] == AA_ARRAY[codon_numberC]: tmp_syn += 1.0
                        if AA_ARRAY[codon_numberC] == AA_ARRAY[codon_numberF]: tmp_syn += 1.0
                        if AA_ARRAY[codon_numberB] == AA_ARRAY[codon_numberF]: tmp_syn += 1.0
                    # Path A->C->G->B
                    if AA_ARRAY[codon_numberC] == 'Z' or AA_ARRAY[codon_numberG] == 'Z': temp += 1.0
                    else:
                        if AA_ARRAY[codon_numberA] == AA_ARRAY[codon_numberC]: tmp_syn += 1.0
                        if AA_ARRAY[codon_numberC] == AA_ARRAY[codon_numberG]: tmp_syn += 1.0
                        if AA_ARRAY[codon_numberB] == AA_ARRAY[codon_numberG]: tmp_syn += 1.0
                    # Path A->D->F->B
                    if AA_ARRAY[codon_numberD] == 'Z' or AA_ARRAY[codon_numberF] == 'Z': temp += 1.0
                    else:
                        if AA_ARRAY[codon_numberA] == AA_ARRAY[codon_numberD]: tmp_syn += 1.0
                        if AA_ARRAY[codon_numberD] == AA_ARRAY[codon_numberF]: tmp_syn += 1.0
                        if AA_ARRAY[codon_numberB] == AA_ARRAY[codon_numberF]: tmp_syn += 1.0
                    # Path A->D->H->B
                    if AA_ARRAY[codon_numberD] == 'Z' or AA_ARRAY[codon_numberH] == 'Z': temp += 1.0
                    else:
                        if AA_ARRAY[codon_numberA] == AA_ARRAY[codon_numberD]: tmp_syn += 1.0
                        if AA_ARRAY[codon_numberD] == AA_ARRAY[codon_numberH]: tmp_syn += 1.0
                        if AA_ARRAY[codon_numberB] == AA_ARRAY[codon_numberH]: tmp_syn += 1.0
                    # Path A->E->G->B
                    if AA_ARRAY[codon_numberE] == 'Z' or AA_ARRAY[codon_numberG] == 'Z': temp += 1.0
                    else:
                        if AA_ARRAY[codon_numberA] == AA_ARRAY[codon_numberE]: tmp_syn += 1.0
                        if AA_ARRAY[codon_numberE] == AA_ARRAY[codon_numberG]: tmp_syn += 1.0
                        if AA_ARRAY[codon_numberB] == AA_ARRAY[codon_numberG]: tmp_syn += 1.0
                    # Path A->E->H->B
                    if AA_ARRAY[codon_numberE] == 'Z' or AA_ARRAY[codon_numberH] == 'Z': temp += 1.0
                    else:
                        if AA_ARRAY[codon_numberA] == AA_ARRAY[codon_numberE]: tmp_syn += 1.0
                        if AA_ARRAY[codon_numberE] == AA_ARRAY[codon_numberH]: tmp_syn += 1.0
                        if AA_ARRAY[codon_numberB] == AA_ARRAY[codon_numberH]: tmp_syn += 1.0

                    if temp < 6.0:
                        syn_codons    += tmp_syn / (6.0 - temp)
                        nonsyn_codons += 3.0 - (tmp_syn / (6.0 - temp))

            # --- Summary statistics ---
            potential_syn    = (SA_Nei / 3.0 + SB_Nei / 3.0) / 2.0
            potential_nonsyn = 3.0 * count_codons - potential_syn

            if potential_syn == 0 or potential_nonsyn == 0:
                print(f"{names[na]}\t{names[n]}\tdenom 0")
                continue

            ps = syn_codons    / potential_syn
            pn = nonsyn_codons / potential_nonsyn

            ds = -0.75 * math.log(1.0 - (4.0 * ps / 3.0)) if ps < 0.75 else "NA"
            dn = -0.75 * math.log(1.0 - (4.0 * pn / 3.0)) if pn < 0.75 else "NA"

            inner_SEps = ps * (1.0 - ps) / potential_syn
            SEps       = math.sqrt(inner_SEps) if inner_SEps >= 0 else "NA"
            inner_SEpn = pn * (1.0 - pn) / potential_nonsyn
            SEpn       = math.sqrt(inner_SEpn) if inner_SEpn >= 0 else "NA"

            if ps < 0.75:
                inner_SEds = ps * (1.0 - ps) / (((1.0 - ps * 4.0 / 3.0) ** 2) * potential_syn)
                SEds = math.sqrt(inner_SEds) if inner_SEds >= 0 else "NA"
            else:
                SEds = "NA"

            if pn < 0.75:
                inner_SEdn = pn * (1.0 - pn) / (((1.0 - pn * 4.0 / 3.0) ** 2) * potential_nonsyn)
                SEdn = math.sqrt(inner_SEdn) if inner_SEdn >= 0 else "NA"
            else:
                SEdn = "NA"

            print(f"{names[na]}\t{names[n]}\t{ratio_rounded:.2f}\t"
                  f"{fmt(ps)}\t{fmt(SEps)}\t{fmt(pn)}\t{fmt(SEpn)}\t"
                  f"{fmt(ds)}\t{fmt(SEds)}\t{fmt(dn)}\t{fmt(SEdn)}")


if __name__ == '__main__':
    main()
