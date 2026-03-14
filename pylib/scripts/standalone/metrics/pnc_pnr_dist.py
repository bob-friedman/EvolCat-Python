"""
pnc_pnr_dist.py — Conservative (pNC) and Radical (pNR) Substitution Rate Calculator
===================================================================================
Implementation of the method by Hughes, Ota, and Nei (1990) to test for positive
selection by comparing the rates of conservative and radical nonsynonymous
substitutions with respect to amino acid properties (e.g., charge).

This method extends the Nei & Gojobori (1986) framework by categorizing
nonsynonymous sites and substitutions as either 'conservative' or 'radical'
based on a provided amino acid property classification.

Reference:
    Hughes, A. L., Ota, T., & Nei, M. (1990). Positive Darwinian selection
    promotes charge profile diversity in the antigen-binding cleft of class I
    major-histocompatibility-complex molecules. Molecular Biology and
    Evolution, 7(6), 515-524.

Usage:
    python pnc_pnr_dist.py <input_fasta> <property_file> [ratio | R]

Arguments:
    input_fasta   : FASTA file containing a pairwise or multiple pairwise
                    sequence alignment after complete gap deletion.
    property_file : Tab-delimited file defining amino acid categories.
                    A sample 'property_charge' file is provided in test_data.
    ratio         : Transition/transversion ratio (default 0.5).
                    Pass R to estimate the ratio from 3rd codon positions.

Input requirements:
    1. FASTA format: each sequence on a SINGLE line following its header.
    2. All sequences must be the same length (complete gap deletion enforced).
    3. Sequence length must be divisible by 3 (complete codons only).
    4. Only the bases A, C, G, T are permitted.
    5. No internal stop codons are permitted.

Output (tab-delimited):
    gene1  gene2  ratio  pnc  SE_pnc  pnr  SE_pnr
    (where pnc is the proportion of conservative nonsynonymous substitutions
    and pnr is the proportion of radical nonsynonymous substitutions).
"""

import sys
import math
import re

# default ratio is 0.5; pass R to estimate from 3rd codon positions

AA_NUMBER = [
    1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 5, 5, 5, 5,
    6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9,
    10, 10, 21, 21, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16,
    17, 17, 21, 18, 19, 19, 19, 19, 6, 6, 19, 19, 20, 20, 20, 20
]

AA_ARRAY = [
    'F', 'F', 'L', 'L', 'L', 'L', 'L', 'L', 'I', 'I', 'I', 'M', 'V', 'V', 'V', 'V',
    'S', 'S', 'S', 'S', 'P', 'P', 'P', 'P', 'T', 'T', 'T', 'T', 'A', 'A', 'A', 'A',
    'Y', 'Y', 'Z', 'Z', 'H', 'H', 'Q', 'Q', 'N', 'N', 'K', 'K', 'D', 'D', 'E', 'E',
    'C', 'C', 'Z', 'W', 'R', 'R', 'R', 'R', 'S', 'S', 'R', 'R', 'G', 'G', 'G', 'G'
]

BASE_NUMBER = {"T": 0, "C": 1, "A": 2, "G": 3}


def codon_conversion(codon):
    return (BASE_NUMBER[codon[1]] * 16) + (BASE_NUMBER[codon[0]] * 4) + BASE_NUMBER[codon[2]]


def syn_site(codon, T, R):
    cn = codon_conversion(codon)
    sites = [
        T*3, T*3, (2*T)*3, (2*T)*3, 3.0, 3.0, (1+T)*3, (1+T)*3,
        ((R+1)/(R+2))*3, ((R+1)/(R+2))*3, (2/(R+2))*3, 0.0,
        3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
        3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
        3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 0.0, 0.0,
        T*3, T*3, T*3, T*3, T*3, T*3, T*3, T*3,
        T*3, T*3, T*3, T*3,
        (R/(R+1))*3, (R/(R+1))*3, 0.0, 0.0,
        3.0, 3.0, 1.5*3, (1+1/(R+2))*3,
        T*3, T*3, (T+1/(R+1))*3, (T+1/(R+2))*3,
        3.0, 3.0, 3.0, 3.0
    ]
    return sites[cn]


def is_transition(b0, b1):
    """True if the two base indices are a transition (A<->G or C<->T)."""
    return (b0 == 0 and b1 == 1) or (b0 == 1 and b1 == 0) or \
           (b0 == 2 and b1 == 3) or (b0 == 3 and b1 == 2)


def compute_cs(u0, v0, w0, codon_numberC, conserv_aa, R):
    """
    Compute the conservative site contribution for one codon (one of the two
    sequences).  Mirrors the Perl inner loop in countsubstitutions.
    """
    cs = 0.0

    # --- first codon position (u) ---
    flag = flag1 = 0
    u = u0
    for _ in range(3):
        u = (u + 1) % 4
        cn = (v0 * 16) + (u * 4) + w0
        if AA_NUMBER[cn] == 21:
            flag += 1
            if is_transition(u0, u):
                flag1 = 1

    u = u0
    for _ in range(3):
        u = (u + 1) % 4
        cn = (v0 * 16) + (u * 4) + w0
        if (conserv_aa[AA_NUMBER[cn]] == conserv_aa[AA_NUMBER[codon_numberC]]
                and AA_NUMBER[cn] != AA_NUMBER[codon_numberC]):
            if flag == 0:
                cs += R / (R + 2) if is_transition(u0, u) else 1 / (R + 2)
            elif flag == 2:
                cs += 1.0
            elif flag == 1 and flag1 == 1:
                cs += 0.5
            elif flag == 1 and flag1 == 0:
                cs += R / (R + 1) if is_transition(u0, u) else 1 / (R + 1)

    # --- second codon position (v) ---
    flag = flag1 = 0
    v = v0
    for _ in range(3):
        v = (v + 1) % 4
        cn = (v * 16) + (u0 * 4) + w0
        if AA_NUMBER[cn] == 21:
            flag += 1
            if is_transition(v0, v):
                flag1 = 1

    v = v0
    for _ in range(3):
        v = (v + 1) % 4
        cn = (v * 16) + (u0 * 4) + w0
        if (conserv_aa[AA_NUMBER[cn]] == conserv_aa[AA_NUMBER[codon_numberC]]
                and AA_NUMBER[cn] != AA_NUMBER[codon_numberC]):
            if flag == 0:
                cs += R / (R + 2) if is_transition(v0, v) else 1 / (R + 2)
            elif flag == 2:
                cs += 1.0
            elif flag == 1 and flag1 == 1:
                cs += 0.5
            elif flag == 1 and flag1 == 0:
                cs += R / (R + 1) if is_transition(v0, v) else 1 / (R + 1)

    # --- third codon position (w) ---
    flag = flag1 = 0
    w = w0
    for _ in range(3):
        w = (w + 1) % 4
        cn = (v0 * 16) + (u0 * 4) + w
        if AA_NUMBER[cn] == 21:
            flag += 1
            if is_transition(w0, w):
                flag1 = 1

    w = w0
    for _ in range(3):
        w = (w + 1) % 4
        cn = (v0 * 16) + (u0 * 4) + w
        if (conserv_aa[AA_NUMBER[cn]] == conserv_aa[AA_NUMBER[codon_numberC]]
                and AA_NUMBER[cn] != AA_NUMBER[codon_numberC]):
            if flag == 0:
                cs += R / (R + 2) if is_transition(w0, w) else 1 / (R + 2)
            elif flag == 2:
                cs += 1.0
            elif flag == 1 and flag1 == 1:
                cs += 0.5
            elif flag == 1 and flag1 == 0:
                cs += R / (R + 1) if is_transition(w0, w) else 1 / (R + 1)

    return cs


def tally_2base(codon_numberA, codon_numberB, codC, codD, conserv_aa):
    """
    Shared logic for 2-base-change cases.
    Returns (tmp_syn, tmp_conserv, temp) without dividing by (2-temp).
    """
    codon_numberC = codon_conversion(codC)
    codon_numberD = codon_conversion(codD)
    tmp_syn = temp = tmp_conserv = 0.0

    for cn_mid in (codon_numberC, codon_numberD):
        if AA_ARRAY[cn_mid] == 'Z':
            temp += 1
        else:
            for cn_end in (codon_numberA, codon_numberB):
                if AA_ARRAY[cn_end] == AA_ARRAY[cn_mid]:
                    tmp_syn += 1
                elif conserv_aa[AA_NUMBER[cn_end]] == conserv_aa[AA_NUMBER[cn_mid]]:
                    tmp_conserv += 1

    return tmp_syn, tmp_conserv, temp


STOP_CODONS = {"TAA", "TAG", "TGA"}


def validate_sequences(names, seqs):
    """
    Enforce biological and mathematical requirements on input sequences.
    Checks performed:
      1. At least two sequences present.
      2. All sequences are the same length (equal-length MSA after complete
         gap deletion required).
      3. Sequence length is divisible by 3 (complete codons only).
      4. Only A, C, G, T characters are present.
      5. No internal stop codons (TAA, TAG, TGA).
      6. Terminal stop codon stripped from all sequences if present.
         (A mixed terminal stop situation cannot occur given checks 2 and 3.)
    """
    if len(seqs) < 2:
        _error("At least two sequences are required.")

    lengths = [len(s) for s in seqs]
    if len(set(lengths)) != 1:
        detail = ", ".join(f"{n}: {l}" for n, l in zip(names, lengths))
        _error(
            "All sequences must be the same length (complete gap deletion "
            f"required).\n  Lengths found — {detail}"
        )
    seq_len = lengths[0]

    if seq_len % 3 != 0:
        _error(
            f"Sequence length ({seq_len}) is not divisible by 3. "
            "Only complete codons are supported."
        )

    valid = set("ACGT")
    for name, seq in zip(names, seqs):
        bad = sorted(set(seq) - valid)
        if bad:
            _error(
                f"Sequence '{name}' contains invalid character(s): "
                f"{bad}. Only A, C, G, T are permitted."
            )

    for name, seq in zip(names, seqs):
        for pos in range(0, seq_len - 3, 3):
            codon = seq[pos:pos + 3]
            if codon in STOP_CODONS:
                codon_num = pos // 3 + 1
                _error(
                    f"Sequence '{name}' contains an internal stop codon "
                    f"({codon}) at codon position {codon_num}. "
                    "Internal stop codons are not permitted for this method."
                )

    terminal_stops = [seq[-3:] in STOP_CODONS for seq in seqs]
    if all(terminal_stops):
        print(
            f"[INFO] Terminal stop codon detected in all {len(seqs)} sequences "
            "— stripped from all sequences before calculation.",
            file=sys.stderr
        )
        seqs = [s[:-3] for s in seqs]

    return seqs


def _error(msg):
    print(f"[ERROR] {msg}", file=sys.stderr)
    sys.exit(1)


def main():
    if len(sys.argv) < 3:
        print("Usage: python pnc_pnr_dist.py <input_fasta> <property_file> [R or ratio]")
        sys.exit(1)

    input_file = sys.argv[1]
    prop_file  = sys.argv[2]
    ratio_arg  = sys.argv[3] if len(sys.argv) > 3 else "0.5"
    base_ratio = 'R' if ratio_arg.upper() == 'R' else float(ratio_arg)

    # --- read property file ---
    # Tab-delimited: header line, then one row per amino acid (20 rows).
    # Column 0 = amino acid name, column 1 = property value.
    # conserv_aa is 1-indexed to match AA_NUMBER (index 21 = stop = 0).
    conserv_aa = [0.0] * 22
    header1 = ""
    with open(prop_file, 'r') as fh:
        header1 = fh.readline().rstrip('\n')
        idx = 0
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            idx += 1
            conserv_aa[idx] = float(parts[1])
    conserv_aa[21] = 0.0

    # --- parse FASTA input ---
    names = []
    seqs  = []
    current_seq = []
    with open(input_file, 'r') as fh:
        for line in fh:
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

    # --- validate all sequences before any calculation ---
    seqs = validate_sequences(names, seqs)

    seq_length = len(seqs[0])
    num_seqs   = len(names)

    # --- estimate global R from all sequences if requested ---
    # R is estimated once from all 3rd codon positions across all pairs,
    # consistent with the stationarity assumption of the K2P / Nei-Gojobori
    # method. A single global estimate is more reliable than per-pair
    # estimation and avoids undefined results from low-divergence pairs.
    if base_ratio == 'R':
        P1 = P2 = Q = 0
        num_codons = seq_length // 3
        for na in range(num_seqs):
            for n in range(na + 1, num_seqs):
                for x in range(0, seq_length, 3):
                    f = seqs[na][x + 2]
                    s = seqs[n][x + 2]
                    if f != s:
                        if   f in 'AG' and s in 'AG': P1 += 1
                        elif f in 'AG' and s in 'CT': Q  += 1
                        elif f in 'CT' and s in 'AG': Q  += 1
                        elif f in 'CT' and s in 'CT': P2 += 1
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

    # --- output header ---
    print(f"gene1\tgene2\tratio\t{header1}\tSE_{header1}\tpnr\tSEpnr")

    R = base_ratio * 2.0
    T = R / (R + 2.0)

    # --- pairwise loop ---
    for na in range(num_seqs):
        seqA = list(seqs[na])
        for n in range(na + 1, num_seqs):
            seqB = list(seqs[n])

            syn_codons = nonsyn_codons = conserv_codons = 0.0
            SA_Nei = SB_Nei = conserv_site = 0.0
            count_codons = 0

            # Perl loop: i=1 to seq_length-1 step 3 (1-based), accessing i-1,i,i+1
            # Python 0-based equivalent: range(0, seq_length-2, 3)
            for i in range(0, seq_length - 2, 3):
                codA = [seqA[i], seqA[i+1], seqA[i+2]]
                codB = [seqB[i], seqB[i+1], seqB[i+2]]
                count_codons += 1

                codon_numberA = codon_conversion(codA)
                codon_numberB = codon_conversion(codB)

                SA_Nei += syn_site(codA, T, R)
                SB_Nei += syn_site(codB, T, R)

                # conservative-site contribution for each of the two codons
                conserv_site += compute_cs(
                    BASE_NUMBER[seqA[i]], BASE_NUMBER[seqA[i+1]], BASE_NUMBER[seqA[i+2]],
                    codon_numberA, conserv_aa, R)
                conserv_site += compute_cs(
                    BASE_NUMBER[seqB[i]], BASE_NUMBER[seqB[i+1]], BASE_NUMBER[seqB[i+2]],
                    codon_numberB, conserv_aa, R)

                # --- no change ---
                if codA == codB:
                    continue

                diffs = sum(1 for j in range(3) if codA[j] != codB[j])

                # --- 1-base change ---
                if diffs == 1:
                    if AA_ARRAY[codon_numberA] == AA_ARRAY[codon_numberB]:
                        syn_codons += 1.0
                    else:
                        if conserv_aa[AA_NUMBER[codon_numberA]] == conserv_aa[AA_NUMBER[codon_numberB]]:
                            conserv_codons += 1.0
                        nonsyn_codons += 1.0

                # --- 2-base changes ---
                elif diffs == 2:
                    if   codA[0] != codB[0] and codA[1] != codB[1]:   # pos 0 & 1
                        codC = [codB[0], codA[1], codA[2]]
                        codD = [codA[0], codB[1], codA[2]]
                    elif codA[0] != codB[0] and codA[2] != codB[2]:   # pos 0 & 2
                        codC = [codB[0], codA[1], codA[2]]
                        codD = [codA[0], codA[1], codB[2]]
                    else:                                               # pos 1 & 2
                        codC = [codA[0], codB[1], codA[2]]
                        codD = [codA[0], codA[1], codB[2]]

                    tmp_syn, tmp_conserv, temp = tally_2base(
                        codon_numberA, codon_numberB, codC, codD, conserv_aa)
                    denom = 2.0 - temp
                    tmp_syn /= denom
                    syn_codons    += tmp_syn
                    nonsyn_codons += 2.0 - tmp_syn
                    conserv_codons += tmp_conserv / denom

                # --- 3-base changes ---
                elif diffs == 3:
                    codC = [codB[0], codA[1], codA[2]];  codon_numberC = codon_conversion(codC)
                    codF = [codB[0], codB[1], codA[2]];  codon_numberF = codon_conversion(codF)
                    codG = [codB[0], codA[1], codB[2]];  codon_numberG = codon_conversion(codG)
                    codE = [codA[0], codA[1], codB[2]];  codon_numberE = codon_conversion(codE)
                    codH = [codA[0], codB[1], codB[2]];  codon_numberH = codon_conversion(codH)
                    codD = [codA[0], codB[1], codA[2]];  codon_numberD = codon_conversion(codD)

                    temp = tmp_syn = tmp_conserv = 0.0

                    # six A->B paths through two intermediates
                    paths = [
                        (codon_numberC, codon_numberF),
                        (codon_numberC, codon_numberG),
                        (codon_numberD, codon_numberF),
                        (codon_numberD, codon_numberH),
                        (codon_numberE, codon_numberG),
                        (codon_numberE, codon_numberH),
                    ]
                    for cn_mid1, cn_mid2 in paths:
                        if AA_ARRAY[cn_mid1] == 'Z' or AA_ARRAY[cn_mid2] == 'Z':
                            temp += 1
                        else:
                            for pair in ((codon_numberA, cn_mid1),
                                         (cn_mid1,       cn_mid2),
                                         (codon_numberB, cn_mid2)):
                                cn_x, cn_y = pair
                                if AA_ARRAY[cn_x] == AA_ARRAY[cn_y]:
                                    tmp_syn += 1
                                elif conserv_aa[AA_NUMBER[cn_x]] == conserv_aa[AA_NUMBER[cn_y]]:
                                    tmp_conserv += 1

                    denom = 6.0 - temp
                    conserv_codons += tmp_conserv / denom
                    syn_codons     += tmp_syn / denom
                    nonsyn_codons  += 3.0 - (tmp_syn / denom)

                else:
                    print(f"Error\t{''.join(codA)}\t{AA_ARRAY[codon_numberA]}"
                          f"\t{''.join(codB)}\t{AA_ARRAY[codon_numberB]}")

            # --- summary statistics ---
            potential_syn    = (SA_Nei / 3.0 + SB_Nei / 3.0) / 2.0
            potential_nonsyn = 3.0 * count_codons - potential_syn
            potential_conserv = conserv_site / 2.0
            potential_rad    = potential_nonsyn - potential_conserv
            rad_codons       = nonsyn_codons - conserv_codons

            ps  = syn_codons    / potential_syn
            pn  = nonsyn_codons / potential_nonsyn
            pnc = conserv_codons / potential_conserv
            pnr = rad_codons    / potential_rad

            SEpnc = math.sqrt(pnc * (1.0 - pnc) / potential_conserv)
            SEpnr = math.sqrt(pnr * (1.0 - pnr) / potential_rad)

            print(f"{names[na]}\t{names[n]}\t{base_ratio:.2f}\t{pnc:.4f}\t{SEpnc:.4f}\t{pnr:.4f}\t{SEpnr:.4f}")


if __name__ == '__main__':
    main()
