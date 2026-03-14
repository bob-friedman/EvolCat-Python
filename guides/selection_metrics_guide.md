# Guide: Selection Metrics in Evolutionary Genetics

This guide provides an overview of the methods used to detect natural selection in protein-coding sequences, focusing on the dN/dS and pNC/pNR metrics implemented in this repository.

---

## 1. Introduction to Selection Metrics

Natural selection leaves a distinct signature on the nucleotide sequences of protein-coding genes. By comparing the rates of different types of substitutions, we can infer the evolutionary pressures acting on a gene or specific regions within it.

- **Purifying (Negative) Selection:** Acts to eliminate deleterious mutations, preserving protein function.
- **Neutral Evolution:** Mutations have no significant effect on fitness and accumulate at a rate determined by mutation and genetic drift.
- **Positive (Darwinian) Selection:** Acts to favor beneficial mutations that increase fitness, promoting adaptation and diversification.

---

## 2. The dN/dS Method (Nei & Gojobori, 1986)

The dN/dS ratio (also denoted as ω or omega) is the classic measure for detecting selection.

### The Concept
- **dS (Synonymous substitutions per synonymous site):** Mutations that do not change the amino acid. These are generally assumed to be neutral.
- **dN (Nonsynonymous substitutions per nonsynonymous site):** Mutations that change the amino acid.

### Interpretation
- **dN/dS < 1:** Purifying selection. Nonsynonymous changes are being suppressed.
- **dN/dS ≈ 1:** Neutral evolution.
- **dN/dS > 1:** Positive selection. Amino acid changes are being favored.

### Assumptions and Limitations
- **Synonymous Neutrality:** Assumes that synonymous mutations have no effect on fitness (ignoring codon usage bias or mRNA stability).
- **Time Scale:** Best suited for comparing sequences that have diverged enough to accumulate multiple substitutions but not so much that they are saturated.
- **Average Effect:** A dN/dS ratio calculated over an entire gene may mask positive selection acting only on a few specific sites.

### Repository Tools
- **`pylib/scripts/standalone/metrics/dsdn_dist.py`**: A standalone tool for pairwise dN/dS estimation using the Nei & Gojobori (1986) method.
- **`pylib/scripts/paml_tools/calculate_dn_ds.py`**: A wrapper for PAML's `yn00` program.
- **`pylib/scripts/paml_tools/calculate_site_specific_ds_dn.py`**: A wrapper for PAML's `codeml` for detecting selection at individual codon sites.

---

## 3. The pNC/pNR Method (Hughes et al., 1990)

Developed by Hughes, Ota, and Nei, this method provides a more nuanced look at nonsynonymous substitutions by categorizing them based on amino acid properties.

### The Concept
Nonsynonymous substitutions are further divided into:
- **pNC (Proportion of conservative substitutions):** Changes between amino acids with similar properties (e.g., both neutral).
- **pNR (Proportion of radical substitutions):** Changes between amino acids with different properties (e.g., neutral to charged).

### Application to MHC Molecules
The original study (Hughes et al., 1990) applied this to the Major Histocompatibility Complex (MHC). They found that in the antigen-binding cleft, radical substitutions (specifically those changing charge) occurred significantly more frequently than expected by chance (**pNR > pNC**). This provided strong evidence for positive selection favoring diversity in the charge profile of the binding cleft to recognize a wider range of foreign peptides.

### Interpretation
- **pNR > pNC:** Evidence that selection is actively promoting radical changes in amino acid properties (positive selection for diversity).
- **pNC > pNR:** Evidence of conservative selection, where the protein's biochemical properties are being maintained even when the amino acid changes.

### Assumptions and Limitations
- **Property Definition:** The results depend entirely on how "conservative" and "radical" are defined (e.g., by charge, polarity, or volume).
- **Baseline Expectation:** Assumes that under random nonsynonymous substitution, pNC and pNR should be equal.

### Repository Tools
- **`pylib/scripts/standalone/metrics/pnc_pnr_dist.py`**: Estimates pNC/pNR metrics based on a user-provided property file.
- **Sample Data**: `pylib/scripts/standalone/test_data/property_charge` provides the classification used in the Hughes et al. (1990) paper.

---

## 4. Practical Workflow

1.  **Alignment:** Obtain a high-quality codon-aware alignment (e.g., using `pal2nal_enforce.py`).
2.  **Gap Removal:** Remove all columns containing gaps (using `nogaps.py`) as these methods require complete codons.
3.  **Basic Selection Test:** Run `dsdn_dist.py` to get an overview of the dN/dS ratio.
4.  **Property-Specific Analysis:** If positive selection is suspected, run `pnc_pnr_dist.py` with relevant properties (like charge or hydrophobicity) to understand the nature of the selective pressure.
5.  **Site-Specific Analysis:** Use `calculate_site_specific_ds_dn.py` (via PAML) to identify the specific residues under selection.

---

## 5. References

- **Nei, M., & Gojobori, T. (1986).** Simple methods for estimating the numbers of synonymous and nonsynonymous nucleotide substitutions. *Molecular Biology and Evolution*, 3(4), 418-426.
- **Hughes, A. L., Ota, T., & Nei, M. (1990).** Positive Darwinian selection promotes charge profile diversity in the antigen-binding cleft of class I major-histocompatibility-complex molecules. *Molecular Biology and Evolution*, 7(6), 515-524.
- **Hughes, A. L., & Nei, M. (1988).** Pattern of nucleotide substitution at major histocompatibility complex loci reveals overdominant selection. *Nature*, 335(6186), 167-170.
