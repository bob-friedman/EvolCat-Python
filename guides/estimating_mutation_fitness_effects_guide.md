<a name="top"></a>
# Guide to Estimating Mutational Fitness Effects from Large-Scale Viral Sequence Data

This guide outlines a conceptual methodology for estimating the fitness effects of mutations in viruses using large-scale sequence data, largely based on the approach described by Bloom and Neher (2023) in their study of SARS-CoV-2. This method leverages public genomic data to create comprehensive maps of mutational effects, which can be invaluable for understanding viral evolution, protein function, and informing public health strategies.

## Introduction

Understanding the impact of mutations on viral fitness is crucial for predicting viral evolution, assessing the risk of new variants, and designing effective antiviral therapies and vaccines. While experimental methods like Deep Mutational Scanning (DMS) provide valuable data, computational approaches analyzing vast numbers of publicly available viral sequences offer a complementary way to estimate fitness effects across entire viral proteomes.

The core idea is to compare the *observed* frequency of each mutation in a large collection of viral sequences against its *expected* frequency under a model of neutral evolution. Deviations from this expectation can indicate selection pressures acting on the mutation.

## Data Source: Mutation-Annotated Phylogenetic Tree

This methodology relies heavily on a large, comprehensive, and accurately constructed mutation-annotated phylogenetic tree. For their SARS-CoV-2 analysis, Bloom and Neher (2023) utilized a publicly available, pre-built UShER tree.

*   **Specific Data Used in Bloom and Neher (2023):**
    *   **UShER Tree File:** `http://vhgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/2023/05/11/public-2023-05-11.all.masked.nextclade.pangolin.pb.gz`
        *   More recent versions of the dataset (McBroome et al. 2021) are available: `http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/`
    *   **Contents:** This compressed file contains a mutation-annotated tree representing approximately 7 million public SARS-CoV-2 sequences available as of May 11, 2023. The tree includes annotations of nucleotide mutations on its branches.
*   **UShER Toolkit:** The [UShER (Ultrafast Sample placement on Existing tRee)](https://usher-wiki.readthedocs.io/en/latest/) toolkit, particularly programs like `matUtils`, is essential for working with these `.pb` (protobuf) tree files. [matUtils](https://usher-wiki.readthedocs.io/en/latest/matUtils.html) can be used to extract information about mutations, samples, and phylogenetic relationships from the tree.

Obtaining and being able to parse such a tree is the foundational step for this type of analysis.

## Conceptual Methodology Outline

The following steps describe the conceptual workflow for estimating mutational fitness effects:

### 1. Obtaining and Understanding the Mutation-Annotated Tree
*   Download the pre-built UShER tree file (e.g., the link provided above for SARS-CoV-2).
*   Familiarize yourself with the UShER toolkit and `matUtils` for tree manipulation and data extraction. This tool allows you to, for instance, extract mutations annotated on each branch of the tree.

### 2. Extracting Mutation Data
*   Use `matUtils` to extract nucleotide mutations from all branches of the tree. The analysis is often done on a per-clade basis to account for different genetic backgrounds and potentially varying mutation spectra.
    *   The paper by Bloom and Neher (2023) subsetted the global tree into major Nextstrain clades.

### 3. Quality Control (QC)
This is a critical step to ensure the reliability of the results. The original paper implemented several QC measures:
*   **Branch-level QC:**
    *   Ignoring branches with an excessive number of nucleotide mutations (e.g., >4 in the paper), which might indicate sequencing errors or atypical evolution (like in chronic infections).
    *   Filtering branches with multiple reversions to the reference or clade founder sequence, which can be artifacts.
    *   Ignoring branches with multiple mutations in the same codon (rare, but simplifies analysis).
*   **Site/Mutation-level QC:**
    *   Excluding known problematic sites prone to sequencing or base-calling errors (often listed in public resources or discovered through analysis).
    *   Excluding specific mutations identified as highly error-prone in certain clades.
    *   Masking reversions from clade founders back to the Wuhan-Hu-1 reference (and their complements) if these are suspected artifacts.

### 4. Calculating Expected Mutation Counts
The goal is to determine how often each type of nucleotide mutation (e.g., A_to_G, C_to_T) should occur neutrally.
*   **Identify Neutral Sites:** Four-fold degenerate synonymous sites are commonly used. These are codon positions where any nucleotide change does not alter the encoded amino acid.
    *   This requires accurate genome annotation to identify coding regions and codon boundaries.
*   **Per-Clade Calculation:** Expected counts are calculated separately for each major viral clade.
*   **Counting Mutations at Neutral Sites:** Count the occurrences of each specific nucleotide mutation type (e.g., C->T, A->G) only at these identified four-fold degenerate sites within each clade.
*   **Calculate Expected Rate:** For each mutation type (e.g., from nucleotide `X` to `Y`), the expected count per site is the total number of `X` to `Y` mutations observed at four-fold degenerate sites, divided by the total number of four-fold degenerate sites that had nucleotide `X` as the parental identity in that clade.

### 5. Counting Actual Observed Mutations
*   For every site in the genome (not just neutral sites), count the actual number of times each specific nucleotide mutation is observed on the branches of the phylogenetic tree (after QC). This is done across all sequences/branches within each clade.

### 6. Estimating Fitness Effects
*   **Convert Nucleotide to Amino Acid Effects:** For coding regions, sum the counts of all non-excluded nucleotide mutations that result in a specific amino acid change (e.g., Alanine to Glycine at site 100).
    *   The analysis typically excludes mutations that are not single-step changes from the clade-founder codon identity at that position.
*   **Aggregate Across Clades (Optional):** For pan-virus estimates, these counts (both actual and expected for the corresponding nucleotide changes) can be summed across all analyzed clades.
*   **Calculate Fitness Effect (Δf):** The fitness effect of each mutation is estimated using the natural logarithm of the ratio of actual to expected counts. A pseudocount (P) is often added to both numerator and denominator to stabilize estimates when counts are low.
    The formula used is:
    `Δf = ln( (n_actual + P) / (n_expected + P) )`
    *   `n_actual`: Sum of actual observed counts for nucleotide mutations encoding the amino acid change.
    *   `n_expected`: Sum of expected counts for those same nucleotide mutations.
    *   `P`: A small pseudocount (e.g., 0.5).
*   **Interpretation:**
    *   `Δf ≈ 0`: Suggests near-neutrality.
    *   `Δf < 0`: Suggests the mutation is deleterious.
    *   `Δf > 0`: May suggest the mutation is beneficial, though these are typically rare and require careful interpretation.

## Considerations and Tools

*   **Computational Resources:** Processing massive phylogenetic trees and sequence datasets requires significant computational resources and specialized tools.
*   **UShER Toolkit:** Essential for handling `.pb` tree files and extracting mutation data.
*   **Bioinformatic Scripting:** Custom scripts (e.g., in Python with libraries like Biopython for sequence manipulation and parsing) are necessary for:
    *   Identifying four-fold degenerate sites from genome annotations.
    *   Implementing the QC logic.
    *   Performing the counting for expected and actual mutations.
    *   Calculating the final fitness ratios.
*   **EvolCat-Python conceptual links:**
    *   While EvolCat-Python does not currently implement this full pipeline, scripts for FASTA manipulation (`clean_fasta_name.py`, `merge_fastas.py`), extracting regions (`extract_region.py`), or parsing sequence data could be adapted for pre-processing or downstream analysis of specific gene regions once mutation data is extracted.
    *   The principles of analyzing selection are related to those in `calculate_site_specific_ds_dn.py`, although the underlying methodology (counts vs. codon models) is different.
*   **Complexity:** This is a sophisticated analytical approach. The details of QC, clade definition, and pseudocount selection can influence the results. Understanding the assumptions and limitations (see "Important Caveats" in the summary within the [Virus Genomics Guide](./virus_genomics_guide.md)) is crucial.

This guide provides a high-level overview. For a detailed implementation, refer to the methods in the original Bloom and Neher (2023) paper and their publicly available code. There is also a related guide on use of this toolkit as it applies to the above documentation: [UShER Toolkit and matUtils: A Condensed Guide](./usher_toolkit_report.md).

## References
*   Bloom, J. D., & Neher, R. A. (2023). Fitness effects of mutations to SARS-CoV-2 proteins. *Virus Evolution*, *9*(2), vead055. [https://doi.org/10.1093/ve/vead055](https://doi.org/10.1093/ve/vead055)
*   McBroome, J., Thornlow, B., Hinrichs, A.S., Kramer, A., De Maio, N., Goldman, N., & Turakhia , Y. (2021).  A Daily-Updated Database and Tools for Comprehensive SARS-CoV-2 Mutation-Annotated Trees. [https://academic.oup.com/mbe/article/38/12/5819/6361626](https://academic.oup.com/mbe/article/38/12/5819/6361626)

[Back to Top](#top)
