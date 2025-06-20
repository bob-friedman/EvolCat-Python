# Guide: Estimating Mutational Fitness Effects from Viral Sequences

> This guide outlines a conceptual methodology for estimating the fitness effects of viral mutations using large-scale public sequence data. It is primarily based on the powerful approach described by **Bloom and Neher (2023)** in their study of SARS-CoV-2.

---

## Table of Contents
1.  [Introduction](#introduction)
2.  [The Foundational Data Source](#the-foundational-data-source-mutation-annotated-phylogenetic-tree)
3.  [The 6-Step Methodology](#the-6-step-methodology)
4.  [Key Considerations and Tools](#key-considerations-and-tools)
5.  [References](#references)

---

## Introduction

Understanding the fitness impact of mutations is crucial for predicting viral evolution and designing effective vaccines. While experimental methods like Deep Mutational Scanning (DMS) are powerful, computational approaches analyzing vast public sequence datasets offer a complementary way to map mutational effects across entire viral proteomes.

> **The Core Idea:** Compare the *observed* frequency of each mutation in a massive collection of viral sequences against its *expected* frequency under a model of neutral evolution. Significant deviations reveal selection pressures.

---

## The Foundational Data Source: Mutation-Annotated Phylogenetic Tree

This methodology hinges on a large, high-quality, mutation-annotated phylogenetic tree. The **UShER toolkit** and its pre-built trees for SARS-CoV-2 are the gold standard for this.

*   **Example Data from Bloom & Neher (2023):**
    *   **UShER Tree File:** `public-2023-05-11.all.masked.nextclade.pangolin.pb.gz`
    *   **Source:** Daily-updated databases are available from the [UShER SARS-CoV-2 repository](http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/).
*   **Essential Toolkit:**
    *   The [UShER (Ultrafast Sample placement on Existing tRee)](https://usher-wiki.readthedocs.io/en/latest/) toolkit is required for processing these `.pb` (protobuf) tree files.
    *   **`matUtils`**, a key part of the toolkit, is used to extract mutation, sample, and phylogenetic data from the tree. See our [condensed guide to matUtils](./usher_toolkit_report.md).

---

## The 6-Step Methodology

Here is a conceptual workflow for estimating mutational fitness effects.

### **Step 1: Obtain and Prepare the Tree**
*   Download a pre-built UShER `.pb` file.
*   Use `matUtils` to parse the tree and extract all nucleotide mutations annotated on its branches. This is often done on a per-clade basis to account for different genetic backgrounds.

### **Step 2: Apply Rigorous Quality Control (QC)**
This is a critical step to filter out noise and potential artifacts.
*   **Filter Branches:** Exclude branches with an excessive number of mutations, multiple reversions, or multiple mutations in the same codon. These often indicate sequencing errors.
*   **Filter Sites/Mutations:** Exclude known problematic sites, error-prone mutations, and suspected reversion artifacts.

### **Step 3: Calculate Expected Neutral Mutation Counts**
The goal is to establish a baseline for how often mutations *should* occur without selection.
*   **Identify Neutral Sites:** Use **four-fold degenerate synonymous sites** as a proxy for neutral evolution. These are codon positions where any nucleotide change results in the same amino acid.
*   **Count Neutral Mutations:** Within each major clade, count the occurrences of each mutation type (e.g., C→T, A→G) *only at these neutral sites*.
*   **Calculate Expected Rate:** For each mutation type (e.g., X→Y), the expected rate is the total number of observed X→Y mutations at neutral sites, divided by the total number of neutral sites that could have made that mutation.

### **Step 4: Count Actual Observed Mutations**
*   For *every site* in the genome (not just neutral ones), count the actual number of times each specific nucleotide mutation is observed across all branches (after QC).

### **Step 5: Aggregate and Convert to Amino Acid Effects**
*   For each amino acid change (e.g., Alanine to Glycine at site 100), sum the *actual* and *expected* counts of all underlying nucleotide mutations that could cause it.
*   For pan-virus estimates, these counts can be summed across all major clades.

### **Step 6: Estimate Fitness Effect (Δf)**
The fitness effect of a mutation is the natural log of the ratio of its actual to expected frequency. A pseudocount (P) is added to prevent division by zero and stabilize estimates.

> **Formula:**
> `Δf = ln( (n_actual + P) / (n_expected + P) )`

*   **`n_actual`**: Total observed counts for the mutation.
*   **`n_expected`**: Total expected counts for the mutation.
*   **`P`**: A small pseudocount (e.g., 0.5).

**Interpretation of `Δf`:**
*   **`Δf ≈ 0`**: Suggests the mutation is **neutral**.
*   **`Δf < 0`**: Suggests the mutation is **deleterious** (negative selection).
*   **`Δf > 0`**: Suggests the mutation may be **beneficial** (positive selection).

---

## Key Considerations and Tools

*   **Computational Scale:** This analysis requires significant computational resources (CPU, memory, storage) and is not a desktop task.
*   **Essential Software:**
    *   **UShER Toolkit (`matUtils`)**: Non-negotiable for handling `.pb` files.
    *   **Custom Scripts (Python/R)**: Required for QC, counting, and calculations. Libraries like `Biopython` are useful for sequence and annotation handling.
*   **Complexity:** This is a sophisticated bioinformatic approach. The choice of QC filters, clade definitions, and pseudocounts can significantly impact the results. A thorough understanding of the underlying assumptions is crucial.
*   **Relation to EvolCat-Python:** While this pipeline is not implemented in EvolCat-Python, the principles of analyzing selection are related to those in `calculate_site_specific_ds_dn.py`, and other scripts could be adapted for pre-processing or downstream analysis.

For a full implementation, please refer to the methods and publicly available code from the original Bloom and Neher (2023) paper.

---

## References

*   Bloom, J. D., & Neher, R. A. (2023). **Fitness effects of mutations to SARS-CoV-2 proteins**. *Virus Evolution*, *9*(2), vead055. [https://doi.org/10.1093/ve/vead055](https://doi.org/10.1093/ve/vead055)
*   McBroome, J., et al. (2021). **A Daily-Updated Database and Tools for Comprehensive SARS-CoV-2 Mutation-Annotated Trees**. *Molecular Biology and Evolution*, *38*(12), 5819–5824. [https://doi.org/10.1093/molbev/msab226](https://doi.org/10.1093/molbev/msab226)
