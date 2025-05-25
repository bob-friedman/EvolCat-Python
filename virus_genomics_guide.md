# Guide to Virus Genomics, Diversity, and Analysis

## Introduction to Virus Genomics
(Placeholder for content: Briefly explain what virus genomics is, its importance, and key characteristics of viral genomes, e.g., size, mutation rates, RNA vs DNA.)

## Measuring Viral Diversity
(Placeholder for content: Discuss common metrics and methods used to quantify viral diversity, e.g., nucleotide diversity, haplotype analysis, quasispecies.)

## Phylogenetic Analysis of Viruses
(Placeholder for content: Explain how phylogenetic methods are applied to viruses, what insights can be gained, e.g., tracking outbreaks, understanding evolutionary relationships. Mention specific EvolCat-Python scripts or external tools if applicable.)

## Tools and Resources
(Placeholder for content: List relevant EvolCat-Python scripts that can be used for viral genomics. Also, list key external databases, software, and web resources for viral genomics research.)

### Site-Specific dN/dS Analysis (`calculate_site_specific_ds_dn.py`)

This script serves as a Python wrapper for the `codeml` program from the PAML (Phylogenetic Analysis by Maximum Likelihood) package. Its primary purpose is to automate and simplify the process of identifying natural selection acting at individual codon sites within a set of aligned coding sequences.

The dN/dS ratio (often denoted as ω or omega) is a key metric in these analyses:
*   **dN/dS > 1**: Suggests positive (Darwinian) selection, where mutations leading to amino acid changes are favored.
*   **dN/dS < 1**: Indicates purifying (negative) selection, where mutations changing amino acids are deleterious and removed.
*   **dN/dS = 1**: Points towards neutral evolution, where amino acid changes have no selective advantage or disadvantage.

**Main Inputs:**
*   **Coding Sequence Alignment (FASTA format):** An alignment of homologous coding sequences (e.g., a specific gene from different viral strains). Sequences must be a multiple of three nucleotides.
*   **Phylogenetic Tree (Newick format):** A tree representing the evolutionary relationships among the sequences in the alignment. Branch lengths are important.
*   **PAML Model:** The specific site model to be used by `codeml` (e.g., M0, M1a, M2a, M3, M7, M8). Different models allow for different assumptions about the distribution of dN/dS ratios across sites.

**Primary Output:**
The script generates several output files prefixed by a user-defined string. The most important for site-specific analysis is a tab-separated values (TSV) file named `<outfile_prefix>_site_analysis.tsv`. This file contains:
*   Site number (codon position in the alignment).
*   The amino acid present at that site in a reference sequence (if applicable from PAML output).
*   The estimated dN/dS ratio for that site (or the dN/dS class the site is assigned to, depending on the model).
*   Posterior probabilities from Bayes Empirical Bayes (BEB) analysis, which are crucial for identifying sites under positive selection with statistical confidence (especially relevant for models like M2a and M8).
*   A "Note" column highlighting sites with significant evidence of positive selection (e.g., BEB probability > 0.95).

**PAML Dependency:**
**Crucially, this script requires PAML (specifically the `codeml` executable) to be installed on your system.** `codeml` must either be accessible in your system's PATH environment variable, or its direct path must be provided to the script using the `--paml_path` argument.

The `calculate_site_specific_ds_dn.py` script automates the creation of `codeml` control files, executes `codeml`, and parses its often complex output files into a more structured and user-friendly TSV format, facilitating the identification of codons potentially under selection.

## Example Workflow: Analyzing a Viral Dataset
(Placeholder for content: Outline a hypothetical step-by-step workflow for analyzing a small viral dataset, from sequence retrieval to diversity calculation and phylogenetic tree construction, mentioning which tools (EvolCat-Python or external) would be used at each step.)
