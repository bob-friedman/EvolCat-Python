# A Guide to the UShER Toolkit and `matUtils`

This guide provides a condensed overview of the UShER toolkit, focusing on its core command-line utility, `matUtils`, for analyzing and manipulating large-scale phylogenetic trees.

> **What is UShER?**
> UShER (Ultrafast Sample placement on Existing tRee) is a high-performance toolkit essential for real-time genomic epidemiology. It excels at placing new pathogen sequences onto massive existing phylogenetic trees, making it invaluable for tracking viral evolution during outbreaks like the COVID-19 pandemic.

---

## Table of Contents
1.  [The Mutation-Annotated Tree (MAT)](#1-the-mutation-annotated-tree-mat)
2.  [`matUtils`: The MAT Power Tool](#2-matutils-the-mat-power-tool)
3.  [Getting Started: The SARS-CoV-2 Dataset](#3-getting-started-the-sars-cov-2-dataset)
4.  [Core `matUtils` Commands: A Practical Cookbook](#4-core-matutils-commands-a-practical-cookbook)
    *   [Summarizing a Tree](#summarizing-a-tree)
    *   [Extracting Subtrees & Converting Formats](#extracting-subtrees--converting-formats)
    *   [Annotating a Tree](#annotating-a-tree)
    *   [Advanced Analysis](#advanced-analysis)
5.  [Advanced Topics and Key Concepts](#5-advanced-topics-and-key-concepts)
6.  [References](#6-references)

---

## 1. The Mutation-Annotated Tree (MAT)

The MAT is the central data structure used by UShER. It's a phylogenetic tree where branches are annotated with the mutations that occurred along them.

*   **Format:** Stored efficiently in a Protocol Buffers (`.pb`) file.
*   **Content:** Compactly represents massive phylogenies (millions of samples) and the specific genetic changes defining each lineage.
*   **Function:** Allows USher to rapidly identify the optimal placement for a new sequence based on shared mutations.

---

## 2. `matUtils`: The MAT Power Tool

`matUtils` is the essential command-line utility for querying, interpreting, and manipulating MAT files. It is the primary tool for any downstream analysis of UShER trees. Its key functions include **extraction**, **conversion**, **annotation**, and **summarization**.

---

## 3. Getting Started: The SARS-CoV-2 Dataset

A daily-updated MAT for public SARS-CoV-2 sequences is maintained by UCSC and serves as the standard public dataset.

*   **Download Page:** [http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/](http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/)
*   **Primary File:** `public-latest.all.masked.pb.gz`

**Download and Decompression in a Linux/macOS terminal:**
```bash
# Download the latest tree
wget http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.all.masked.pb.gz

# Decompress the file
gunzip public-latest.all.masked.pb.gz
```
The resulting `public-latest.all.masked.pb` file is the input for all `matUtils` commands.

---

## 4. Core `matUtils` Commands: A Practical Cookbook

All commands follow the pattern: `matUtils <subcommand> -i your_tree.pb [options]`.

### Summarizing a Tree

Use `matUtils summary` to get a high-level overview of your MAT file.

```bash
# Get basic counts (nodes, samples, parsimony score)
matUtils summary -i your_tree.pb

# List all samples and their parsimony scores to a file
matUtils summary -i your_tree.pb --samples samples.txt

# List all clades and their sample counts to a file
matUtils summary -i your_tree.pb --clades clades.tsv
```

### Extracting Subtrees & Converting Formats

`matUtils extract` is the most versatile command. It allows you to create smaller, focused datasets and convert them to standard formats.

```bash
# --- Subtree Extraction ---

# Extract all samples of a specific clade into a new MAT file
matUtils extract -i your_tree.pb -c "B.1.1.7" -o B.1.1.7_subtree.pb

# Extract samples containing a specific mutation (e.g., S:E484K)
matUtils extract -i your_tree.pb -m S:E484K -o E484K_subtree.pb


# --- Format Conversion ---

# Convert an entire MAT or a subtree to Newick format (.nwk) for other phylogenetic tools
matUtils extract -i your_tree.pb -t output_tree.nwk

# Convert an entire MAT to VCF (Variant Call Format)
matUtils extract -i your_tree.pb -v all_mutations.vcf

# Extract mutations for just one clade to a VCF
matUtils extract -i your_tree.pb -c "B.1.1.7" -v B.1.1.7_mutations.vcf

# Extract a subtree for interactive visualization in Auspice (JSON format)
# This gets "sample_X" and its 25 nearest neighbors
matUtils extract -i your_tree.pb -k "sample_X:25" -j sample_X_context.json
```

> **Extracting Branch-Specific Mutations:**
> For a detailed list of mutations on *every branch* within a clade (not just the final sample genotypes), use the `--all-paths` flag. This is useful for detailed evolutionary studies.
> ```bash
> matUtils extract -i your_tree.pb -c "B.1.1.7" --all-paths B.1.1.7_branch_muts.txt
> ```

### Annotating a Tree

Use `matUtils annotate` to add or update clade definitions on your tree from a TSV file.

```bash
# The TSV file should have [clade_name] [representative_sample_id]
# Example: B.1.1.7    England/204820464/2020

matUtils annotate -i your_tree.pb -c clade_assignments.tsv -o annotated_tree.pb
```

### Advanced Analysis

`matUtils` also includes subcommands for more specialized analyses.

```bash
# Assess placement certainty for a list of samples
matUtils uncertainty -i your_tree.pb -s query_samples.txt -e uncertainty_scores.tsv

# Infer phylogeographic introduction events (experimental)
matUtils introduce -i your_tree.pb -s sample_regions.tsv -o introductions_report.tsv
```

---

## 5. Advanced Topics and Key Concepts

Mastering `matUtils` involves understanding these key ideas:

*   **Placement Uncertainty:** The `uncertainty` command is crucial for assessing confidence. High **Equally Parsimonious Placements (EPPs)** for a sample means its position on the tree is ambiguous. This is vital for avoiding over-interpretation of transmission events.
*   **Functional Context:** The `summary --translate` option links nucleotide mutations to amino acid changes, providing a window into the functional impact of evolution.
*   **Recurrent Mutations:** The **RoHo (Ratio of Homoplasic Offspring)** score, also from `summary`, helps identify mutations that have appeared independently multiple times across the tree, a potential signal of positive selection.
*   **Bridging Disciplines:** By integrating phylogenetic structure, mutations, and metadata, `matUtils` serves as a critical bridge between genomics, phylogenetics, and epidemiology.

---

## 6. References

*   **UShER Main Paper:** Turakhia, Y., et al. (2021). Ultrafast Sample placement on Existing tRees (UShER) enables real-time phylogenetics... *Nature Genetics*. [doi:10.1038/s41588-021-00862-7](https://doi.org/10.1038/s41588-021-00862-7)
*   **`matUtils` Paper:** McBroome, J., et al. (2021). A daily-updated database and tools for comprehensive SARS-CoV-2 mutation-annotated trees. *Molecular Biology and Evolution*. [doi:10.1093/molbev/msab264](https://doi.org/10.1093/molbev/msab264)
*   **Official Wiki:** [https://usher-wiki.readthedocs.io/](https://usher-wiki.readthedocs.io/)
