# UShER Toolkit and matUtils: A Condensed Guide

## 1. Introduction to UShER

UShER (Ultrafast Sample placement on Existing tRee) is a phylogenetic toolkit crucial for real-time genomic epidemiology. It enables rapid placement of new pathogen sequences (like SARS-CoV-2) onto massive existing phylogenetic trees, helping to track viral evolution and spread. Its speed and scalability are vital for timely public health responses.

## 2. The Mutation Annotated Tree (MAT) Format

The Mutation Annotated Tree (MAT) is a phylogenetic tree format where branches are annotated with the mutations inferred to have occurred along them. This format, typically stored as a Protocol Buffers (`.pb`) file, efficiently represents large phylogenies and the genetic changes defining each lineage. It's central to UShER's operations, allowing for rapid searching and updating of the tree. Key characteristics include:

*   **Efficiency:** Stores large trees and extensive mutation data compactly.
*   **Annotation Rich:** Each node and branch can carry detailed mutation information.
*   **Foundation for UShER:** UShER processes query sequences by identifying the optimal placement based on shared mutations within the MAT.

## 3. matUtils: Overview and Core Functions

`matUtils` is a command-line toolkit designed for querying, interpreting, and manipulating MAT files. It is an essential companion to UShER for downstream analysis. Core functions include:

*   **Annotation:** Adding metadata or clade definitions (e.g., Pango lineages) to the tree.
*   **Extraction:** Creating subtrees based on various criteria (e.g., specific clades, samples, or proximity to a query sequence). This is useful for focusing on particular parts of the phylogeny.
*   **Conversion:** Transforming MAT files into other standard phylogenetic or variant formats, such as:
    *   Newick (`.nh`): For use with many standard phylogenetic viewers.
    *   VCF (Variant Call Format): For analysis of mutations.
*   **Summarization:** Generating summary statistics about the tree, such as clade sizes or mutation distributions.
*   **Introduction:** Facilitating the integration of new mutations or samples into an existing MAT, although `usher` is the primary tool for de novo placement.
*   **Optimization:** While `matOptimize` is a separate tool, `matUtils` can help prepare or analyze trees for such parsimony optimization steps.

## 4. Accessing the Latest SARS-CoV-2 Dataset

A daily-updated, pre-processed mutation-annotated tree (MAT) for public SARS-CoV-2 sequences is maintained by UCSC. The key file to download is `public-latest.all.masked.pb.gz`.

*   **Download Link:** [http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/](http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/)
*   **Specific File:** `public-latest.all.masked.pb.gz` (This is the main MAT file, compressed)
*   This directory also contains associated metadata and VCF files.

**Download and Decompression:**
You can download the file using `wget`:
```bash
wget http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.all.masked.pb.gz
```
After downloading, decompress the `.gz` file using `gunzip`:
```bash
gunzip public-latest.all.masked.pb.gz
```
This will result in the `public-latest.all.masked.pb` file, which is used as input for `matUtils`.

For visualization of this global tree and its genotypes, [https://cov2tree.org/](https://cov2tree.org/) (with treenome enabled) is a useful resource.

## 5. Practical Procedures for Using matUtils

`matUtils` offers a suite of subcommands for common tasks. All commands require an input MAT file (e.g., `-i path/to/your_mat.pb`). Here are some practical examples:

*   **Summarizing Tree Statistics:**
    *   Get basic counts (nodes, samples, parsimony): `matUtils summary -i path/to/your_mat.pb`
    *   List all samples and their parsimony scores: `matUtils summary -i path/to/your_mat.pb --samples samples.txt`
    *   Get amino acid translations for mutations: `matUtils summary -i path/to/your_mat.pb -t translations.tsv -g genes.gtf -f reference.fasta`

*   **Extracting Subtrees and Converting Formats:**
    *   Extract all samples belonging to a specific clade (e.g., "B.1.1.7") into a new MAT file, using the decompressed global tree as input:
        `matUtils extract -i public-latest.all.masked.pb -c B.1.1.7 -o B.1.1.7_subtree.pb`
    *   Convert an entire MAT (e.g., the decompressed global tree) or a subtree to a Newick tree file:
        `matUtils extract -i public-latest.all.masked.pb -t output_tree.nwk`
        This uses the decompressed `public-latest.all.masked.pb` as input. The Newick format (.nwk) is a standard text-based format ideal for use with other phylogenetic software like MEGA, IQ-TREE, RAxML, or visualization tools like FigTree.
    *   Extract samples containing a specific mutation (e.g., "S:E484K") and output as VCF:
        `matUtils extract -i public-latest.all.masked.pb -m S:E484K -v E484K_samples.vcf`
    *   Extract all mutations for all samples from the MAT into a VCF file:
        `matUtils extract -i public-latest.all.masked.pb -v all_mutations.vcf`
        This VCF lists all mutations for every sample relative to the reference genome used in the MAT's construction. The MAT format stores mutations on branches; the resulting VCF effectively translates these paths of mutations from the root to each sample into a standard genotype format.
        *Note:* For highly specific analyses requiring explicit per-branch mutation lists in a custom format, one might need custom scripting (e.g., Python with dendropy/ete3 and pyVCF). However, for most applications, the MAT itself (queried with `matUtils`) or the comprehensive VCF output is sufficient.
    *   Extract a subtree around a specific sample ("sample_X") including 25 nearest neighbors into a JSON for Auspice:
        `matUtils extract -i public-latest.all.masked.pb -k "sample_X:25" -j sample_X_context.json`

*   **Annotating Trees:**
    *   Add clade annotations to your MAT file based on a TSV file (`clade_assignments.tsv`) where column 1 is clade name and column 2 is a representative sample ID for that clade:
        `matUtils annotate -i path/to/your_mat.pb -c clade_assignments.tsv -o annotated_tree.pb`
    *   Clear existing annotations before applying new ones:
        `matUtils annotate -i path/to/your_mat.pb -c clade_assignments.tsv -o annotated_tree.pb --clear-current`

*   **Assessing Placement Uncertainty:**
    *   Calculate Equally Parsimonious Placements (EPPs) and Neighborhood Size Score (NSS) for a list of samples (`query_samples.txt`):
        `matUtils uncertainty -i path/to/your_mat.pb -s query_samples.txt -e uncertainty_scores.tsv`

*   **Analyzing Introductions (Experimental Heuristic):**
    *   Infer introduction events into geographic regions based on a TSV file (`sample_regions.tsv`) mapping samples to regions:
        `matUtils introduce -i path/to/your_mat.pb -s sample_regions.tsv -o introductions_report.tsv`

**General Workflow:**
1.  **Obtain MAT:** Download the latest public tree (e.g., for SARS-CoV-2) or generate one using `usher`.
2.  **Explore:** Use `matUtils summary` to understand its basic properties.
3.  **Manipulate/Query:**
    *   Use `matUtils extract` to get subtrees of interest (by clade, mutation, sample list, etc.) and convert them to VCF, Newick, or JSON.
    *   Use `matUtils annotate` to add custom clade definitions.
    *   Use `matUtils uncertainty` to check placement robustness for key samples.
    *   Use `matUtils introduce` for phylogeographic insights.
4.  **Visualize/Downstream:** Use outputs with tools like Auspice (via JSON), Nextclade, MicrobeTrace, or standard phylogenetic software.

## 6. Key Understandings and Advanced Applications

Beyond routine tasks, `matUtils` enables more sophisticated analyses and provides deeper insights into pathogen evolution and epidemiology:

*   **Understanding Tree Structure and Certainty:**
    *   The MAT format itself (a protobuf) is highly efficient for storing massive phylogenies with full mutation annotations. This is key to handling datasets like the global SARS-CoV-2 tree.
    *   `matUtils uncertainty` (calculating EPPs and NSS) is crucial for assessing the confidence of individual sample placements. High EPPs or large NSS values indicate ambiguity, which is important when interpreting potential transmission links or identifying novel variants. This helps avoid over-interpreting tree topology.

*   **Customized Data Extraction and Integration:**
    *   The `extract` command is powerful. Its various filtering options (by clade, mutation, metadata, parsimony score, branch length) allow for highly specific dataset generation for targeted research questions.
    *   Ability to output to VCF allows integration with variant analysis pipelines. Newick output facilitates use with a wide array of phylogenetic tools. JSON output is tailored for Nextstrain/Auspice, enabling rich, interactive visualizations.

*   **Advanced Annotation and Interpretation:**
    *   `matUtils annotate` allows for dynamic updating of clade definitions on the tree, essential as new variants of concern emerge or classification systems (like Pango lineages) are updated. It can infer clade roots from representative samples or assign them directly.
    *   The `--translate` option in `summary` links nucleotide mutations to amino acid changes, providing functional context.
    *   The RoHo (Ratio of Homoplasic Offspring) score calculation in `summary` can help identify recurrent mutations potentially under positive selection, flagging mutations of interest for further investigation.

*   **Phylogeographic Analysis:**
    *   `matUtils introduce`, while partly heuristic, provides methods (including Association Index and Max Monophyletic Clade size) to estimate viral introductions into geographic regions. This can help quantify transmission dynamics and assess the impact of travel or public health interventions. The confidence metrics help gauge the strength of these inferences.

*   **Scalability and Automation:**
    *   `matUtils` is designed to work with very large trees (millions of samples) efficiently, often on standard computing hardware.
    *   Its command-line nature makes it suitable for scripting and integration into automated bioinformatics pipelines for routine surveillance reporting.

*   **Bridging Phylogenetics and Epidemiology:**
    *   By allowing mutations, clade assignments, sample metadata (dates, locations), and phylogenetic structure to be jointly analyzed and manipulated, `matUtils` serves as a critical bridge between raw sequence data, phylogenetic inference, and epidemiological interpretation.

*   **Considerations for Advanced Use:**
    *   **Input Data Quality:** The quality of the input MAT (and the underlying sequence alignments and tree construction) significantly impacts `matUtils` results. Masking problematic sites in the reference genome is crucial.
    *   **Parameter Choice:** Many `matUtils` commands have numerous options. Understanding their specific effects is important for obtaining meaningful results (e.g., allele frequency thresholds in `annotate`, confidence cutoffs in `introduce`).
    *   **Interpretation Context:** Phylogenetic results, even with `matUtils`, require careful interpretation in the context of sampling biases, sequencing coverage, and available epidemiological data.

By mastering `matUtils`, researchers and public health professionals can more effectively navigate and utilize the vast genomic data generated during pandemics like COVID-19.

## 7. Advanced: Extracting and Processing Branch-Specific Mutations

While `matUtils extract -v` provides VCF files detailing mutations for each sample (leaf node) relative to the reference, sometimes a more explicit list of mutations occurring on each internal branch of the tree is desired. UShER's MAT file inherently contains this information, as it maps parsimony-inferred mutations to each branch.

### Prioritizing `matUtils dump -m` for Branch Mutations

The most direct and computationally efficient method to obtain a list of mutations for each branch is using the `matUtils dump` subcommand with the `-m` option. This command directly accesses the branch mutation data stored within the MAT file.

**Command:**
```bash
matUtils dump -i path/to/your_mat.pb -m usher_branch_mutations.tsv
```
*   `-i path/to/your_mat.pb`: Specifies your input MAT file (e.g., `public-latest.all.masked.pb`).
*   `-m usher_branch_mutations.tsv`: Specifies the output file where branch mutations will be written.

**Expected Output (`usher_branch_mutations.tsv`):**
The output is typically a tab-separated file. Each line often represents a node (or the branch leading to it) and lists the mutations inferred on that specific branch relative to its parent. The exact format can vary slightly with `matUtils` versions, but it generally includes:
*   A node identifier.
*   A list of mutations (e.g., `A123T,G456C`) that occurred on the branch leading to this node.

This output provides a clear, direct list of mutations per branch as inferred by UShER's parsimony algorithm.

### Using the `process_usher_branch_mutations.py` Script

The `usher_branch_mutations.tsv` file from `matUtils dump -m` is often sufficient. However, you might want to integrate this mutation data with a tree structure in Python for custom analyses, visualizations, or to associate mutations with specific nodes in a traversable tree object. The following conceptual script demonstrates how to do this using `dendropy` for tree manipulation and `pandas` for parsing the mutation file.

An executable Python script version of this conceptual procedure is available in this repository: [`process_usher_branch_mutations.py`](scripts/process_usher_branch_mutations.py). The script includes comments outlining the necessary `matUtils` commands to generate its input files (`usher_tree.nwk` and `usher_branch_mutations.tsv`). The 'Required Prior Steps' and 'Note on Node ID Mapping' described below provide essential context for using the script.

**Required Prior Steps (using `matUtils`):**

1.  **Download and decompress the latest MAT file:**
    ```bash
    wget http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.all.masked.pb.gz
    gunzip public-latest.all.masked.pb.gz
    ```
    (This provides `public-latest.all.masked.pb`)

2.  **Export the Newick tree:**
    ```bash
    matUtils extract -i public-latest.all.masked.pb -t usher_tree.nwk
    ```

3.  **Dump the branch mutations:**
    ```bash
    matUtils dump -i public-latest.all.masked.pb -m usher_branch_mutations.tsv
    ```

**Note on Node ID Mapping:** The critical part of the script is correctly mapping the `Node_ID` from `matUtils dump -m` to the nodes in the `dendropy` tree object. Leaf nodes are usually straightforward using `node.taxon.label`. Internal nodes can be more challenging if their labels are not consistent between `matUtils extract`'s Newick output and `matUtils dump`'s node identifiers. Further customization of the script might be needed based on the specifics of these identifiers in your dataset.

### Alternative (More Complex): Ancestral Reconstruction from Newick + VCF

If `matUtils dump -m` were unavailable, or if one wished to use a different ancestral state reconstruction algorithm, it's conceptually possible to infer branch mutations by:
1.  Parsing a Newick tree (`matUtils extract -t`).
2.  Parsing a VCF file containing tip (sample) genotypes (`matUtils extract -v`).
3.  Using libraries like `pysam` (for VCF) and `dendropy` (for tree traversal and parsimony algorithms like Fitch's algorithm) to reconstruct ancestral states at each internal node for each polymorphic site.
4.  Comparing parent and child node genotypes along the tree to identify mutations on each branch.

This approach is significantly more complex to implement correctly and is computationally far more intensive than using UShER's pre-computed parsimony inferences via `matUtils dump -m`. For large trees like the global SARS-CoV-2 phylogeny, re-implementing this in Python would be impractical for routine use compared to leveraging UShER's optimized C++ functionalities.

---

## 8. References

Primary UShER package papers:

*   Turakhia, Y., De Maio, N., Thornlow, B. et al. Ultrafast Sample placement on Existing tRees (UShER) enables real-time phylogenetics for structured pathogen surveillance. *Nat Genet* **53**, 839–846 (2021). [https://doi.org/10.1038/s41588-021-00862-7](https://doi.org/10.1038/s41588-021-00862-7)
*   McBroome, J., Thornlow, B., Hinrichs, A.S. et al. A daily-updated database and tools for comprehensive SARS-CoV-2 mutation-annotated trees. *Mol Biol Evol* **38**, 5469–5474 (2021). [https://doi.org/10.1093/molbev/msab264](https://doi.org/10.1093/molbev/msab264) (matUtils)

Further information:

*   UShER Website: [https://usher.bio/](https://usher.bio/)
*   UShER GitHub Repository & Wiki: [https://github.com/yatisht/usher](https://github.com/yatisht/usher)
