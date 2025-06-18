# Ancestral Reconstruction Tutorial

Modern virology research often involves analyzing vast datasets, encompassing thousands to millions of viral sequences. This scale presents significant challenges for traditional phylogenetic methods. This tutorial will guide you through the process of ancestral reconstruction, focusing on a hybrid strategy. This strategy combines the capabilities of Python libraries like DendroPy, BioPython, and ETE Toolkit for data manipulation and basic phylogenetic tasks, with specialized high-performance tools for computationally intensive steps like large-scale tree building and the ancestral reconstruction itself (e.g., using TreeTime or DendroPy's own ASR functions).

Ancestral reconstruction is a powerful technique used to infer the genetic sequences or other character states of ancestral organisms, providing insights into viral evolution, adaptation, and the emergence of new traits.

## Foundational Concepts and Tool Selection

Performing ancestral reconstruction, especially at scale, requires careful tool selection based on:

*   **Scalability:** Methods must efficiently handle large alignments and numerous trees. This involves both algorithmic efficiency and practical implementation for parallel processing or handling large data structures.
*   **Accuracy:** The reliability of inferred ancestral states is paramount. This depends on appropriate evolutionary models, quality of input data (alignment, tree), and the reconstruction algorithm.
*   **Flexibility and Programmability:** The ability to script, customize, and integrate different tools into a cohesive workflow is crucial for complex analyses. Python libraries play a significant role here.

**Key Components of a Hybrid Strategy:**

*   **Python Libraries for Phylogenetics:**
    *   **BioPython:** A cornerstone for bioinformatics in Python. Used for sequence manipulation (reading/writing various formats, translation, etc.), accessing online databases (e.g., NCBI Entrez for sequence retrieval), and basic phylogenetic tree operations.
        ```python
        # Example: Fetching sequences with Bio.Entrez
        from Bio import Entrez
        from Bio import SeqIO
        Entrez.email = "Your.Name.Here@example.org" # Always tell NCBI who you are
        handle = Entrez.efetch(db="nucleotide", id=["AY884001", "AY884002"], rettype="fasta", retmode="text")
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()
        # print(records[0].id, records[0].seq[:10])
        ```
    *   **DendroPy:** A Python library specifically for phylogenetic computing. It offers sophisticated capabilities for reading, writing, and manipulating phylogenetic trees and character matrices. Critically, DendroPy includes its own functions for ancestral state reconstruction (both maximum parsimony and maximum likelihood), which can be an alternative or complement to tools like TreeTime, especially for custom analyses or when working entirely within a Python environment.
    *   **ETE Toolkit (Environment for Tree Exploration):** Another powerful Python library for tree manipulation, annotation, and visualization. ETE is excellent for programmatic tree rendering and comparing different tree topologies or annotations.

*   **Specialized Command-Line Tools:**
    *   **Sequence Alignment:** MAFFT, Clustal Omega, MUSCLE. For very large datasets, consider tools that support iterative refinement or adding sequences to existing alignments.
    *   **Phylogenetic Inference:** IQ-TREE, RAxML, FastTree. The choice often depends on the dataset size, desired accuracy, and available computational resources. For extremely large datasets, tools like UShER (for placing new sequences onto a backbone tree) can be invaluable.
    *   **Ancestral Reconstruction & Time-Scaling:**
        *   **TreeTime:** Excellent for joint inference of ancestral sequences, time-scaled trees (when dates are provided), and evolutionary rates. Its command-line interface is convenient for integration into scripts.
        *   **DendroPy's ASR:** Provides functions like `dendropy.Tree.reconstruct_ancestral_states()` which can be used programmatically for maximum parsimony or likelihood ASR. This is useful for deeper integration into Python-based workflows.

*   **Workflow Management:** For large-scale, multi-step analyses, workflow managers like Snakemake or Nextflow help automate, parallelize, and ensure reproducibility.

## Prerequisites

Before you begin, ensure you have the following installed:

*   **Python:** Version 3.7 or higher.
*   **Core Python Libraries:**
    *   `pip install biopython dendropy ete3`
*   **TreeTime (Optional, but a focus of this tutorial for command-line ASR):**
    *   `pip install phylotreetime`
*   **MAFFT:** For multiple sequence alignment (see MAFFT website for installation).
*   **IQ-TREE or other tree building software:** (see their respective websites).
*   **(Optional) For very large datasets (subsampling/clustering):**
    *   **MMseqs2:** [https://github.com/soedinglab/MMseqs2](https://github.com/soedinglab/MMseqs2)
    *   **CD-HIT:** [http://weizhongli-lab.org/cd-hit/](http://weizhongli-lab.org/cd-hit/)

You will also need:

*   A set of sequences (e.g., FASTA format). This tutorial will show an example of fetching them using BioPython.
*   A phylogenetic tree file in Newick or Nexus format if you are not inferring it as part of the workflow.
*   **Metadata (Highly Recommended for TreeTime):** A CSV or TSV file containing sampling dates for your sequences if you intend to create a time-scaled tree with TreeTime. Format:
    ```
    name,date
    sequence_id1,2020.12
    sequence_id2,2021.03
    # Or for more complex date formats TreeTime can parse:
    # sequence_id3,2020-03-15
    ```

## Step 1: Preparing your data (Bulk Data Preparation)

Efficient data handling and preparation are critical, especially for large datasets.

**1. Acquiring Sequences (Example with BioPython):**

If your sequences are not already in a file, you can fetch them from NCBI using BioPython's Entrez module.

```python
from Bio import Entrez
from Bio import SeqIO
import os

# Always provide your email to NCBI
Entrez.email = "your.email@example.com"

# List of GenBank accession numbers
accession_numbers = ["MN908947.3", "MT291826.1", "LR890198.1"] # Example: SARS-CoV-2 sequences
output_fasta = "retrieved_sequences.fasta"

print(f"Fetching {len(accession_numbers)} sequences...")
sequences_to_write = []
for accno in accession_numbers:
    try:
        handle = Entrez.efetch(db="nucleotide", id=accno, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        sequences_to_write.append(record)
        handle.close()
        print(f"Fetched {accno}")
    except Exception as e:
        print(f"Error fetching {accno}: {e}")

# Write to a single FASTA file
SeqIO.write(sequences_to_write, output_fasta, "fasta")
print(f"All fetched sequences saved to {output_fasta}")
```
*Remember to replace `"your.email@example.com"` with your actual email address.*

**2. Subsampling/Clustering for Very Large Datasets (Optional):**

If you have tens of thousands of sequences or more, aligning them all and building a comprehensive tree can be computationally prohibitive or unnecessary for certain analyses.

*   **Strategy:** Cluster sequences by similarity and select representatives from each cluster.
*   **Tools:**
    *   **MMseqs2:** Fast and sensitive clustering and searching of sequence sets.
        ```bash
        # Example: Cluster sequences at 99% identity, keep representatives
        mmseqs easy-cluster retrieved_sequences.fasta cluster_results tmp_mmseqs --min-seq-id 0.99 -c 0.8 --cov-mode 1
        # The representative sequences will be in cluster_results_rep_seq.fasta
        # Check MMseqs2 documentation for up-to-date commands and options.
        ```
    *   **CD-HIT:** Another popular tool for sequence clustering.
        ```bash
        cd-hit-est -i retrieved_sequences.fasta -o representative_sequences.fasta -c 0.99 -n 10 -d 0 -M 16000
        # -c 0.99: 99% identity threshold
        # -n 10: word size
        ```
    Choose one set of representative sequences (e.g., `representative_sequences.fasta`) for the next steps. Let's assume this is now `aligned_sequences.fasta` for subsequent steps if subsampling was done, otherwise `retrieved_sequences.fasta` would be used.

**3. Multiple Sequence Alignment (MSA):**

Align your (potentially subsampled) sequences. MAFFT is a common choice.

*   Standard MAFFT example (using `retrieved_sequences.fasta` or your subsampled FASTA):
    ```bash
    # Using --auto selects appropriate strategy (e.g. FFT-NS-2 for speed with many sequences)
    # --thread -1 uses all available cores
    mafft --auto --thread -1 representative_sequences.fasta > aligned_sequences.fasta
    ```
*   Ensure your MSA is in FASTA format and that sequence names are consistent.

**4. Phylogenetic Tree Inference:**

Infer a tree from your MSA. IQ-TREE is robust and feature-rich.

*   Example with IQ-TREE:
    ```bash
    # Using a model like GTR+F+G4 is common for viruses.
    # -nt AUTO lets IQ-TREE determine optimal threads.
    # Use a distinct prefix for output files, e.g., 'iqtree_run'
    iqtree -s aligned_sequences.fasta -m GTR+F+G4 -nt AUTO --threads-max AUTO -mem 16G -pre iqtree_run
    # The main tree file will be iqtree_run.treefile
    ```
*   The leaf names in the tree **must exactly match** the sequence names in `aligned_sequences.fasta` and your metadata file (if used).
*   For TreeTime, providing sampling dates via a metadata file to get a time-scaled tree is highly recommended. If dates are not available, TreeTime will produce a tree scaled by divergence (substitutions per site).
The output tree (e.g., `iqtree_run.treefile`) will be used in the next step.

## Step 2: Ancestral Sequence Reconstruction (ASR)

Once you have your alignment and phylogenetic tree (e.g., `aligned_sequences.fasta` and `iqtree_run.treefile`), you can perform ASR. We'll cover two main approaches: TreeTime (command-line) and DendroPy (Python library).

**A. Using TreeTime (Command-Line for ASR and Time-Scaling)**

TreeTime is particularly strong when you also want to infer a time-scaled tree along with ancestral sequences.

**Basic TreeTime Command:**
Provide your alignment, tree, and (optionally but recommended) a metadata file with dates.

```bash
# Ensure your metadata file (e.g., dates.csv) has 'name,date' columns
# matching sequence names in aligned_sequences.fasta
treetime ancestral --aln aligned_sequences.fasta \
                   --tree iqtree_run.treefile \
                   --dates dates.csv \
                   --outdir treetime_asr_output
```

*   `--aln`: Path to your MSA file (e.g., `aligned_sequences.fasta`).
*   `--tree`: Path to your Newick tree file (e.g., `iqtree_run.treefile`).
*   `--dates`: Path to your CSV/TSV metadata file with sampling dates. This enables time-scaling.
*   `--outdir`: Directory where TreeTime will save the results.

**Key Outputs from TreeTime:**

*   `annotated_tree.nexus` (in `--outdir`): The tree with inferred ancestral sequences and mutations annotated on branches. Viewable in FigTree.
*   `ancestral_sequences.fasta` (in `--outdir`): FASTA file of inferred sequences for internal nodes.
*   `molecular_clock.txt` (in `--outdir`): Estimated clock rate, TMRCA, etc.
*   `gtr.txt` (or similar, in `--outdir`): Parameters of the substitution model.
*   (Optional) `report.html`: An HTML report if you add the `--report` flag.

**B. Using DendroPy (Programmatic ASR in Python)**

DendroPy allows you to perform ASR directly within your Python scripts, offering fine-grained control. It supports both maximum parsimony and maximum likelihood ASR.

```python
import dendropy
from dendropy.model.discrete import Jc69 # Example substitution model

# Load the character matrix (alignment)
# Ensure this alignment corresponds to the tips of your tree
char_matrix = dendropy.DnaCharacterMatrix.get(
    path="aligned_sequences.fasta",
    schema="fasta"
)

# Load the tree
# Make sure taxon namespace of tree matches character matrix
tree = dendropy.Tree.get(
    path="iqtree_run.treefile",
    schema="newick",
    taxon_namespace=char_matrix.taxon_namespace # Important!
)

# DendroPy can perform both Maximum Parsimony (MP) and Maximum Likelihood (ML) ASR.
# ML ASR is generally more accurate but computationally more intensive and requires
# specifying a substitution model. MP ASR is faster and simpler.

# Example: Maximum Parsimony ASR with DendroPy
# The `reconstruct_ancestral_states()` method annotates the tree nodes in place
# with the inferred ancestral states.
print("Running Maximum Parsimony ASR with DendroPy...")
tree.reconstruct_ancestral_states(
    character_matrix=char_matrix,
    method="parsimony",  # Use "likelihood" for ML ASR
    # For ML ASR, you'd also need to pass a `submodel` object:
    # submodel=dendropy.model.discrete.Jc69() # Or other models like Hky85, Gtr
)

# After running ASR, ancestral states are stored in node annotations.
# You can iterate through nodes and access these states.
# For parsimony, this might be a set of equally parsimonious states.
# For ML, this would typically be the state with the highest probability.

print("\nAncestral states (parsimony) at internal nodes (first 10 sites for brevity):")
for node in tree.internal_nodes():
    if node.parent_node is None: # Root node
        label = "Root"
    else:
        label = node.label if node.label else "Internal"

    # The actual annotation attribute might depend on DendroPy version and ASR method.
    # It's often `node.ancestral_state_annotations` or similar.
    # This is a conceptual way to access and print:
    # We'll assume states are on `node.annotations.get('ancestral_states')` or similar
    # and that it's a list of character state objects for each site.

    # The following is a placeholder for how you might access character data.
    # Actual extraction depends on how DendroPy stores these results.
    # Typically, you'd access node.annotations or a specific attribute set by ASR.
    # For this example, we'll just print a message.
    # To get a full sequence string, you'd concatenate these character states.

    # Example: Accessing parsimony result (often a set of states)
    # This is illustrative. Check DendroPy docs for precise API for state extraction.
    # reconstructed_states_at_node = []
    # if hasattr(node, 'your_asr_annotation_attribute'): # Replace with actual attribute
    #    for site_idx in range(min(10, char_matrix.sequence_size)):
    #       state_at_site = node.your_asr_annotation_attribute.get_state(site_idx) # Conceptual
    #       reconstructed_states_at_node.append(str(state_at_site))
    #    print(f"Node {label}: {''.join(reconstructed_states_at_node)}...")
    pass # Pass for now as detailed extraction is complex for a brief example.

print("\nNote: Extracting and formatting full ancestral sequences from DendroPy's ASR results")
print("requires iterating through node annotations and character states.")
print("Please consult the DendroPy documentation for the precise API and examples,")
print("especially for `reconstruct_ancestral_states_likelihood` for ML ASR.")

# Save the tree with ASR annotations (DendroPy might add them as custom annotations)
# tree.write(path="dendropy_annotated_tree_parsimony.newick", schema="newick")

# For ML ASR with DendroPy:
# 1. Define a substitution model (e.g., Jc69, Hky85, Gtr).
#    from dendropy.model.discrete import Hky85
#    sub_model = Hky85(kappa=2.0) # Example: HKY model with kappa = 2.0
# 2. Run `tree.reconstruct_ancestral_states(char_matrix, method="likelihood", submodel=sub_model)`
# 3. Process results similarly, but expect probability distributions or most likely states.
# ML ASR is significantly more complex to configure correctly than parsimony.
# For robust, out-of-the-box ML ASR, dedicated tools like TreeTime, PAML, PhyML, or FastML are often preferred
# if a command-line interface is suitable. DendroPy provides the building blocks for ML ASR
# if you need it within a Python-driven workflow and are prepared for more detailed setup.
```

**Handling Ambiguities in ASR:**

*   **Input Data:** Ambiguous characters (N, R, Y, etc.) in your input MSA can influence ASR. Most tools will try to resolve them based on the phylogenetic context or treat them as partial information.
*   **Output:**
    *   **TreeTime:** Typically outputs the most probable character at each ancestral site. Confidence scores for these reconstructions can often be found in supplementary output files or tree annotations.
    *   **DendroPy (Parsimony):** Fitch's algorithm (common in parsimony) can result in multiple equally parsimonious states at a node (e.g., a site could be A or G). The output might be a set of states. For ML, it usually gives probabilities for each state.
    *   **Interpretation:** Be aware of sites with low confidence or high ambiguity in ancestral sequences. These might be regions of rapid evolution, recombination, or insufficient data.

## Step 3: Strategies for Very Large Datasets & Scalability

Standard ASR can be slow on trees with tens of thousands of tips or more.

*   **Subsampling/Clustering:** Already discussed in Step 1 (MMseqs2/CD-HIT). This is the primary way to reduce dataset size before tree building and ASR.
*   **Divide and Conquer:**
    1.  Build a tree for the full dataset (or a large representative subset).
    2.  Identify major clades or subtrees.
    3.  Extract these subtrees. Both DendroPy and ETE Toolkit can do this.
        *   **DendroPy Example (conceptual):**
            ```python
            # Assuming 'full_tree' is a DendroPy Tree object and 'clade_node_label' is the label of the MRCA of your desired clade
            # mrca_node = full_tree.find_node_with_label(clade_node_label)
            # if mrca_node:
            #    clade_tree = dendropy.Tree(seed_node=mrca_node, taxon_namespace=full_tree.taxon_namespace)
            #    clade_tree.is_rooted = True # or False, depending on context
            #    # Prune the character matrix to match the clade_tree
            #    clade_char_matrix = char_matrix.extract_taxa(taxon_labels=clade_tree.infer_taxa().labels())
            #    # Now perform ASR on clade_tree and clade_char_matrix
            ```
        *   **ETE Toolkit Example:**
            ```python
            from ete3 import Tree
            # full_ete_tree = Tree("iqtree_run.treefile", format=1) # format=1 for internal node names
            # Assuming 'NODE_X' is the name of your MRCA node in Newick file
            # clade_ancestor_node = full_ete_tree.search_nodes(name="NODE_X")[0]
            # clade_ete_tree = clade_ancestor_node.detach() # Detaches the subtree
            # # Get list of leaf names in the new subtree
            # clade_leaf_names = clade_ete_tree.get_leaf_names()
            # # You would then need to filter your alignment to include only these leaf names
            # # and then run ASR on clade_ete_tree and the filtered alignment.
            # clade_ete_tree.write(outfile="clade_subtree.newick")
            ```
    4.  Perform ASR separately on each subtree (and its corresponding filtered alignment). This can be parallelized.
    5.  Optionally, ASR can be run on a "backbone" tree (the main tree with major clades collapsed or represented by single sequences) to understand deeper relationships.
    *   **Considerations:** Defining appropriate clades (monophyletic groups), ensuring correct extraction of corresponding alignment data, and potentially re-optimizing models for sub-analyses are key. This approach is more complex but can make very large analyses tractable.

*   **Hardware & Parallelization:**
    *   Utilize multi-core processors. Most tree-building tools (IQ-TREE) and some ASR steps can be parallelized.
    *   Ensure sufficient RAM, especially for large alignments and trees.

*   **Tool-Specific Scalability:**
    *   **TreeTime:** Reasonably efficient for its joint inference. The most time-consuming part is often the initial tree building, not TreeTime itself if the tree is pre-calculated.
    *   **DendroPy:** Performance depends on the chosen algorithm (parsimony is generally faster than ML) and the size of the tree/matrix. As a Python library, very large operations might be slower than optimized C/C++ command-line tools unless specific parts are Cythonized or call underlying C libraries.

## Step 4: Visualizing and Interpreting Ancestral Sequences

Visualization is key to understanding ASR results.

1.  **Using FigTree (or similar like iTOL, IcyTree):**
    *   Open FigTree.
    *   Go to `File > Open` and select `annotated_tree.nexus` from your TreeTime output (or a similarly annotated tree from DendroPy if you save it in Nexus format with annotations).
    *   **Displaying Mutations/States:**
        *   In FigTree, look under "Branch Labels" or "Node Labels". You might need to enable display of specific annotations (e.g., "mutations", "states", or custom labels you/the tool added).
        *   TreeTime often annotates mutations on branches (e.g., `A123G` meaning A at site 123 changed to G).
        *   Ancestral states (the full sequence or important sites) can be displayed at nodes if the tree file contains them.
    *   **Coloring:** Use FigTree's features to color branches or nodes by inferred states, mutations, or confidence scores to highlight patterns.

2.  **Interpreting `ancestral_sequences.fasta` (from TreeTime or custom script):**
    *   This file contains the most probable sequence for each internal node (e.g., `NODE_0000001`).
    *   These can be used for:
        *   Identifying specific amino acid or nucleotide changes along evolutionary pathways.
        *   Synthesizing ancestral genes for experimental validation (though this is highly advanced).
        *   Comparing ancestral states at different points in the tree.

3.  **Programmatic Analysis (DendroPy, ETE Toolkit, BioPython):**
    *   Use these libraries to parse tree files with ASR annotations.
    *   Traverse the tree, extract ancestral states/sequences for specific nodes.
    *   Automate the counting of specific types of mutations or changes along particular lineages.
    *   Generate custom plots or reports.

    ```python
    # Example: Basic tree traversal with ETE Toolkit to access node names
    # (Assuming you have a Newick tree, perhaps annotated)
    from ete3 import Tree
    # t = Tree("iqtree_run.treefile") # Load your tree
    # for node in t.traverse("postorder"): # or "preorder"
    #   if not node.is_leaf():
    #       print("Node:", node.name)
    #       # If your tree has ASR information stored as features, access them:
    #       # if hasattr(node, "ancestral_sequence"):
    #       #    print(node.ancestral_sequence[:10])
    ```

## Step 5: Advanced Options (TreeTime Specific) & Further Considerations

This section largely retains the TreeTime advanced options from the previous version, as they are specific to that tool. If using DendroPy for ML ASR, analogous options would involve selecting substitution models, optimizing parameters, etc., within DendroPy's API.

TreeTime offers many options to customize your ancestral reconstruction:
*   **Time-scaled trees (`--dates`):** As highlighted, crucial for temporal inference.
*   **Substitution Model (`--model`):** (e.g., `HKY`, `JC69`, `GTR` for nucleotides; various protein models with `--aa`).
    ```bash
    treetime ancestral ... --model HKY
    treetime ancestral --aa ... --model WAG
    ```
*   **Clock Model and Rate (`--clock-rate`, `--clock-std-dev`, `--vary-rate`):** Control clock assumptions.
*   **Coalescent Skyline (`--coalescent skyline` or other priors):** For estimating effective population size over time with TreeTime. This requires a time-scaled tree.
*   **Confidence Estimates (`--confidence`):** In TreeTime, this typically refers to confidence intervals for dates on the time-scaled tree. For ancestral sequence confidence, TreeTime's output might include probabilities for each character at each site in supplementary files or directly in the annotated Nexus tree (e.g., as `posterior_probabilities` which can be parsed).
*   **Keep polytomies (`--keep-polytomies`):** Preserves polytomies if desired, otherwise TreeTime tries to resolve them.
*   **Rooting (`--root`):** Specify a root for the tree (e.g., `best`, `least-squares`, a specific sequence name, or an outgroup).
*   **Filter by number of mutations/branch length:** (`--max-branch-length`, etc.) Useful for QC.

Refer to the TreeTime documentation (`treetime ancestral --help` or the online documentation) for a full list of options.

**Further Considerations (General for ASR):**
*   **Model Selection:** Choosing an appropriate substitution model (for nucleotides or amino acids) is critical for accurate ML ASR. Tools like ModelFinder (often integrated into IQ-TREE: `iqtree -s alignment -m TEST` or `MFP`) can help select the best-fitting model based on your alignment.
*   **Uncertainty:** ASR results always have inherent uncertainty. Be cautious about over-interpreting a single reconstructed state, especially if confidence is low or multiple states are nearly equally probable (check posterior probabilities if available). Visualizing this uncertainty can be important.
*   **Data Quality:** The quality of your multiple sequence alignment and phylogenetic tree heavily impacts ASR. "Garbage in, garbage out." Ensure your alignment is reliable and your tree is well-supported.
*   **Iterative Refinement:** Large-scale phylogenetic analysis, including ASR, is often an iterative process. You might:
    1.  Build an initial tree and perform ASR.
    2.  Identify potential issues (e.g., outlier sequences/dates from TreeTime's clock analysis, poorly supported ancestral states).
    3.  Refine your dataset (e.g., remove problematic sequences, correct metadata) or tree inference parameters.
    4.  Re-run the analysis.
    Workflow management systems (see Further Reading) can help manage these iterations.

## Conclusion

Ancestral sequence reconstruction is a powerful lens for viewing evolution. This tutorial has outlined a hybrid approach, leveraging versatile Python libraries like BioPython, DendroPy, and ETE Toolkit for data management and programmatic analysis, alongside specialized tools such as MAFFT for alignment, IQ-TREE for tree building, and TreeTime or DendroPy's own functions for ASR.

Successfully navigating large-scale viral datasets requires thoughtful strategies for data subsampling (e.g., with MMseqs2/CD-HIT), efficient tree inference, and potentially divide-and-conquer approaches for ASR. Whether using command-line tools like TreeTime for integrated time-scaling and ASR, or delving into Python-based ASR with DendroPy for finer control, the key is to choose tools and methods appropriate for your research question, dataset size, and computational resources. Remember to critically assess the uncertainty inherent in ASR and to ground your interpretations in robust evolutionary models and high-quality data.

## Further Reading/Resources

*   **Python Libraries:**
    *   **BioPython:** [https://biopython.org/](https://biopython.org/) (Tutorial: [https://biopython.org/DIST/docs/tutorial/Tutorial.html](https://biopython.org/DIST/docs/tutorial/Tutorial.html))
    *   **DendroPy:** [https://dendropy.org/](https://dendropy.org/) (Documentation: [https://dendropy.org/library/](https://dendropy.org/library/))
    *   **ETE Toolkit:** [http://etetoolkit.org/](http://etetoolkit.org/) (Tutorials: [http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html](http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html))
*   **ASR & Phylodynamics Tools:**
    *   **TreeTime:** [https://treetime.readthedocs.io/](https://treetime.readthedocs.io/)
    *   Sagulenko et al. (2018), TreeTime: Maximum-likelihood phylodynamic analysis. Virus Evolution.
*   **Sequence Alignment & Clustering:**
    *   **MAFFT:** [https://mafft.cbrc.jp/alignment/software/](https://mafft.cbrc.jp/alignment/software/)
    *   **MMseqs2:** [https://github.com/soedinglab/MMseqs2](https://github.com/soedinglab/MMseqs2)
    *   **CD-HIT:** [http://weizhongli-lab.org/cd-hit/](http://weizhongli-lab.org/cd-hit/)
*   **Phylogenetic Inference:**
    *   **IQ-TREE:** [http://www.iqtree.org/](http://www.iqtree.org/)
    *   **FastTree:** [http://www.microbesonline.org/fasttree/](http://www.microbesonline.org/fasttree/)
    *   **UShER:** [https://usher-wiki.readthedocs.io/](https://usher-wiki.readthedocs.io/)
*   **Tree Visualization:**
    *   **FigTree:** [http://tree.bio.ed.ac.uk/software/figtree/](http://tree.bio.ed.ac.uk/software/figtree/)
    *   **iTOL (Interactive Tree Of Life):** [https://itol.embl.de/](https://itol.embl.de/)
*   **Workflow Management:**
    *   **Nextflow:** [https://www.nextflow.io/](https://www.nextflow.io/)
    *   **Snakemake:** [https://snakemake.readthedocs.io/](https://snakemake.readthedocs.io/)

This tutorial provides a foundation. The field of phylogenetics is dynamic; continued exploration of documentation and new methods is encouraged to enhance your evolutionary analyses.

---

## Credits
This work is made possible by the collaborative efforts of Jules and Gemini 2.5 Pro (Google).