<a name="top"></a>
# EvolCat-Python: A Python Suite for Evolutionary and Comparative Genomics

EvolCat-Python is a collection of Python-based scripts designed to facilitate common tasks in bioinformatics, particularly in the fields of evolutionary biology and comparative genomics. This project is a Python-based bioinformatics library converted from an existing collection of Perl scripts. It provides a set of tools for sequence manipulation, format conversion, analysis, and interaction with NCBI databases.

This README provides an overview, installation instructions, and conceptual guides on how you might use EvolCat-Python scripts.

**User Responsibility:**
The scripts and library components provided here are for research and informational purposes. Users are responsible for validating the results obtained using this software, interpreting them correctly, and ensuring they are appropriate for their specific application. The original authors and the converters of this code disclaim any liability for its use or misuse. It is recommended to test the tools with known datasets and compare results with other established bioinformatics software where appropriate. Users may need to adapt or modify the code to suit their specific research needs and computational environment.

## Table of Contents

1.  [Overview](#overview)
2.  [Dependencies](#dependencies)
3.  [Installation](#installation)
    *   [General Python Environment Setup](#general-python-environment-setup)
    *   [Installing EvolCat-Python](#installing-evolcat-python)
    *   [Access in a Windows OS with WSL](#access-in-a-windows-os-with-wsl-windows-subsystem-for-linux)
4.  [General Script Usage](#general-script-usage)
5.  [NCBI Tools](#ncbi-tools)
6.  [Relationship with Biopython and Scope of Provided Scripts](#relationship-with-biopython-and-scope-of-provided-scripts)
7.  [Workflow Examples](#workflow-examples)
    *   [A. Building a Local Sequence Database from NCBI](#a-building-a-local-sequence-database-from-ncbi)
    *   [B. Performing a Phylogenetic Tree Analysis](#b-performing-a-phylogenetic-tree-analysis)
    *   [C. Performing a Standard Pairwise Sequence Alignment with Biopython](#c-performing-a-standard-pairwise-sequence-alignment-with-biopython)
    *   [D. Basic Motif Scanning](#d-basic-motif-scanning)
    *   [H. VCF File Analysis and Filtering](#h-vcf-file-analysis-and-filtering)
    *   [I. Extracting and Analyzing CDS Regions](#i-extracting-and-analyzing-cds-regions)
    *   [J. SARS-CoV-2 Lineage Classification Pipeline](#j-sars-cov-2-lineage-classification-pipeline)
8.  [Virus Genomics, Diversity, and Analysis](#virus-genomics-diversity-and-analysis)
    *   [Special Topic: Virus Genomics, Diversity, and Analysis](#g-special-topic-virus-genomics-diversity-and-analysis)
9.  [Technical Guides](#technical-guides)
    *   [Guide to Accessing MHC Sequence Databases](#e-guide-to-accessing-mhc-sequence-databases)
    *   [Guide to Interpreting Phylogenetic Trees with Python](#f-guide-to-interpreting-phylogenetic-trees-with-python)
10. [Detailed Script Usage](#detailed-script-usage)
11. [Development and Contributions](#development-and-contributions)
12. [Citation](#citation)


[Back to Top](#top)

## Overview

The library is organized into:

*   `pylib/utils/`: Contains core utility modules for tasks like sequence parsing.
*   `pylib/scripts/`: Contains executable Python scripts that replicate and extend the functionality of original bioinformatics command-line tools. Many of these scripts depend on the `pylib/utils/` core utility modules. The scripts are designed to find these modules by default when EvolCat-Python is structured with `pylib/utils/` as a subdirectory. Some scripts may have their own detailed `README.md` files within this directory (e.g., `extract_cds_region.py`).
    *   `pylib/scripts/ncbi/`: Contains tools specifically for interacting with NCBI.
    *   `pylib/scripts/paml_tools/`: Contains tools specifically for PAML genomics analysis.


[Back to Top](#top)

## Dependencies

The primary dependencies for this library are listed in the `requirements.txt` file and `setup.py`. Key dependencies include:

*   Python 3.7 or higher (Python 3.10.0 specifically tested with WSL setup below)
*   Biopython
*   Matplotlib
*   Requests (specifically for the NCBI tools module)
*   (Other specific dependencies might be required by individual scripts or future additions. External phylogenetic software like ClustalW/MUSCLE/MAFFT, RAxML/IQ-TREE/PhyML will be needed for parts of the phylogenetic workflow).

Please refer to `requirements.txt` or `setup.py` for a complete list of dependencies and their versions.


[Back to Top](#top)

## Installation

### General Python Environment Setup
It is highly recommended to use a virtual environment (e.g., `venv` or `conda`) to manage dependencies for this project.

```bash
# Using venv (standard Python)
python3 -m venv evolcat_env
source evolcat_env/bin/activate # On Unix/macOS
# evolcat_env\Scripts\activate # On Windows

# Or using conda
# conda create -n evolcat_env python=3.10
# conda activate evolcat_env
```
Once your environment is activated, you can install dependencies.


[Back to Top](#top)

### Installing EvolCat-Python

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/bob-friedman/EvolCat-Python.git
    cd EvolCat-Python
    ```
2.  **Install the package and its dependencies:**
    If you are installing the package from the root directory of this repository (i.e., after cloning and `cd EvolCat-Python`), you can install it along with all dependencies by running:
    ```bash
    pip install .
    ```
    This command reads `setup.py` and `requirements.txt` to install EvolCat-Python and its listed dependencies like Biopython, Matplotlib, and Requests.

    *(Note: If the package were published to PyPI, e.g., as `evolcat_python`, you would install it via `pip install evolcat_python`.)*


[Back to Top](#top)

### Access in a Windows OS with WSL (Windows Subsystem for Linux)

To access these scripts from a Linux environment in supported versions of Windows, the first step is to verify the installation of WSL. The subsequent steps involve using "pip" to manage packages in Python. However, the default Python version shipped with some Linux distributions (like Ubuntu) available via WSL might have compatibility issues with `pip`'s ability to install packages in user-writable locations or may lead to conflicts with system-managed Python packages.

To work around such issues and manage Python versions more effectively within WSL, `pyenv` can be used. The following steps describe a third-party procedure that has been tested on one system. **These steps modify your shell environment and install software; proceed with caution and understand what each command does.**

In your Ubuntu/WSL shell:

1.  **Install build dependencies for `pyenv` and Python:**
    ```bash
    sudo apt update
    sudo apt install -y gcc make build-essential libssl-dev libffi-dev zlib1g-dev libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm libncurses5-dev xz-utils tk-dev libxml2-dev libxmlsec1-dev liblzma-dev
    ```
2.  **Install `pyenv`:**
    ```bash
    curl https://pyenv.run | bash
    ```
3.  **Configure your shell environment for `pyenv`:**
    Add the following lines to your `~/.bashrc` file (or `~/.zshrc` if using Zsh):
    ```bash
    export PYENV_ROOT="$HOME/.pyenv"
    command -v pyenv >/dev/null || export PATH="$PYENV_ROOT/bin:$PATH"
    eval "$(pyenv init -)"
    ```
    Then, apply the changes by sourcing your shell configuration file (e.g., `source ~/.bashrc`) or by opening a new terminal.

4.  **Install a specific Python version using `pyenv`:**
    (e.g., Python 3.10.0, which is relatively recent as of May 2025 context from original notes)
    ```bash
    pyenv install 3.10.0
    ```
5.  **Set the Python version for your environment:**
    You can set it globally or locally for the current project directory.
    ```bash
    pyenv global 3.10.0  # Set 3.10.0 as the default version of Python for your user
    # OR
    # pyenv local 3.10.0   # Set 3.10.0 as the version when running within the current folder
    ```
    After setting the Python version with `pyenv`, `python3` and `pip` commands will use this version. You can then proceed with cloning the EvolCat-Python repository and installing it using `pip install .` as described above within your WSL environment.


[Back to Top](#top)

## General Script Usage

*   All executable scripts are located in the `pylib/scripts/` directory (and subdirectories like `pylib/scripts/ncbi/`).
*   If you have installed the package using `pip install .`, some scripts might be available directly on your PATH (depending on `setup.py` configuration). Otherwise, run scripts from the command line using `python3 path/to/script.py [arguments]`.
    For example, from the root `EvolCat-Python` directory:
    ```bash
    python3 pylib/scripts/gb2fasta.py my_input.gb
    ```
*   Use the `-h` or `--help` flag with any script to see its specific command-line options and a brief description.
    ```bash
    python3 pylib/scripts/gb2fasta.py -h
    ```


[Back to Top](#top)

## NCBI Tools

This library includes specific tools for interacting with NCBI services. These are located in `pylib/scripts/ncbi/` and include:

*   `blast_ncbi_seq_def_py` (likely `pylib/scripts/ncbi/blast_ncbi_seq_def.py`): Retrieves sequence definitions from NCBI based on BLAST results.
*   `query_ncbi_gi_py` (likely `pylib/scripts/ncbi/query_ncbi_gi.py`): Queries NCBI using GenInfo Identifiers (GIs) to retrieve sequence data.

These tools require the `requests` library, which will be installed if you use `pip install .`. For more detailed information on these specific NCBI tools, please see `pylib/scripts/ncbi/README.md` (if present) or use the `--help` flag with the scripts themselves.


[Back to Top](#top)

## Relationship with Biopython and Scope of Provided Scripts

EvolCat-Python is built upon and requires [Biopython](https://biopython.org) as a core dependency. Many scripts in this suite act as convenient command-line wrappers or implement common workflows by intelligently utilizing Biopython's underlying capabilities. This approach allows for rapid execution of specific, often complex, tasks directly from the command line, making them powerful tools for bioinformatics analyses. The [Workflow Examples](#workflow-examples) section of this README further illustrates how these scripts can be combined to streamline multi-step processes.

While EvolCat-Python offers these targeted solutions, Biopython itself is a comprehensive library with vast functionalities. For highly advanced, programmatic, or uniquely customized analyses, users are encouraged to leverage Biopython's modules directly. The following areas highlight how EvolCat-Python's scripts provide useful functionalities, and also where Biopython offers significantly more depth for those wishing to delve deeper:

*   **Performing General Pairwise Alignments:**
    *   While EvolCat-Python includes `approximate_string_match.py` for calculating edit distance (a measure of sequence dissimilarity), Biopython's `Bio.Align.PairwiseAligner` provides robust tools for standard biological global (Needleman-Wunsch) and local (Smith-Waterman) alignments with full control over substitution matrices (e.g., BLOSUM, PAM) and affine gap penalties.

*   **Advanced Multiple Sequence Alignment (MSA) Objects and Analysis:**
    *   EvolCat-Python offers `nogaps.py` for basic MSA curation and `analyze_msa.py` for obtaining alignment statistics and consensus sequences.
    *   For more intricate MSA manipulations, Biopython's `Bio.Align` module (and the older `Bio.AlignIO`) provides rich `Alignment` objects for parsing numerous MSA formats, deriving detailed conservation scores, and performing complex operations directly on MSA objects.

*   **Phylogenetic Analysis:**
    *   EvolCat-Python provides a suite of tools to facilitate phylogenetic analysis. Scripts like `calculate_dna_distances.py` and `calculate_k2p.py` generate distance matrices, and `fas2phy.py` aids in format conversion for external tree-building programs. Furthermore, `build_tree_from_distances.py` allows for direct construction of phylogenetic trees using NJ or UPGMA methods from a distance matrix, offering a convenient way to generate initial tree topologies.
    *   For more comprehensive phylogenetic tasks, Biopython's `Bio.Phylo` module offers extensive capabilities, including parsing and manipulating various tree formats (Newick, Nexus, NeXML), implementing a wider array of tree construction algorithms, advanced tree operations (rerooting, pruning, comparing topologies), and sophisticated tree visualization.

*   **Sequence Analysis and Statistics:**
    *   EvolCat-Python includes scripts for specific analyses such as `calculate_nucleotide_diversity.py`, which directly computes population genetics statistics from sequence data. This provides a ready-to-use tool for assessing genetic variation.
    *   While such specific metrics can be computed using Biopython's foundational objects with custom scripting, Biopython's strength lies in providing the building blocks for a broader range of sequence analyses.

*   **Sequence Motif Discovery and Analysis:**
    *   EvolCat-Python includes `scan_sequences_for_motif.py` for searching for user-defined motifs (including IUPAC strings) within sequences, a valuable tool for identifying potential functional sites. It also offers `find_tandem_repeats.py` for a specific type of repeat.
    *   Biopython's `Bio.motifs` module provides a more extensive toolkit for creating, representing (e.g., position-weight matrices - PWMs), and scanning for sequence motifs, as well as interfacing with motif databases and external tools like MEME.

*   **Evolutionary Rate (dN/dS) Estimation:**
    *   EvolCat-Python provides excellent script-based wrappers for PAML's `yn00` and `codeml` programs (via `calculate_dn_ds.py` and `calculate_site_specific_ds_dn.py`), simplifying the process of estimating dN/dS ratios. These are powerful additions for studying molecular evolution.
    *   For users who prefer or need to implement alternative methods within Python, Biopython's `Bio.Align.analysis` module allows for the direct calculation of dN/dS ratios using several established algorithms (e.g., Nei-Gojobori, Li-Wu-Luo, Yang-Nielsen).

*   **Flexible Parsing of Diverse Sequence Search Outputs (`Bio.SearchIO`):**
    *   EvolCat-Python scripts such as `parse_blast_text.py` and `blast_to_table.py` are useful for processing standard text BLAST output.
    *   For broader compatibility, Biopython's `Bio.SearchIO` module provides a modern, unified interface for parsing various output formats (including XML and tabular) from a wide range of sequence search tools (BLAST, HMMER, Exonerate, etc.), offering a more robust and extensible approach for diverse search result inputs.

Users familiar with Python programming can readily combine the convenience and focused power of EvolCat-Python scripts with the extensive library functionalities of Biopython to build sophisticated and efficient bioinformatics pipelines.


[Back to Top](#top)

## Workflow Examples

This section outlines general workflows, highlighting where EvolCat-Python scripts can be utilized.


[Back to Top](#top)

### A. Building a Local Sequence Database from NCBI

This workflow describes creating a local sequence database using NCBI resources.

#### Step 1: Identifying and Retrieving Initial Sequences
Often, building a local database starts with a set of query sequences or keywords to find related entries in NCBI.

*   **Using NCBI Website or E-utilities:** Perform initial searches directly on the NCBI website (e.g., BLAST, Entrez search) or use E-utilities (e.g., `esearch`, `efetch`) via command-line tools like `efetch` from NCBI's Entrez Direct, or by scripting with Biopython's `Bio.Entrez` module.
*   **EvolCat-Python for BLAST output processing:** If you run a BLAST search (Altschul et al. 1997) and get text output:
    *   `pylib/scripts/parse_blast_text.py`: Parses standard text BLAST output.
    *   `pylib/scripts/blast_to_table.py`: Converts BLAST text output to a tab-delimited table.
    *   `pylib/scripts/find_blast_top_pairs.py`: Helps identify top-scoring unique pairs.
    From these results, extract accession numbers for your local database.

#### Step 2: Fetching Full Records for Retrieved IDs
Once you have a list of NCBI accession numbers, download their full records.

*   **Using NCBI E-utilities (`efetch`):** Use NCBI's Entrez Direct `efetch` or Biopython's `Bio.Entrez.efetch` to download records (e.g., GenBank, FASTA).
    ```bash
    # Example using Entrez Direct (assuming you have a file of accessions: accs.txt)
    # efetch -db nucleotide -format gb -input accs.txt > my_sequences.gb
    # efetch -db protein -format fasta -input accs.txt > my_proteins.fasta
    ```

#### Step 3: Converting Formats
NCBI often provides GenBank files; FASTA is often preferred for downstream tasks.

*   **EvolCat-Python `gb2fasta.py`:** Converts GenBank files to FASTA.
    ```bash
    python3 pylib/scripts/gb2fasta.py my_sequences.gb > my_sequences.fasta
    ```
*   Other EvolCat-Python conversion scripts like `phy2fas.py` can be used if your source data is in other formats.

#### Step 4: Extracting Relevant Information
You might only need specific features, like Coding DNA Sequences (CDS).

*   **EvolCat-Python `gbCDS.py`:** Extracts CDS information (Locus ID, nucleotide sequence, translation) from GenBank files.
    ```bash
    python3 pylib/scripts/gbCDS.py my_genome.gb > my_cds_info.txt
    # The output of gbCDS.py might need further parsing to create a FASTA file of just CDS sequences.
    ```
*   **EvolCat-Python `translate_seq.py`:**
    If you have nucleotide FASTA files (e.g., extracted CDS), you can translate them to protein sequences.
    ```bash
    python3 pylib/scripts/translate_seq.py --input_file my_cds.fna > my_proteins.faa
    ```

#### Step 5: Cleaning and Standardizing Sequence Data
Sequence headers from public databases can be long and contain characters problematic for some software.

*   **EvolCat-Python `clean_fasta_name.py`:**
    Standardizes FASTA headers by replacing spaces and periods with underscores and converting to uppercase.
    ```bash
    python3 pylib/scripts/clean_fasta_name.py my_sequences.fasta > cleaned_sequences.fasta
    ```
*   **EvolCat-Python `nogaps.py`:**
    If your initial dataset comes from alignments, this script removes columns containing gaps or ambiguous characters.
    ```bash
    python3 pylib/scripts/nogaps.py aligned_sequences.fasta > aligned_nogaps.fasta
    ```

#### Step 6: Merging and Organizing Your Local Database
If you have retrieved sequences from multiple searches or in multiple files, you'll want to combine them.

*   **EvolCat-Python `merge_fastas.py`:**
    Merges multiple FASTA files. If sequences have identical headers, their sequence parts are concatenated.
    ```bash
    python3 pylib/scripts/merge_fastas.py part1.fasta part2.fasta part3.fasta > local_database.fasta
    ```
    After merging, you will have a single FASTA file representing your initial local sequence database. You might want to organize sequences into sub-databases based on organism, gene family, etc., by further splitting this master file if needed (which could be done with custom scripting or standard command-line tools). You can then format this `local_database.fasta` using NCBI's `makeblastdb` (external tool) if you wish to BLAST against it locally.


[Back to Top](#top)

### B. Performing a Phylogenetic Tree Analysis

This workflow outlines the steps involved in phylogenetic analysis, indicating where EvolCat-Python scripts and external tools are typically used. The `EvolTree.py` module mentioned in earlier design documents is conceptual at this stage for direct tree building within EvolCat-Python; the current scripts primarily support pre- and post-processing for external tree builders.

#### Step 1: Sequence Preparation and Alignment

Phylogenetic analysis begins with a set of homologous sequences, typically in FASTA format. These might come from your [local database](#a-building-a-local-sequence-database-from-ncbi) or other sources.

*   **Sequence Retrieval & Formatting:** Use scripts like `gb2fasta.py`, `clean_fasta_name.py`, or `merge_fastas.py` as needed to prepare your input FASTA file.
*   **Multiple Sequence Alignment (MSA) (External Tools):** This is a critical step. EvolCat-Python does not perform MSA itself but prepares data for, and can process outputs from, common MSA programs. You will need to use external software like:
    *   Clustal Omega / ClustalW (Thompson et al. 1994)
    *   MUSCLE
    *   MAFFT
    These programs will take your multi-FASTA file as input and produce an aligned FASTA file (or other alignment formats).
    ```bash
    # Example using an external MSA tool (conceptual)
    # muscle -in unaligned_sequences.fasta -out aligned_sequences.afa
    ```

#### Step 2: Alignment Curation and Basic Analysis

Raw alignments often need refinement.

*   **EvolCat-Python `nogaps.py`:** Removes columns containing mostly gaps or ambiguous characters from your aligned FASTA file. This can improve the quality of phylogenetic inference.
    ```bash
    python3 pylib/scripts/nogaps.py aligned_sequences.afa > aligned_sequences_nogaps.afa
    ```
*   **EvolCat-Python `analyze_msa.py`:** After generating or obtaining an MSA, this script can be used to:
    *   Calculate and view a consensus sequence.
    *   Get basic statistics about the alignment (length, number of sequences, gap information).
    *   Convert the MSA to other formats if needed for downstream tools.
    ```bash
    # Example: Get stats and consensus for an alignment in Clustal format
    python3 pylib/scripts/analyze_msa.py my_alignment.aln --informat clustal --get_stats --get_consensus
    ```
*   **Manual Inspection/Other Tools:** You might also use alignment viewers (e.g., Jalview, AliView) to manually inspect and edit the alignment.

#### Step 3: Phylogenetic Inference (External Tools)

Once you have a curated alignment, you can infer a phylogenetic tree. EvolCat-Python currently focuses on preparing data for these tools and can help convert formats.

*   **Format Conversion for Tree Builders:** Many tree-building programs require specific input formats, like PHYLIP.
    *   **EvolCat-Python `fas2phy.py`:** Converts your curated aligned FASTA file to PHYLIP sequential format. **Important:** All sequences must be of the same length.
        ```bash
        python3 pylib/scripts/fas2phy.py aligned_sequences_nogaps.afa > aligned_sequences.phy
        ```
*   **Calculating Distance Matrices (Optional, for distance-based tree methods):**
    If you plan to use distance-based tree building methods (like Neighbor-Joining), you first need a distance matrix.
    *   **EvolCat-Python `calculate_dna_distances.py` or `calculate_k2p.py`:** These scripts can take your aligned FASTA file and produce a table of pairwise distances. This table would then need to be formatted appropriately for input into a distance-based tree program (e.g., Phylip `neighbor`).
        ```bash
        python3 pylib/scripts/calculate_dna_distances.py aligned_sequences_nogaps.afa > dna_distances.tsv
        ```
    *   **Building a Tree from Distances (EvolCat-Python):**
        *   Once you have a distance matrix (e.g., in PHYLIP format, potentially generated from the scripts above and formatted), you can use **EvolCat-Python's `pylib/scripts/build_tree_from_distances.py`** script to construct a tree using Neighbor-Joining (NJ) or UPGMA methods. This script outputs the tree in Newick format.
            ```bash
            # Example assuming my_distances.phy is your PHYLIP distance matrix
            python3 pylib/scripts/build_tree_from_distances.py my_distances.phy --method nj --outfile my_tree.nwk
            ```
            This provides a direct way to generate a tree topology from your calculated distances within the EvolCat-Python suite.
*   **Tree Building (Alternative External Tools for Advanced Methods):**
    For more advanced phylogenetic models or methods (e.g., Maximum Likelihood, Bayesian Inference), or if you prefer other specialized distance-based programs, you can use external software. Common choices include:
    *   **Distance-based (alternative):** Phylip's `neighbor` program (Felsenstein 1989).
    *   **Maximum Likelihood (ML):** RAxML, IQ-TREE, PhyML.
    *   **Bayesian Inference (BI):** MrBayes, BEAST.
    These programs will typically take your curated alignment (e.g., PHYLIP or FASTA) or a distance matrix and produce a tree file (e.g., `my_tree.nwk`).
    ```bash
    # Example using an external ML tool (conceptual)
    # raxml-ng --msa aligned_sequences.phy --model GTR+G --prefix my_analysis --threads 2
    # iqtree -s aligned_sequences.phy -m GTR+G -bb 1000
    ```

#### Step 4: Tree Visualization and Basic Manipulation (Conceptual)

Once you have a tree file (e.g., in Newick format), you can visualize and analyze it.

*   **Visualization (External Tools or Libraries):**
    *   Dedicated tree viewers: FigTree, iTOL, Dendroscope.
    *   Python libraries: Biopython's `Bio.Phylo` can parse and draw trees (often with Matplotlib for better graphics). The `EvolPlot.py` script in EvolCat-Python might incorporate some of these capabilities or provide wrappers.
        ```python
        # Example using Biopython directly (conceptual, not an EvolCat script)
        # from Bio import Phylo
        # import matplotlib.pyplot as plt
        # tree = Phylo.read("my_tree.nwk", "newick")
        # Phylo.draw(tree)
        # plt.show()
        ```
*   **Basic Tree Manipulation (Biopython):** Biopython's `Bio.Phylo` module allows for programmatic tree manipulation (e.g., re-rooting, extracting clades, calculating distances between tips). While `EvolCat-Python` might not have dedicated high-level scripts for all these manipulations yet, the underlying Biopython library provides these functions.

This phylogenetic workflow highlights how EvolCat-Python scripts primarily serve as helper tools for preparing data for, and converting formats between, specialized external programs for MSA and tree inference.


[Back to Top](#top)

### C. Performing a Standard Pairwise Sequence Alignment with Biopython

A fundamental task in bioinformatics is aligning two sequences to identify regions of similarity, which can imply functional, structural, or evolutionary relationships. While EvolCat-Python provides `approximate_string_match.py` for edit distance, it does not currently offer a standalone script for performing standard biological pairwise alignments (e.g., Needleman-Wunsch global or Smith-Waterman local alignments) using comprehensive scoring systems (e.g., BLOSUM/PAM matrices for proteins, or specific match/mismatch scores for DNA, along with affine gap penalties).

However, since EvolCat-Python depends on Biopython, you can easily perform these alignments by writing a short Python script that utilizes Biopython's `Bio.Align.PairwiseAligner`.

**Objective:** To find the optimal alignment(s) between two sequences.

**Steps:**

1.  **Prepare Input Sequences:**
    *   Verify that two sequences (such as `sequence1` and `sequence2`) are in separate FASTA files (e.g., `seq1.fasta`, `seq2.fasta`). Note that these files are also compatible with the EvolCat-Python scripts `gb2fasta.py` and `extract_region.py`.

2.  **Create a Python Script for Alignment:**
    *   Create a new Python file (e.g., `run_pairwise_alignment.py`).
    *   The script will use Biopython modules. Below is a conceptual outline:

    ```python
    from Bio import SeqIO
    from Bio import Align
    from Bio.Align import substitution_matrices # If using predefined matrices like BLOSUM/PAM

    # --- 1. Load your sequences ---
    try:
        seq_record1 = SeqIO.read("seq1.fasta", "fasta")
        seq_record2 = SeqIO.read("seq2.fasta", "fasta")
    except FileNotFoundError as e:
        print(f"Error: One of the sequence files was not found: {e}")
        exit()

    sequence1 = seq_record1.seq
    sequence2 = seq_record2.seq

    # --- 2. Initialize the PairwiseAligner ---
    aligner = Align.PairwiseAligner()

    # --- 3. Configure Aligner Parameters ---
    # Choose mode: 'global' (Needleman-Wunsch) or 'local' (Smith-Waterman)
    aligner.mode = 'global'  # Or 'local'

    # Example for Protein Alignment (using BLOSUM62):
    # aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    # aligner.open_gap_score = -11  # Typical BLOSUM62 gap open penalty
    # aligner.extend_gap_score = -1   # Typical BLOSUM62 gap extend penalty

    # Example for DNA Alignment:
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2.5 # Example gap open
    aligner.extend_gap_score = -1 # Example gap extend

    # For more nuanced gap penalties (e.g., different end gap scores):
    # aligner.target_end_gap_score = 0.0
    # aligner.query_end_gap_score = 0.0

    # --- 4. Perform the Alignment ---
    print(f"Aligning {seq_record1.id} and {seq_record2.id}...")
    alignments = aligner.align(sequence1, sequence2)

    # --- 5. Print and Interpret Results ---
    if not alignments:
        print("No alignments found.")
    else:
        print(f"Number of optimal alignments: {len(alignments)}")
        print(f"Alignment score: {alignments.score}") # All alignments in the iterator share the same score
        
        # Print the first optimal alignment (or loop through all alignments)
        # For very large numbers of alignments, you might only want to print one or a few.
        # Use list(alignments) with caution if len(alignments) is huge.
        print("First alignment:")
        print(alignments[0]) 
        
        # Example: To print all alignments if there are few:
        # for i, alignment in enumerate(alignments):
        #     print(f"Alignment {i+1}:")
        #     print(alignment)
        #     # You can also access coordinates, etc.
        #     # print(alignment.coordinates)

    print("Alignment complete.")
    ```

3.  **Run the Python Script:**
    ```bash
    python3 run_pairwise_alignment.py
    ```

4.  **Interpret Results:**
    *   The script will output the alignment score. Higher scores generally indicate better similarity.
    *   The alignment itself will be printed, showing how the sequences align, with dashes (`-`) indicating gaps.
    *   For local alignments, only the highest-scoring similar region(s) will be shown. For global alignments, the entire length of both sequences will be represented.

This conceptual workflow demonstrates how to leverage the power of Biopython for a core bioinformatics task that is not currently available as a dedicated script in the EvolCat-Python suite. By creating simple scripts like the one outlined, users can easily extend the capabilities of their analysis environment.


[Back to Top](#top)

### D. Basic Motif Scanning

Identifying known motifs (e.g., transcription factor binding sites, short functional patterns) within a set of sequences is a common requirement.

*   **Prepare Input Sequences:** Ensure your sequences are in a FASTA file (e.g., `promoter_regions.fasta`).
*   **Scan for Motif:** Use **EvolCat-Python's `pylib/scripts/scan_sequences_for_motif.py`** to search for occurrences of a known motif (IUPAC string or simple sequence) on both strands.
    ```bash
    # Example: Scan for the motif 'TGACGTCA' in promoter_regions.fasta
    python3 pylib/scripts/scan_sequences_for_motif.py promoter_regions.fasta --motif TGACGTCA --output_report found_motifs.tsv
    ```
    The script will output a tab-delimited file listing the sequence ID, motif, start, end, strand, and the matched sequence for each occurrence.


[Back to Top](#top)

## Technical Guides <a name="technical-guides"></a>
### Guide to Accessing MHC Sequence Databases <a name="e-guide-to-accessing-mhc-sequence-databases"></a>

[MHC Database Guide](guides/mhc-database-guide.md)


[Back to Top](#top)

### Guide to Interpreting Phylogenetic Trees with Python <a name="f-guide-to-interpreting-phylogenetic-trees-with-python"></a>

[Phylogenetic Tree Interpretation](guides/phylogenetic-tree-interpretation.md)


[Back to Top](#top)

## Virus Genomics, Diversity, and Analysis <a name="virus-genomics-diversity-and-analysis"></a>
### Special Topic: Virus Genomics, Diversity, and Analysis <a name="g-special-topic-virus-genomics-diversity-and-analysis"></a>

[Guide to Virus Genomics, Diversity, and Analysis](virology_tools/virus_genomics_guide.md)


[Back to Top](#top)

### H. VCF File Analysis and Filtering

Variant Call Format (VCF) files are essential in genomics for storing data about genetic variations. The process of identifying these genetic differences between a sample genome and a reference is known as variant calling. A common and crucial step after obtaining VCF files is to filter them to retain only high-quality and reliable variants for downstream analysis.

The `analyze_vcf.py` script, part of the EvolCat-Python suite, provides functionality for basic filtering of VCF files. This script allows users to filter variants based on metrics such as variant quality (QUAL) and read depth (DP).

For a more comprehensive understanding of VCF files, the variant calling process, its diverse research applications, and a detailed look at the `analyze_vcf.py` script's role and capabilities, please refer to the dedicated README in the `vcf_analysis_tools` directory. (See [Background on VCF Analysis](vcf_analysis_tools/README.md) for more details).

For specific command-line options and usage examples for the script, see [`analyze_vcf.py` usage](docs/USAGE.md#pylibscriptsanalyze_vcfpy).

### I. Extracting and Analyzing CDS Regions

The `pylib/scripts/extract_cds_region.py` script is a versatile tool for working with coding sequences (CDS) in GenBank files. It allows for two main types of operations:

1.  **Extracting a specific genomic region and translating overlapping CDS:** If you have a particular region of interest in a genome (defined by start and end coordinates), this script can extract the nucleotide sequence for that region. Furthermore, if any annotated CDS features overlap with your specified region, the script will identify the precise portion of the CDS that overlaps, extract its nucleotide sequence (correctly handling strand and exon information), and translate it into the corresponding protein sequence. This is useful for examining the protein products of specific genomic segments.

2.  **Analyzing a single nucleotide position to find its codon and amino acid:** If you are interested in a specific nucleotide position within a genome, this script can determine if that position falls within a CDS. If it does, the script will identify the codon that includes this nucleotide, indicate the nucleotide's position within the codon (1st, 2nd, or 3rd base), and provide the translated amino acid. This is helpful for detailed analysis of specific sites, such as investigating the impact of a single nucleotide variant (SNV) if it falls within a coding region.

The script handles complexities like forward and reverse strands, compound CDS locations (exons), and `codon_start` annotations to ensure accurate extraction and translation.

For detailed usage instructions, command-line options, input/output formats, and examples, please refer to its dedicated README file:
[extract_cds_region.py README](pylib/scripts/README.md)

**Example Scenarios:**

*   **Range Extraction:**
    ```bash
    python pylib/scripts/extract_cds_region.py my_genome.gb --start_pos 1000 --end_pos 1500
    ```
    This command would analyze the region from base 1000 to 1500 in `my_genome.gb`.

*   **Single Position Analysis:**
    ```bash
    python pylib/scripts/extract_cds_region.py my_genome.gb --single_position 1234
    ```
    This command would provide details about the codon and amino acid related to nucleotide position 1234 in `my_genome.gb`.

[Back to Top](#top)

### J. SARS-CoV-2 Lineage Classification Pipeline

[Documentation for the SARS-CoV-2 Lineage Classification Pipeline](pipelines/sars_cov2_lineage_classification/README.md)

[Back to Top](#top)

## Detailed Script Usage

For detailed command-line options and examples for each script, please refer to:
`docs/USAGE.md`

You can also use the `-h` or `--help` flag with any script:
```bash
python3 pylib/scripts/script_name.py -h
```


[Back to Top](#top)

## Development and Contributions

This library was primarily converted from its original Perl source using AI-assisted tooling (Model: Gemini 2.5 Pro, Agent: Jules AI). The process involved:

*   Analyzing original Perl scripts.
*   Translating Perl logic to Python, utilizing standard libraries like Biopython and Matplotlib.
*   Structuring the code into a Python package with scripts and utility modules.
*   Creating command-line interfaces using `argparse`.
*   Setting up packaging with `setup.py`.
*   Developing initial documentation.

Human oversight and review are crucial for ensuring the accuracy and robustness of the converted code. This library is currently under development. Contributions, bug reports, and feature requests are welcome via GitHub Issues and Pull Requests.


[Back to Top](#top)

## Citation

If you use EvolCat-Python in your research or find the repository useful, please cite the following paper:

Friedman, R. EvolCat-Python: A Python Suite for Evolutionary and Comparative Genomics. Preprints 2025, 2025052059. https://www.preprints.org/manuscript/202505.2059/v1


[Back to Top](#top)
