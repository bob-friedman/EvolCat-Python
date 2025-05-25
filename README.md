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
6.  [Workflow Examples](#workflow-examples)
    *   [A. Building a Local Sequence Database from NCBI](#a-building-a-local-sequence-database-from-ncbi)
        *   [Step 1: Identifying and Retrieving Initial Sequences](#step-1-identifying-and-retrieving-initial-sequences)
        *   [Step 2: Fetching Full Records for Retrieved IDs](#step-2-fetching-full-records-for-retrieved-ids)
        *   [Step 3: Converting Formats](#step-3-converting-formats)
        *   [Step 4: Extracting Relevant Information](#step-4-extracting-relevant-information)
        *   [Step 5: Cleaning and Standardizing Sequence Data](#step-5-cleaning-and-standardizing-sequence-data)
        *   [Step 6: Merging and Organizing Your Local Database](#step-6-merging-and-organizing-your-local-database)
    *   [B. Performing a Phylogenetic Tree Analysis](#b-performing-a-phylogenetic-tree-analysis)
        *   [Step 1: Sequence Preparation and Alignment](#step-1-sequence-preparation-and-alignment)
        *   [Step 2: Alignment Curation](#step-2-alignment-curation)
        *   [Step 3: Phylogenetic Inference (External Tools)](#step-3-phylogenetic-inference-external-tools)
        *   [Step 4: Tree Visualization and Basic Manipulation (Conceptual)](#step-4-tree-visualization-and-basic-manipulation-conceptual)
    *   [C. Guide to Accessing MHC Sequence Databases](#c-guide-to-accessing-mhc-sequence-databases)
    *   [D. Guide to Interpreting Phylogenetic Trees with Python](#d-guide-to-interpreting-phylogenetic-trees-with-python)
    *   [E. Special Topic: Virus Genomics, Diversity, and Analysis](#e-special-topic-virus-genomics-diversity-and-analysis)
7.  [Testing](#testing)
8.  [Development and Contributions](#development-and-contributions)

## Overview

The library is organized into:

*   `pylib/utils/`: Contains core utility modules for tasks like sequence parsing.
*   `pylib/scripts/`: Contains executable Python scripts that replicate and extend the functionality of original bioinformatics command-line tools. Many of these scripts depend on the `pylib/utils/` core utility modules. The scripts are designed to find these modules by default when EvolCat-Python is structured with `pylib/utils/` as a subdirectory.
    *   `pylib/scripts/ncbi/`: Contains tools specifically for interacting with NCBI.
    *   `pylib/scripts/viral_tools/`: Contains tools specifically for viral genomics analysis.

## Dependencies

The primary dependencies for this library are listed in the `requirements.txt` file and `setup.py`. Key dependencies include:

*   Python 3.7 or higher (Python 3.10.0 specifically tested with WSL setup below)
*   Biopython
*   Matplotlib
*   Requests (specifically for the NCBI tools module)
*   (Other specific dependencies might be required by individual scripts or future additions. External phylogenetic software like ClustalW/MUSCLE/MAFFT, RAxML/IQ-TREE/PhyML will be needed for parts of the phylogenetic workflow).

Please refer to `requirements.txt` or `setup.py` for a complete list of dependencies and their versions.

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

## NCBI Tools

This library includes specific tools for interacting with NCBI services. These are located in `pylib/scripts/ncbi/` and include:

*   `blast_ncbi_seq_def_py` (likely `pylib/scripts/ncbi/blast_ncbi_seq_def.py`): Retrieves sequence definitions from NCBI based on BLAST results.
*   `query_ncbi_gi_py` (likely `pylib/scripts/ncbi/query_ncbi_gi.py`): Queries NCBI using GenInfo Identifiers (GIs) to retrieve sequence data.

These tools require the `requests` library, which will be installed if you use `pip install .`. For more detailed information on these specific NCBI tools, please see `pylib/scripts/ncbi/README.md` (if present) or use the `--help` flag with the scripts themselves.

## Workflow Examples

This section outlines general workflows, highlighting where EvolCat-Python scripts can be utilized.

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

#### Step 2: Alignment Curation

Raw alignments often need refinement.

*   **EvolCat-Python `nogaps.py`:** Removes columns containing mostly gaps or ambiguous characters from your aligned FASTA file. This can improve the quality of phylogenetic inference.
    ```bash
    python3 pylib/scripts/nogaps.py aligned_sequences.afa > aligned_sequences_nogaps.afa
    ```
*   **Manual Inspection/Other Tools:** You might also use alignment viewers (e.g., Jalview, AliView) to manually inspect and edit the alignment.

#### Step 3: Phylogenetic Inference (External Tools)

Once you have a curated alignment, you can infer a phylogenetic tree. EvolCat-Python currently focuses on preparing data for these tools and can help convert formats.

*   **Format Conversion for Tree Builders:** Many tree-building programs require specific input formats, like PHYLIP.
    *   **EvolCat-Python `fas2phy.py`:** Converts your curated aligned FASTA file to PHYLIP sequential format. **Important:** All sequences must be of the same length.
        ```bash
        python3 pylib/scripts/fas2phy.py aligned_sequences_nogaps.afa > aligned_sequences.phy
        ```
*   **Tree Building (External Tools):** Use specialized software for tree inference. Common choices include:
    *   **Distance-based:** Phylip's `neighbor` program (Felsenstein 1989) (can use distance matrices, which EvolCat-Python can help calculate - see below).
    *   **Maximum Likelihood (ML):** RAxML, IQ-TREE, PhyML (Yang 1997 for PAML which includes ML methods; Adachi and Hasegawa 1996 for Molphy).
    *   **Bayesian Inference (BI):** MrBayes, BEAST.
    These programs will take your PHYLIP (or FASTA) alignment and produce a tree file, typically in Newick format (e.g., `my_tree.nwk`).
    ```bash
    # Example using an external ML tool (conceptual)
    # raxml-ng --msa aligned_sequences.phy --model GTR+G --prefix my_analysis --threads 2
    # iqtree -s aligned_sequences.phy -m GTR+G -bb 1000
    ```
*   **Calculating Distance Matrices (Optional, for distance-based tree methods):**
    If you plan to use distance-based tree building methods (like Neighbor-Joining), you first need a distance matrix.
    *   **EvolCat-Python `calculate_dna_distances.py` or `calculate_k2p.py`:** These scripts can take your aligned FASTA file and produce a table of pairwise distances. This table would then need to be formatted appropriately for input into a distance-based tree program (e.g., Phylip `neighbor`).
        ```bash
        python3 pylib/scripts/calculate_dna_distances.py aligned_sequences_nogaps.afa > dna_distances.tsv
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

### C. Guide to Accessing MHC Sequence Databases

[MHC Database Guide](https://github.com/bob-friedman/EvolCat-Python/blob/main/mhc-database-guide.md)

### D. Guide to Interpreting Phylogenetic Trees with Python

[Phylogenetic Tree Interpretation](https://github.com/bob-friedman/EvolCat-Python/blob/main/phylogenetic-tree-interpretation.md)

### E. Special Topic: Virus Genomics, Diversity, and Analysis

[Guide to Virus Genomics, Diversity, and Analysis](virus_genomics_guide.md)

## Testing

This project uses Python's built-in `unittest` framework for testing. Unit tests are located alongside the scripts they test within the `pylib/scripts/` directory (e.g., tests for NCBI scripts are in `pylib/scripts/ncbi/`).

To discover and run all tests from the root directory of the project, use the following command:

```bash
python -m unittest discover -s pylib/scripts/ -p "test_*.py"
```

This command will find and execute all files named `test_*.py` within the `pylib/scripts/` directory and its subdirectories (like `pylib/scripts/ncbi/`), ensuring that all components of the library are checked.

For more specific testing, such as running tests for a particular module or script, unit tests can target specific files or directories. For example, to run only the NCBI-related tests, use the following command:
```bash
python -m unittest discover pylib/scripts/ncbi/
```

## Development and Contributions

This library was primarily converted from its original Perl source using AI-assisted tooling (Model: Gemini 2.5 Pro, Agent: Jules AI). The process involved:

*   Analyzing original Perl scripts.
*   Translating Perl logic to Python, utilizing standard libraries like Biopython and Matplotlib.
*   Structuring the code into a Python package with scripts and utility modules.
*   Creating command-line interfaces using `argparse`.
*   Setting up packaging with `setup.py`.
*   Developing initial documentation.

Human oversight and review are crucial for ensuring the accuracy and robustness of the converted code. This library is currently under development. Contributions, bug reports, and feature requests are welcome via GitHub Issues and Pull Requests.
