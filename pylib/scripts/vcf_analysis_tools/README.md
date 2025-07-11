# VCF Analysis Tools

This directory contains tools for analyzing VCF (Variant Call Format) files. The `analyze_vcf.py` script provides basic filtering capabilities for VCF data.

## Understanding VCF Files and Variant Calling

### What is a VCF file?

A **Variant Call Format (VCF)** file is a standardized text file format used in bioinformatics to store information about genetic sequence variations. It is the common output format from variant calling pipelines.

Key characteristics of VCF files:
*   **Text-based:** VCF files are human-readable and can be opened with standard text editors, although they are typically processed using specialized command-line tools or libraries.
*   **Structured Format:** Each VCF file consists of:
    *   **Meta-information lines (header):** These lines start with `##` and provide information about the VCF file itself, such as the file format version, reference genome used, and definitions of annotations used in the INFO, FILTER, and FORMAT fields.
    *   **Column header line:** A single line starting with `#CHROM` that defines the data columns.
    *   **Data lines:** Each subsequent line represents a variant call and contains information about a specific genetic variation.
*   **Standard Columns:** The main data columns include:
    *   `#CHROM`: The chromosome or contig where the variant is located.
    *   `POS`: The 1-based starting position of the variant on the chromosome.
    *   `ID`: An identifier for the variant (e.g., an rsID from dbSNP if available, otherwise often '.').
    *   `REF`: The reference base(s) at the given position.
    *   `ALT`: The alternate allele(s) observed in the sample(s). This can be a single base, multiple bases (for insertions/deletions), or a list of comma-separated alternate alleles.
    *   `QUAL`: A Phred-scaled quality score for the assertion made in ALT. Higher scores indicate higher confidence in the variant call.
    *   `FILTER`: A flag indicating whether the variant call passed quality filters (e.g., `PASS`) or failed specific filters.
    *   `INFO`: An extensible field containing additional information about the variant, often as key-value pairs (e.g., `DP` for read depth, `AF` for allele frequency).
    *   `FORMAT` (optional): If genotype information for samples is present, this column specifies the format of the sample data.
    *   Sample Data (optional): Subsequent columns contain genotype information for each sample included in the VCF file, corresponding to the `FORMAT` specification.

**Purpose:** The primary purpose of VCF files is to store and exchange information about genetic variations, including SNPs (Single Nucleotide Polymorphisms), indels (insertions and deletions), and structural variants.

### What is Variant Calling?

**Variant calling** is the process of identifying differences (variants) between the genome of an individual or sample and a reference genome. This process is fundamental to many areas of genomics research.

The general workflow for variant calling typically involves:
1.  **Sequencing:** Obtaining raw DNA sequence data from a sample using technologies like Next-Generation Sequencing (NGS).
2.  **Alignment:** Mapping the raw sequence reads to a reference genome. This step identifies where each read originated from in the reference sequence.
3.  **Variant Identification:** Analyzing the aligned reads to detect positions where the sample's sequence differs from the reference. Sophisticated algorithms are used to distinguish true genetic variants from sequencing errors or alignment artifacts. The output of this step is often a VCF file.

### Research Purposes and Applications of VCF-based Variant Calling

VCF files and the variant data they contain are crucial for a wide range of research and clinical applications:

*   **Population Genetics:**
    *   Studying genetic diversity within and between populations.
    *   Inferring population structure, migration patterns, and demographic history.
    *   Identifying regions of the genome under selection.
*   **Disease Association Studies:**
    *   Identifying genetic variants associated with common or rare diseases (e.g., Genome-Wide Association Studies - GWAS).
    *   Discovering variants that contribute to disease susceptibility or drug response.
*   **Personal Genomics:**
    *   Determining an individual's genetic ancestry.
    *   Assessing predisposition to certain traits or conditions.
    *   Pharmacogenomics: predicting how an individual will respond to specific drugs based on their genetic makeup.
*   **Evolutionary Studies:**
    *   Comparing genomes across different species to understand evolutionary relationships.
    *   Identifying genetic changes that drive adaptation and speciation.
*   **Cancer Genomics:**
    *   Identifying somatic mutations (variants that occur in tumor cells but not in normal cells) to understand cancer development.
    *   Guiding targeted cancer therapies based on the tumor's genetic profile.

### The `analyze_vcf.py` Script

The `analyze_vcf.py` script provided in the scripts/ directory is a basic tool for working with VCF files. It allows users to filter variants based on two common criteria:
*   **Variant Quality (QUAL):** The minimum confidence score for a variant call.
*   **Read Depth (DP):** The minimum number of sequence reads covering a variant position.

Filtering is often one of the first steps in VCF analysis pipelines to remove low-confidence calls or variants with insufficient evidence, thereby improving the reliability of downstream analyses. This script can output a new, filtered VCF file or a summary report of the variants that pass the specified filters.
