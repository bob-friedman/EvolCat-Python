<a name="top"></a>
# Guide to Virus Genomics, Diversity, and Analysis

This guide provides a comprehensive overview of the core concepts, databases, and bioinformatic workflows used in the study of viral genomics, with a focus on diversity, evolution, and practical analysis.

### Table of Contents
*   [Introduction to Virus Genomics](#introduction-to-virus-genomics)
*   [Key Databases & Resources](#key-databases-and-web-resources-for-viral-sequence-data)
*   [Measuring Viral Diversity](#measuring-viral-diversity)
*   [Phylogenetic Analysis of Viruses](#phylogenetic-analysis-of-viruses)
*   [Core Bioinformatic Tasks in Python](#core-bioinformatic-tasks-and-python-implementations)
*   [Viral Evolution Recombination and Reassortment](#viral-evolution-recombination-and-reassortment)
*   [Advanced Topics and Modeling](#advanced-topics-and-modeling)
    *   [Population Dynamics and Immunogenetics](#population-dynamics-and-immunogenetics)
    *   [Predictive Modeling of Viral Evolution](#predictive-modeling-of-viral-evolution)
    *   [Estimating Mutational Fitness Effects](#estimating-mutational-fitness-effects-from-large-scale-sequence-data)
    *   [Ancestral Sequence Based Viral Evolution Modeling](#ancestral-sequence-based-viral-evolution-modeling)
*   [Tools & Example Workflow](#tools-and-example-workflow)
*   [Acknowledgements and References](#acknowledgements-and-references)

---

## Introduction to Virus Genomics
Virus genomics is the field of virology dedicated to studying the complete genetic material (genome) of viruses. This comprehensive analysis involves sequencing viral genomes and subsequently examining their **structure** (organization and architecture), **function** (how genes are expressed and proteins operate), and **evolution** (how viruses change over time and adapt to their hosts and environments).

Understanding these interactions extends from the molecular to the population level, where ecological models such as predator-prey dynamics can offer insights into virus-host population fluctuations [1].

### Importance of Virus Genomics
The study of viral genomes is crucial for numerous reasons, impacting both basic research and applied sciences:
*   **Understanding Viral Evolution:** Genomics allows scientists to track how viruses change over time, infer their origins (e.g., zoonotic spillovers), and elucidate the evolutionary relationships between different viral strains or species. This is fundamental to understanding viral diversity and adaptation, including the mechanisms of mutation and recombination that drive viral change [1].
*   **Epidemiology and Public Health:** Viral genomics is a cornerstone of modern epidemiology. It enables the monitoring of viral outbreaks in real-time, identification of infection sources, and tracking of transmission pathways within populations (a field known as phylodynamics). This information is vital for informing public health responses.
*   **Pathogenesis:** By comparing genomes of virulent and attenuated strains, or by analyzing mutations that arise during infection, researchers can investigate how viruses cause disease. This includes identifying specific genes or mutations associated with increased virulence or immune evasion. The interaction with the host immune system, particularly the presentation of viral peptides by MHC molecules, is a key area informed by genomics [2].
*   **Antiviral Drug Development:** Genomics helps identify potential viral targets for new antiviral drugs and is used to monitor the emergence of drug-resistant variants.
*   **Vaccine Development:** Designing effective vaccines relies on understanding viral antigens. Viral genomics helps identify these antigens, track their evolution, and evaluate vaccine efficacy. Predicting which viral peptides are immunogenic is a significant computational challenge [2].
*   **Virus Discovery:** High-throughput sequencing and metagenomic approaches allow for the identification of novel viruses.

### Key Characteristics of Viral Genomes
Viral genomes are remarkably diverse:
*   **Nature of Genetic Material:** DNA or RNA; single-stranded (ss) or double-stranded (ds).
*   **Genome Size:** Typically small (kb to a few Mb), influencing evolutionary constraints like overlapping genes [1].
*   **Genome Structure:** Linear, circular, or segmented.
*   **Mutation Rates:** High, especially in RNA viruses due to polymerases often lacking proofreading. This drives rapid evolution and the formation of quasispecies.
*   **Genome Compactness:** Achieved through strategies like overlapping ORFs, polycistronic mRNAs, and alternative splicing.

[Back to Top](#top)

---
## Key Databases and Web Resources for Viral Sequence Data
Access to comprehensive and well-curated databases is fundamental for viral genomics research.

<details>
<summary><b>Click to expand the full list of major databases and resources</b></summary>

### Primary Sequence Repositories (General)
These are large, international collaborations that archive and provide access to a vast range of nucleotide and protein sequences, including viral data.

*   **NCBI GenBank:** [https://www.ncbi.nlm.nih.gov/genbank/](https://www.ncbi.nlm.nih.gov/genbank/)
    *   A comprehensive public database of nucleotide sequences, part of the International Nucleotide Sequence Database Collaboration (INSDC).
---
*   **European Nucleotide Archive (ENA):** [https://www.ebi.ac.uk/ena/browser/home](https://www.ebi.ac.uk/ena/browser/home)
    *   Another member of the INSDC, providing a comprehensive record of the world's nucleotide sequencing information. Data is synchronized with GenBank and DDBJ.
---
*   **DNA Data Bank of Japan (DDBJ):** [https://www.ddbj.nig.ac.jp/index-e.html](https://www.ddbj.nig.ac.jp/index-e.html)
    *   The third member of the INSDC, collecting nucleotide sequence data primarily from Japanese researchers but is globally comprehensive due to data exchange.

### NCBI Specialized Viral Resources
NCBI offers several portals specifically tailored for viral data, enhancing discoverability and providing analysis tools.

*   **NCBI Virus:** [https://www.ncbi.nlm.nih.gov/labs/virus/vssi/](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/)
    *   A central resource that provides access to viral sequence data from GenBank and RefSeq, along with tools for sequence analysis, genome browsing, and phylogenetic tree generation.
---
*   **NCBI RefSeq (Reference Sequence Collection):** [https://www.ncbi.nlm.nih.gov/refseq/](https://www.ncbi.nlm.nih.gov/refseq/)
    *   Provides a curated, non-redundant set of sequences, including many viral reference genomes.
---
*   **Sequence Read Archive (SRA):** [https://www.ncbi.nlm.nih.gov/sra](https://www.ncbi.nlm.nih.gov/sra)
    *   A primary archive for raw sequencing data from next-generation sequencing (NGS) platforms.

### Pathogen-Specific and Specialized Databases
These resources often focus on particular viruses, providing highly curated data and specialized tools.

*   **GISAID (Global Initiative on Sharing All Influenza Data):** [https://www.gisaid.org/](https://www.gisaid.org/)
    *   A critical global science initiative providing rapid access to genomic and associated metadata for influenza viruses, SARS-CoV-2, and other epidemic-prone pathogens.
---
*   **Virus Pathogen Resource (ViPR):** [https://www.viprbrc.org/](https://www.viprbrc.org/)
    *   An NIAID-funded resource dedicated to research on viral pathogens of public health concern (e.g., Influenza, Flaviviruses, Coronaviruses, Ebola).
---
*   **Los Alamos HIV Sequence Database:** [https://www.hiv.lanl.gov/](https://www.hiv.lanl.gov/)
    *   A comprehensive, curated database of HIV and SIV sequences with analysis tools.
---
*   **International Committee on Taxonomy of Viruses (ICTV):** [https://ictv.global/](https://ictv.global/)
    *   The official body responsible for the naming and classification of viruses. Essential for correctly categorizing viral sequences.
---
*   **ViralZone:** [https://viralzone.expasy.org/](https://viralzone.expasy.org/)
    *   Provides concise, expert-reviewed fact sheets on all known viral families/genera.

### Metagenomic and Environmental Viral Databases
*   **IMG/VR (Integrated Microbial Genomes and Microbiomes - Viruses):** [https://img.jgi.doe.gov/vr/](https://img.jgi.doe.gov/vr/)
    *   A large, integrated database of viral sequences identified from metagenomic and metatranscriptomic datasets.
---
*   **MGnify (formerly EBI Metagenomics):** [https://www.ebi.ac.uk/metagenomics/](https://www.ebi.ac.uk/metagenomics/)
    *   Provides analysis and archiving of metagenomic data, including viral sequences identified within these datasets.

### Other Useful Resources
*   **Nextstrain:** [https://nextstrain.org/](https://nextstrain.org/)
    *   A powerful open-source project for real-time tracking of pathogen evolution using publicly available genomic data, providing interactive phylogenetic visualizations.

</details>

[Back to Top](#top)

---
## Measuring Viral Diversity
Measuring viral diversity is crucial for understanding adaptability, immune escape, drug resistance, and evolution.

### Common Metrics for Viral Diversity
*   **Nucleotide Diversity (π):**
    > The average number of nucleotide differences per site between any two sequences chosen randomly from the sample population. A higher π value indicates greater genetic diversity.
    >
    > **EvolCat-Python Tool:** The script `pylib/scripts/calculate_nucleotide_diversity.py` can be used to compute π from a FASTA alignment.

*   **Haplotype Diversity (Hd):**
    > Measures the uniqueness of haplotypes (unique sequences) in a sample. It ranges from 0 (all sequences are identical) to 1 (all sequences are unique).

*   **Genetic Distance:**
    > Quantifies the number of nucleotide or amino acid differences between two sequences.
    > *   **p-distance:** The simplest measure; the proportion of sites at which two sequences differ.
    > *   **Kimura 2-Parameter (K2P) distance:** A model that corrects for multiple substitutions by distinguishing between transitions and transversions.
    >
    > **EvolCat-Python Tools:** `calculate_dna_distances.py` and `calculate_k2p.py`.

### Key Concepts in Viral Diversity
*   **Quasispecies:** As described in [1], viral populations, especially RNA viruses, exist as a dynamic cloud of related variants, which is central to their adaptability.
*   **Viral Evolution Rates:** Influenced by polymerase fidelity, genome type, and host selective pressures.
*   **Population Dynamics:** Interactions between virus and host populations can be modeled using ecological frameworks [1].

[Back to Top](#top)

---
## Phylogenetic Analysis of Viruses
Phylogenetic analysis uses genetic sequences to infer the evolutionary history and relationships of viruses, depicted in an evolutionary tree.

### Applications in Virology
*   **Tracking Outbreak Origins and Spread (Phylodynamics):** A major application for monitoring epidemics, identifying origins, and tracking transmission in real-time. However, reconstructing these origins can be highly challenging due to stochastic events and incomplete data, a problem central to the study of viral genesis [3].
*   **Identifying Transmission Chains:** Helps in contact tracing and targeted public health measures.
*   **Understanding Evolutionary Relationships:** Clarifies relationships between viral strains and aids in classification [1].
*   **Studying Evolution of Virulence or Drug Resistance:** Correlates genetic changes with phenotypic traits.
*   **Vaccine Strain Selection:** Helps predict which variants are likely to become dominant, crucial for seasonal vaccines like influenza.
*   **Recombination Detection:** Phylogenetic incongruence (conflicting trees from different genome parts) is a strong indicator of recombination.

### Workflow and Tools
A typical viral phylogenetic analysis involves:
1.  **Data Preparation:** Gathering high-quality sequences and accurate metadata.
2.  **Multiple Sequence Alignment (MSA):** Using tools like MAFFT, MUSCLE, or Clustal Omega.
3.  **Alignment Curation:** Manual inspection and trimming (e.g., using `nogaps.py`).
4.  **Phylogenetic Tree Inference:** Using methods like Maximum Likelihood (ML) or Bayesian Inference (BI) with software like IQ-TREE, RAxML, or BEAST.
5.  **Tree Visualization and Interpretation:** Using tools like FigTree or iTOL. Refer to the [Interpreting Phylogenetic Trees with Python](./phylogenetic-tree-interpretation.md) guide.

[Back to Top](#top)

---
## Core Bioinformatic Tasks and Python Implementations
This section provides Python examples for several common tasks in viral genomics, primarily using Biopython.

<details>
<summary><b>Click to expand Python code examples</b></summary>

### Data Acquisition from NCBI
```python
from Bio import Entrez, SeqIO

# Always provide your email to NCBI
Entrez.email = "your_email@example.com" 

def fetch_viral_sequences_ncbi(search_term, retmax=5, output_filename="viral_sequences.fasta"):
    """Searches NCBI and fetches viral sequences."""
    print(f"[NCBI] Searching for: {search_term}")
    try:
        handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=retmax, idtype="acc")
        record = Entrez.read(handle)
        ids = record["IdList"]
        print(f"[NCBI] Found {len(ids)} IDs. Fetching up to {retmax}.")
        if not ids: return

        handle = Entrez.efetch(db="nucleotide", id=ids, rettype="fasta", retmode="text")
        with open(output_filename, "w") as f:
            f.write(handle.read())
        print(f"[NCBI] Sequences saved to {output_filename}")
    except Exception as e:
        print(f"[NCBI] Error: {e}")

# Example usage:
# fetch_viral_sequences_ncbi(
#    search_term='"Zaire ebolavirus"[Organism] AND complete genome',
#    retmax=3,
#    output_filename="ebolavirus_genomes.fasta"
# )
```

### Basic Sequence Analysis with Python
```python
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio import SeqIO

def analyze_viral_genome_details(fasta_file):
    """Performs basic analysis on viral sequences in a FASTA file."""
    print(f"\n--- Analyzing sequences in {fasta_file} ---")
    for record in SeqIO.parse(fasta_file, "fasta"):
        print(f"\n> ID: {record.id}")
        print(f"  Length: {len(record.seq)} bp")
        print(f"  GC Content: {GC(record.seq):.2f}%")
```

### Variant Calling and Analysis
After aligning reads and calling variants to produce a VCF file, you can use EvolCat-Python for filtering.
```bash
# Filter input.vcf, keeping variants with QUAL >= 50 and DP >= 20
python pylib/scripts/analyze_vcf.py \
    --vcf_file input.vcf \
    --min_qual 50 \
    --min_dp 20 \
    --output_vcf filtered_variants.vcf \
    --summary_report variant_summary.tsv
```

### Multiple Sequence Alignment (MSA) with Python Wrapper
```python
from Bio.Align.Applications import MafftCommandline
import os

def run_mafft_alignment_wrapper(input_fasta, output_alignment):
    """A simple wrapper to run MAFFT alignment."""
    if not os.path.exists(input_fasta):
        print(f"Error: Input file not found: {input_fasta}")
        return
    
    # Assumes 'mafft' is in your system's PATH
    mafft_cline = MafftCommandline(input=input_fasta)
    print(f"[MAFFT] Running alignment for {input_fasta}...")
    try:
        stdout, stderr = mafft_cline()
        with open(output_alignment, "w") as f:
            f.write(stdout)
        print(f"[MAFFT] Alignment saved to {output_alignment}")
    except Exception as e:
        print(f"[MAFFT] Error: {e}")
```
</details>

[Back to Top](#top)

---
## Viral Evolution Recombination and Reassortment
While point mutations provide constant variation, genetic exchange via recombination and reassortment can cause significant, abrupt evolutionary shifts.

### Recombination: Creating Mosaic Viral Genomes
Recombination occurs when the viral polymerase switches templates during replication, creating a chimeric genome.
> **Case Study: SARS-CoV-2 Evolution**
> The origin of SARS-CoV-2 likely involved complex recombination events in animal reservoirs. Throughout the pandemic, ongoing recombination between circulating lineages (e.g., the emergence of the Omicron XBB subvariant from two other BA.2 lineages) has demonstrated how this mechanism can create novel combinations of mutations affecting transmissibility and immune evasion [4, 5].

> **Case Study: HIV-1 Circulating Recombinant Forms (CRFs)**
> Extensive recombination in HIV-1 has led to numerous CRFs (e.g., CRF01_AE, CRF02_AG), which are mosaic viruses derived from different parental subtypes. These CRFs now account for a significant proportion of global infections and impact diagnostics and drug resistance patterns [6].

### Reassortment: Segment Shuffling in Segmented Viruses
Reassortment is the shuffling of entire genome segments between segmented viruses co-infecting the same host cell.
> **Case Study: Origin of Pandemic Influenza Strains**
> Reassortment is the primary driver of antigenic shift in influenza A viruses, which can lead to pandemics. The 2009 H1N1 pandemic virus was a complex "quadruple reassortant" containing segments derived from swine, human, and avian influenza viruses, likely assembled in swine "mixing vessels" [7].

### Detecting Recombination and Reassortment Events
Identifying these events requires careful bioinformatic analysis, often looking for **phylogenetic incongruence**—where trees from different genome parts show conflicting histories.
*   **Software Tools:** RDP, GARD (HyPhy), SimPlot are commonly used.
*   **EvolCat-Python's Role:** Scripts for sequence manipulation (`clean_fasta_name.py`, etc.) and format conversion (`fas2phy.py`) are crucial for preparing data for these external tools.

[Back to Top](#top)

---
## Advanced Topics and Modeling

### Population Dynamics and Immunogenetics
*   **Virus-Host Interactions:** Can be conceptualized using ecological models like the Lotka-Volterra predator-prey system, where the virus is the "predator" and the susceptible host is the "prey" [1].
*   **Predicting Immunogenic Peptides:** A key challenge is predicting which viral peptides, when presented by MHC molecules, will elicit a T-cell response. This involves modeling TCR-pMHC interactions using tools like **TCRBuilder2** and **PanPep** [2].
*   **Deep Learning:** Models like AlphaFold are used for protein structure prediction, while others are being applied to model viral evolution and predict immunogenicity.

### Predictive Modeling of Viral Evolution
<details>
<summary><b>Click for a detailed overview of methods for predicting pathogen evolution</b></summary>

A review by [8] explores the evolving landscape of predicting pathogen evolution and immune evasion, with a focus on Artificial Intelligence (AI). These data-driven approaches, which combine phylogenetic trees with AI, are at the forefront of efforts to outmaneuver viral evolution [9].

#### Methods Described
*   **Phylogenetic Analysis:** Using evolutionary trees to identify mutations associated with successful lineages.
*   **Deep Mutational Scans (DMS):** High-throughput experiments to systematically evaluate the functional impact of numerous mutations.
*   **Variational Autoencoders (VAEs):** Deep learning models that encode complex data (like viral genomes) into lower-dimensional representations.
*   **Protein Language Models (PLMs):** Models that treat protein sequences as "sentences" to learn the underlying rules of protein function and evolution. The application of these models, particularly those based on the Transformer architecture, is central to modern biological sequence analysis [10].

#### Key Takeaways
*   The review emphasizes that the optimal approach depends on the specific context and available resources, often advocating for a **multi-strategic approach** that combines different methods.
*   **Data-driven decision-making** is paramount; the choice of method is guided by the availability and quality of data.
*   The paper highlights the transformative potential of AI in moving beyond traditional methods to tackle the complex problem of viral evolution.

*Critique*: While an excellent overview, the review could benefit from more technical depth on the algorithms, a more direct comparison of the methods' strengths and weaknesses, and a deeper discussion of ethical considerations.
</details>

### Estimating Mutational Fitness Effects from Large-Scale Sequence Data
<details>
<summary><b>Click for a detailed overview of fitness effect estimation</b></summary>

A powerful approach to understand the fitness consequences of mutations involves leveraging vast amounts of public sequence data, as exemplified by [11] for SARS-CoV-2.

#### Core Methodology
1.  **Calculate Expected Mutation Counts:** Determine how many times each possible mutation is *expected* to occur along a phylogeny assuming no selection.
2.  **Count Observed Mutations:** Count how many times each mutation is *actually observed* on a large phylogenetic tree (~7 million sequences in the study).
3.  **Estimate Fitness:** The ratio of observed to expected counts estimates the mutation's fitness effect.
    *   `Ratio < 1`: Deleterious (purifying selection)
    *   `Ratio ≈ 1`: Neutral
    *   `Ratio > 1`: Potentially beneficial

#### Key Findings & Significance
*   **Comprehensive Fitness Maps:** This computational framework generates detailed maps of mutational fitness effects across all viral proteins, which correlate well with experimental data.
*   **Informing Public Health:** These maps are valuable for assessing the impact of new variants and guiding the design of antivirals or vaccines by targeting highly constrained regions.
*   **Broad Applicability:** The method can be applied to any virus with sufficient sequence data.
*   **Interactive Resource:** The study's findings for SARS-CoV-2 are explorable at [jbloomlab.github.io/SARS2-mut-fitness/](https://jbloomlab.github.io/SARS2-mut-fitness/).
</details>

### Ancestral Sequence Based Viral Evolution Modeling
<details>
<summary><b>Click for a detailed overview on this topic</b></summary>

An end-to-end guided pipeline for training a Transformer-based sequence-to-sequence model on viral evolution data. It demonstrates data acquisition from UShER and NCBI, ancestral state reconstruction (ASR) with IQ-TREE, data preparation, model training with TensorFlow/Keras, and inference. The guide also discusses the computational and biological limitations encountered on standard cloud-based hardware and serves as a documented baseline for future research. (Credits: Jules AI for presentation of this work and Gemini Pro for assistance in pipeline design, code development and debugging, scientific explanation, and the drafting of this technical guide.

Refer to the technical guides for further information: [Ancestral Sequence-Based Viral Evolution Modeling](./PIPELINE.md) and [Interpretable Viral Evolution Modeling](./PIPELINE_2.md).
</details>

[Back to Top](#top)

---
## Tools and Example Workflow
<a name="tools-and-example-workflow"></a>
### EvolCat-Python Scripts for Viral Genomics
*   **`calculate_site_specific_ds_dn.py`**: A wrapper for PAML `codeml` to detect natural selection at individual codon sites by comparing non-synonymous (dN) and synonymous (dS) substitution rates.
*   **Diversity & Distance:** `calculate_nucleotide_diversity.py`, `calculate_dna_distances.py`, `calculate_k2p.py`.
*   **File Preparation:** `clean_fasta_name.py`, `merge_fastas.py`, `nogaps.py`, `fas2phy.py`.
*   **Variant Filtering:** `analyze_vcf.py` for basic filtering of VCF files.

### Example Workflow: Analyzing a Viral Dataset

This hypothetical workflow illustrates how to use these tools in practice.
1.  **Data Acquisition:** Obtain viral sequences and metadata.
2.  **Sequence Preparation:** Standardize headers and merge files.
    ```bash
    python3 pylib/scripts/clean_fasta_name.py raw_sequences.fasta > cleaned_sequences.fasta
    ```
3.  **Multiple Sequence Alignment (MSA):** Use an external tool like MAFFT.
    ```bash
    mafft --auto cleaned_sequences.fasta > aligned_sequences.afa
    ```
4.  **Alignment Curation:** Remove gappy columns.
    ```bash
    python3 pylib/scripts/nogaps.py aligned_sequences.afa > curated_alignment.afa
    ```
5.  **Calculate Genetic Diversity:**
    ```bash
    python3 pylib/scripts/calculate_nucleotide_diversity.py curated_alignment.afa > diversity_report.txt
    ```
6.  **Phylogenetic Tree Construction:** Use an external tool like IQ-TREE.
    ```bash
    iqtree -s curated_alignment.afa -m MFP -B 1000 -T AUTO --prefix viral_phylogeny
    ```
7.  **Test for Selection:** Use the EvolCat-Python wrapper for PAML.
    ```bash
    python3 pylib/scripts/paml_tools/calculate_site_specific_ds_dn.py \
        --alignment curated_alignment.afa \
        --tree viral_phylogeny.treefile \
        --model M8 \
        --outfile_prefix viral_selection_m8
    ```
8.  **Interpretation and Reporting:** Synthesize all results to draw conclusions.

[Back to Top](#top)

---
## Acknowledgements and References
<a name="acknowledgements-and-references"></a>
### Acknowledgements
This guide has been developed with significant assistance from the AI language model **Gemini (v2.5 Pro)**. The AI's assistance was invaluable in structuring content, drafting text, formatting for readability, and integrating research. While the AI provided substantial support, the overall direction, core scientific insights, specific script development (for EvolCat-Python), and final editorial decisions for this guide remain the work of the human author.

### References
1.  Friedman, R. (2022). A Hierarchy of Interactions between Pathogenic Virus and Vertebrate Host. *Symmetry*, 14(11), 2274. [https://doi.org/10.3390/sym14112274](https://doi.org/10.3390/sym14112274)
2.  Friedman, R. (2024). Techniques for Theoretical Prediction of Immunogenic Peptides. *Encyclopedia*, 4(1), 600-621. [https://doi.org/10.3390/encyclopedia4010038](https://doi.org/10.3390/encyclopedia4010038)
3.  Friedman, R. (2025). The Elusive Genesis: Stochasticity and the Challenge of Reconstructing Viral Origins. *Preprints*. [https://doi.org/10.20944/preprints202505.2277.v2 ](https://doi.org/10.20944/preprints202505.2277.v2 )
4.  Andersen, K. G., Rambaut, A., Lipkin, W. I., Holmes, E. C., & Garry, R. F. (2020). The proximal origin of SARS-CoV-2. *Nature Medicine*, 26, 450-452.
5.  Uriu, K., Ito, J., Zahradnik, J., Fujita, S., Kosugi, Y., Schreiber, G. (2023). Enhanced transmissibility, infectivity, and immune resistance of the SARS-CoV-2 omicron XBB.1.5 variant. *Lancet Infectious Diseases*, 23, 280-281.
6.  Robertson, D. L., Anderson, J. P., Bradac, J. A., Carr, J. K., Foley, B., Funkhouser, R. K., ... & Korber, B. (2000). HIV-1 nomenclature proposal. *Science*, 288, 55-56.
7.  Garten, R. J., Davis, C. T., Russell, C. A., Shu, B., Lindstrom, S., Balish, A., ... & Cox, N. J. (2009). Antigenic and Genetic Characteristics of Swine-Origin 2009 A(H1N1) Influenza Viruses Circulating in Humans. *Science*, 325, 197-201.
8.  Hamelin, D.J., Scicluna, M., Saadie, I., Mostefai, F., Grenier, J.C., Baron, C., ... & Hussin, J.G. (2024). Predicting pathogen evolution and immune evasion in the age of artificial intelligence. *Computational and Structural Biotechnology Journal*, *23*, 1370-1382. [https://doi.org/10.1016/j.csbj.2024.03.044](https://doi.org/10.1016/j.csbj.2024.03.044)
9.  Friedman, R. (2025). The Viral Chase: Outsmarting Evolution with Data Trees and AI Predictions. *Preprints*. [https://doi.org/10.20944/preprints202506.0456.v1 ](https://doi.org/10.20944/preprints202506.0456.v1 )
10. Friedman, R. (2025). Anatomy of a Transformer: An Essay on the Code and Concepts for Biological Sequence Analysis. *Preprints*. [https://doi.org/10.20944/preprints202506.0623.v1 ](https://doi.org/10.20944/preprints202506.0623.v1 )
11. Bloom, J. D., & Neher, R. A. (2023). Fitness effects of mutations to SARS-CoV-2 proteins. *Virus Evolution*, *9*(2), vead055. [https://doi.org/10.1093/ve/vead055](https://doi.org/10.1093/ve/vead055)

[Back to Top](#top)
