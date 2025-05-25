# Guide to Virus Genomics, Diversity, and Analysis

**Table of Contents**
*   [Introduction to Virus Genomics](#introduction-to-virus-genomics)
    *   [Importance of Virus Genomics](#importance-of-virus-genomics)
    *   [Key Characteristics of Viral Genomes](#key-characteristics-of-viral-genomes)
*   [Key Databases and Web Resources for Viral Sequence Data](#key-databases-and-web-resources-for-viral-sequence-data)
    *   [Primary Sequence Repositories General](#primary-sequence-repositories-general)
    *   [NCBI Specialized Viral Resources](#ncbi-specialized-viral-resources)
    *   [Pathogen Specific and Specialized Databases](#pathogen-specific-and-specialized-databases)
    *   [Metagenomic and Environmental Viral Databases](#metagenomic-and-environmental-viral-databases)
    *   [Other Useful Resources](#other-useful-resources)
*   [Measuring Viral Diversity](#measuring-viral-diversity)
    *   [Common Metrics for Viral Diversity](#common-metrics-for-viral-diversity)
    *   [Key Concepts in Viral Diversity](#key-concepts-in-viral-diversity)
    *   [Methods for Assessing Diversity](#methods-for-assessing-diversity)
*   [Phylogenetic Analysis of Viruses](#phylogenetic-analysis-of-viruses)
    *   [Applications in Virology](#applications-in-virology)
    *   [Workflow and Tools Conceptual Overview](#workflow-and-tools-conceptual-overview)
    *   [Considerations for Viral Phylogenetics](#considerations-for-viral-phylogenetics)
*   [Core Bioinformatic Tasks and Python Implementations](#core-bioinformatic-tasks-and-python-implementations)
    *   [Data Acquisition from NCBI](#data-acquisition-from-ncbi)
    *   [Sequence Quality Control QC and Assembly](#sequence-quality-control-qc-and-assembly)
    *   [Genome Annotation](#genome-annotation)
    *   [Basic Sequence Analysis with Python](#basic-sequence-analysis-with-python)
    *   [Variant Calling and Analysis](#variant-calling-and-analysis)
    *   [Multiple Sequence Alignment MSA with Python Wrapper](#multiple-sequence-alignment-msa-with-python-wrapper)
    *   [Comparative Genomics](#comparative-genomics)
*   [Advanced Modeling Population Dynamics and Immunogenetics](#advanced-modeling-population-dynamics-and-immunogenetics)
    *   [Population Level Virus Host Interactions](#population-level-virus-host-interactions)
    *   [Predicting Immunogenic Peptides and TCR Interactions](#predicting-immunogenic-peptides-and-tcr-interactions)
    *   [Deep Learning in Viral Evolution and Immunogenicity](#deep-learning-in-viral-evolution-and-immunogenicity)
*   [Tools and Resources](#tools-and-resources)
    *   [EvolCat Python Scripts for Viral Genomics](#evolcat-python-scripts-for-viral-genomics)
    *   [Key External Databases and Resources Recap](#key-external-databases-and-resources-recap)
    *   [Key External Software Recap](#key-external-software-recap)
*   [Example Workflow Analyzing a Viral Dataset](#example-workflow-analyzing-a-viral-dataset)
    *   [Step 1 Data Acquisition and Initial Assessment](#step-1-data-acquisition-and-initial-assessment)
    *   [Step 2 Sequence Preparation and Cleaning](#step-2-sequence-preparation-and-cleaning)
    *   [Step 3 Multiple Sequence Alignment MSA](#step-3-multiple-sequence-alignment-msa)
    *   [Step 4 Alignment Curation](#step-4-alignment-curation)
    *   [Step 5 Calculate Genetic Diversity](#step-5-calculate-genetic-diversity)
    *   [Step 6 Phylogenetic Tree Construction](#step-6-phylogenetic-tree-construction)
    *   [Step 7 Test for Selection using Site Specific dN dS Analysis](#step-7-test-for-selection-using-site-specific-dn-ds-analysis)
    *   [Step 8 Interpretation and Reporting](#step-8-interpretation-and-reporting)
*   [Concluding Thoughts](#concluding-thoughts)
*   [References](#references)

## Introduction to Virus Genomics
Virus genomics is the field of virology dedicated to studying the complete genetic material (genome) of viruses. This comprehensive analysis involves sequencing viral genomes and subsequently examining their structure (organization and architecture), function (how genes are expressed and proteins operate), and evolution (how viruses change over time and adapt to their hosts and environments). Understanding these interactions extends from the molecular to the population level, where ecological models such as predator-prey dynamics can offer insights into virus-host population fluctuations [Friedman, 2022](#references).

### Importance of Virus Genomics
The study of viral genomes is crucial for numerous reasons, impacting both basic research and applied sciences:
*   **Understanding Viral Evolution:** Genomics allows scientists to track how viruses change over time, infer their origins (e.g., zoonotic spillovers), and elucidate the evolutionary relationships between different viral strains or species. This is fundamental to understanding viral diversity and adaptation, including the mechanisms of mutation and recombination that drive viral change [Friedman, 2022](#references).
*   **Epidemiology and Public Health:** Viral genomics is a cornerstone of modern epidemiology. It enables the monitoring of viral outbreaks in real-time, identification of infection sources, and tracking of transmission pathways within populations (a field known as phylodynamics). This information is vital for informing public health responses.
*   **Pathogenesis:** By comparing genomes of virulent and attenuated strains, or by analyzing mutations that arise during infection, researchers can investigate how viruses cause disease. This includes identifying specific genes, genetic markers, or mutations associated with increased virulence, host tropism, or immune evasion. The interaction with the host immune system, particularly the presentation of viral peptides by MHC molecules and their recognition by T-cells, is a key area informed by genomics [Friedman, 2024](#references).
*   **Antiviral Drug Development:** Genomics helps identify potential viral targets for new antiviral drugs and is used to monitor the emergence of drug-resistant variants.
*   **Vaccine Development:** Designing effective vaccines relies on understanding viral antigens. Viral genomics helps identify these antigens, track their evolution, and evaluate vaccine efficacy. Predicting which viral peptides are immunogenic is a significant computational challenge [Friedman, 2024](#references).
*   **Virus Discovery:** High-throughput sequencing and metagenomic approaches allow for the identification of novel viruses.

### Key Characteristics of Viral Genomes
Viral genomes are remarkably diverse:
*   **Nature of Genetic Material:** DNA or RNA; single-stranded (ss) or double-stranded (ds); positive-sense, negative-sense, or ambisense RNA.
*   **Genome Size:** Typically small (kb to a few Mb), influencing evolutionary constraints like overlapping genes for packaging efficiency [Friedman, 2022](#references).
*   **Genome Structure:** Linear, circular, or segmented.
*   **Mutation Rates:** High, especially in RNA viruses due to polymerases often lacking proofreading. This drives rapid evolution and the formation of quasispecies.
*   **Genome Compactness:** Achieved through strategies like overlapping ORFs, polycistronic mRNAs, and alternative splicing.

## Key Databases and Web Resources for Viral Sequence Data

Access to comprehensive and well-curated databases is fundamental for viral genomics research. Below is a list of major resources providing viral sequences, associated metadata, and often, analysis tools.

### Primary Sequence Repositories General

These are large, international collaborations that archive and provide access to a vast range of nucleotide and protein sequences, including viral data.

1.  **NCBI GenBank**
    *   **URL:** [https://www.ncbi.nlm.nih.gov/genbank/](https://www.ncbi.nlm.nih.gov/genbank/)
    *   **Description:** A comprehensive public database of nucleotide sequences for over 400,000 formally described species, including a massive collection of viral genomes and partial sequences. It's part of the International Nucleotide Sequence Database Collaboration (INSDC).
    *   **Viral Data Access:** Sequences can be searched directly via Entrez Nucleotide or Protein, or through specialized NCBI viral resources (see below).

2.  **European Nucleotide Archive (ENA)**
    *   **URL:** [https://www.ebi.ac.uk/ena/browser/home](https://www.ebi.ac.uk/ena/browser/home)
    *   **Description:** Another member of the INSDC, ENA provides a comprehensive record of the world's nucleotide sequencing information, covering raw reads, assembled sequences, and annotations. Data is synchronized with GenBank and DDBJ.
    *   **Viral Data Access:** Searchable interface, programmatic access, and specialized data portals.

3.  **DNA Data Bank of Japan (DDBJ)**
    *   **URL:** [https://www.ddbj.nig.ac.jp/index-e.html](https://www.ddbj.nig.ac.jp/index-e.html)
    *   **Description:** The third member of the INSDC, DDBJ collects nucleotide sequence data primarily from Japanese researchers but is globally comprehensive due to data exchange.
    *   **Viral Data Access:** Similar search and retrieval capabilities as GenBank and ENA.

### NCBI Specialized Viral Resources

NCBI offers several portals specifically tailored for viral data, enhancing discoverability and providing analysis tools.

4.  **NCBI Virus**
    *   **URL:** [https://www.ncbi.nlm.nih.gov/labs/virus/vssi/](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/) (Viral Sequence Search Interface)
    *   **URL (Main Portal):** [https://www.ncbi.nlm.nih.gov/genome/viruses/](https://www.ncbi.nlm.nih.gov/genome/viruses/)
    *   **Description:** A central resource that provides access to viral sequence data from GenBank and RefSeq, along with tools for sequence analysis, genome browsing, and phylogenetic tree generation. It allows users to easily find and download complete viral genomes, segments, and proteins.
    *   **Features:** Faceted search, reference sequences (RefSeq), pre-computed alignments, and tree viewers.

5.  **NCBI RefSeq (Reference Sequence Collection)**
    *   **URL:** [https://www.ncbi.nlm.nih.gov/refseq/](https://www.ncbi.nlm.nih.gov/refseq/)
    *   **Description:** While not exclusively viral, RefSeq provides a curated, non-redundant set of sequences, including many viral reference genomes. These are often used as standards for annotation and comparative genomics.
    *   **Viral Data Access:** Searchable via Entrez, and heavily utilized by NCBI Virus.

6.  **Sequence Read Archive (SRA)**
    *   **URL:** [https://www.ncbi.nlm.nih.gov/sra](https://www.ncbi.nlm.nih.gov/sra)
    *   **Description:** A primary archive for raw sequencing data from next-generation sequencing (NGS) platforms. Essential for studies involving viral metagenomics, discovery, or intra-host diversity analysis where re-analysis of raw reads is needed.
    *   **Viral Data Access:** Contains raw reads from many viral sequencing projects.

### Pathogen Specific and Specialized Databases

These resources often focus on particular viruses or groups of viruses, providing highly curated data and specialized tools.

7.  **GISAID (Global Initiative on Sharing All Influenza Data)**
    *   **URL:** [https://www.gisaid.org/](https://www.gisaid.org/)
    *   **Description:** A critical global science initiative providing rapid access to genomic and associated metadata for influenza viruses, and more recently, SARS-CoV-2, RSV, MERS-CoV, and other epidemic-prone pathogens. Emphasizes data sharing with attribution.
    *   **Features:** Sequence data, phylogenetic analysis tools, epidemiological tracking. (Registration is typically required for data submission and detailed access).

8.  **Virus Pathogen Resource (ViPR)**
    *   **URL:** [https://www.viprbrc.org/](https://www.viprbrc.org/)
    *   **Description:** An NIAID-funded Bioinformatics Resource Center (BRC) dedicated to supporting research on viral pathogens of public health concern (e.g., Influenza, Flaviviruses like Zika and Dengue, Coronaviruses, Ebola).
    *   **Features:** Curated sequence data, genome annotations, comparative genomics tools, sequence alignment, phylogenetic tree construction, epitope prediction, 3D protein structure visualization. Integrates data from IRD (Influenza Research Database) and DRD (Dengue Research Database).

9.  **Influenza Research Database (IRD)**
    *   **URL:** [https://www.fludb.org/](https://www.fludb.org/) (Now largely integrated into ViPR)
    *   **Description:** A resource for influenza virus sequence and surveillance data, offering tools for sequence analysis, gene annotation, and epitope prediction.

10. **Los Alamos HIV Sequence Database**
    *   **URL:** [https://www.hiv.lanl.gov/](https://www.hiv.lanl.gov/)
    *   **Description:** A comprehensive, curated database of HIV and SIV sequences. Provides tools for sequence retrieval, alignment, phylogenetic analysis, epitope mapping, and drug resistance analysis. A long-standing and highly valued resource in HIV research.

11. **International Committee on Taxonomy of Viruses (ICTV)**
    *   **URL:** [https://ictv.global/](https://ictv.global/)
    *   **Description:** While not a sequence database per se, the ICTV is the official body responsible for the naming and classification of viruses. Their website provides the latest official viral taxonomy, which is essential for correctly identifying and categorizing viral sequences. They also link to exemplar genomes.

12. **ViralZone**
    *   **URL:** [https://viralzone.expasy.org/](https://viralzone.expasy.org/)
    *   **Description:** A web resource from the Swiss Institute of Bioinformatics (SIB) providing concise, expert-reviewed fact sheets on all known viral families/genera, covering molecular biology, replication cycle, host range, and links to sequence data. Excellent for understanding viral biology.

### Metagenomic and Environmental Viral Databases

These focus on viruses discovered through metagenomic sequencing of environmental or host-associated samples.

13. **IMG/VR (Integrated Microbial Genomes and Microbiomes - Viruses)**
    *   **URL:** [https://img.jgi.doe.gov/vr/](https://img.jgi.doe.gov/vr/)
    *   **Description:** A large, integrated database of viral sequences identified from metagenomic and metatranscriptomic datasets, primarily from the Joint Genome Institute (JGI). Focuses on uncultivated viruses.
    *   **Features:** Searchable database, tools for comparative analysis, and contextual metadata.

14. **MGnify (formerly EBI Metagenomics)**
    *   **URL:** [https://www.ebi.ac.uk/metagenomics/](https://www.ebi.ac.uk/metagenomics/)
    *   **Description:** Provides analysis and archiving of metagenomic data, including viral sequences identified within these datasets. Offers tools for functional and taxonomic analysis of metagenomes.

### Other Useful Resources

15. **Nextstrain**
    *   **URL:** [https://nextstrain.org/](https://nextstrain.org/)
    *   **Description:** Not a primary sequence database, but an incredibly powerful open-source project for real-time tracking of pathogen evolution using publicly available genomic data. Provides interactive visualizations of phylogenies integrated with geographic and temporal data for viruses like SARS-CoV-2, Influenza, Ebola, Zika, etc. Often links back to source databases for sequences.

## Measuring Viral Diversity
Measuring viral diversity is crucial for understanding the adaptability of viruses, their mechanisms of immune escape, the development of drug resistance, and for tracking their evolution and transmission patterns. Viral populations with high genetic diversity are more likely to possess variants that can adapt to new hosts, evade pre-existing or vaccine-induced immunity, or become resistant to antiviral therapies.

### Common Metrics for Viral Diversity

*   **Nucleotide Diversity (π):**
    *   **Definition:** Nucleotide diversity (often represented by the Greek letter pi, π) is defined as the average number of nucleotide differences per site between any two sequences chosen randomly from the sample population.
    *   **Interpretation:** It is a common measure of the extent of genetic variation within a population at the nucleotide level. A higher π value indicates greater genetic diversity.
    *   **EvolCat-Python Tool:** The script `pylib/scripts/viral_tools/calculate_nucleotide_diversity.py` can be used to compute nucleotide diversity from a FASTA alignment file.

*   **Haplotype Diversity (Hd):**
    *   **Definition:** A haplotype is a unique sequence in a given dataset. Haplotype diversity measures the uniqueness of these haplotypes in a sample. It is calculated based on the frequencies of these unique sequences and ranges from 0 (all sequences in the sample are identical) to 1 (all sequences in the sample are unique).
    *   **Interpretation:** This metric provides insight into the richness (number of different haplotypes) and evenness (distribution of haplotype frequencies) of the viral population.

*   **Genetic Distance:**
    *   **Definition:** Genetic distance quantifies the number of nucleotide or amino acid differences between two sequences. It can be calculated for pairs of sequences or as an average over a population.
    *   **Examples:**
        *   **p-distance:** This is the simplest measure, representing the proportion of sites at which two sequences differ. It is calculated by dividing the number of differing sites by the total length of the sequence compared.
        *   **Kimura 2-Parameter (K2P) distance:** This model corrects for the fact that multiple substitutions might have occurred at the same site but are observed as a single difference. It specifically distinguishes between transitional (purine to purine, or pyrimidine to pyrimidine) and transversional (purine to pyrimidine, or vice versa) substitution rates.
    *   **EvolCat-Python Tools:**
        *   `pylib/scripts/calculate_dna_distances.py` can compute pairwise distances using various substitution models, including p-distance and others.
        *   `pylib/scripts/calculate_k2p.py` is specifically designed to calculate pairwise Kimura 2-Parameter distances.

### Key Concepts in Viral Diversity

*   **Quasispecies:** As described in [Friedman, 2022](#references), viral populations, especially RNA viruses, exist as a dynamic cloud of related variants. This heterogeneity is central to their adaptability.
*   **Viral Evolution Rates:** Influenced by polymerase fidelity, genome type, replication rate, and host selective pressures. High mutation rates contribute significantly to genetic diversity.
*   **Population Dynamics:** The interaction between virus and host populations can be modeled using ecological frameworks (e.g., Lotka-Volterra predator-prey models), where factors like host immunity decay and viral evolution influence population cycles [Friedman, 2022](#references). Spatial heterogeneity in host populations can also impact viral persistence.

### Methods for Assessing Diversity

*   **Sequence Alignment:**
    *   A critical prerequisite for most diversity analyses is an accurate multiple sequence alignment (MSA). The alignment places homologous residues (nucleotides or amino acids) from different sequences into the same columns, allowing for meaningful site-by-site comparisons. Errors in alignment can lead to grossly inaccurate diversity estimates.

*   **Software Tools:**
    *   While EvolCat-Python provides scripts for calculating specific diversity metrics (as mentioned above), more comprehensive population genetics and viral diversity analyses often require specialized software packages. Examples include:
        *   **DnaSP (DNA Sequence Polymorphism):** Widely used for analyzing nucleotide polymorphism from population sequence data, calculating various diversity indices, and performing neutrality tests.
        *   **Arlequin:** A comprehensive population genetics software package that can compute diversity indices, population structure, and perform tests of selection.
        *   **VAPiD (Viral Analysis Pipeline for Intrahost Diversity):** Specifically designed for analyzing viral intra-host population diversity from next-generation sequencing data.

## Phylogenetic Analysis of Viruses
Phylogenetic analysis is a powerful bioinformatic approach used to study the evolutionary history and relationships of organisms, including viruses. In virology, it involves comparing the genetic sequences of different viral isolates to infer an evolutionary tree (phylogeny) that depicts their relatedness. These trees show how viruses have diverged from common ancestors over time.

### Applications in Virology
Phylogenetic methods have a wide array of applications in understanding viral biology, evolution, and epidemiology:

*   **Tracking Outbreak Origins and Spread (Phylodynamics):** This is a major application, especially for rapidly evolving RNA viruses. By analyzing sequences collected during an epidemic, researchers can:
    *   Determine the likely geographic origin of an outbreak.
    *   Track how the virus spreads geographically and chronologically.
    *   Identify key transmission events and patterns (e.g., superspreading events).
    *   Estimate the effective reproductive number (R0/Re) of the virus.
    This field, often called phylodynamics, is crucial for informing public health interventions.

*   **Identifying Transmission Chains:** In specific outbreak investigations (e.g., within a hospital or community), phylogenetics can help understand who likely infected whom, providing valuable information for contact tracing and targeted public health measures.

*   **Understanding Evolutionary Relationships:** Phylogenetics clarifies the evolutionary relationships between different viral strains, species, or even higher taxonomic groups (families, orders). This aids in viral classification, understanding the origins of new viruses (e.g., through zoonotic spillover), and mapping out broader evolutionary patterns. It's important to note that constructing virus phylogenies can be more complex than for cellular organisms due to high evolution rates and lack of universally conserved genes across all viral types [Friedman, 2022](#references).

*   **Studying the Evolution of Virulence or Drug Resistance:** By mapping phenotypic traits (like disease severity, host range, or resistance to antiviral drugs) onto a phylogenetic tree, researchers can correlate these traits with specific genetic changes (mutations or acquisition of genes). This helps identify genetic determinants of virulence or resistance.

*   **Vaccine Strain Selection:** For viruses like influenza that undergo rapid antigenic evolution, phylogenetic analysis of currently circulating strains helps predict which variants are likely to become dominant in the near future. This information is critical for selecting appropriate strains to include in seasonal vaccine formulations.

*   **Recombination Detection:** Recombination (exchange of genetic material between different viral genomes co-infecting the same cell) can be a significant evolutionary mechanism for some viruses. Phylogenetic incongruence, where trees built from different parts of the viral genome show conflicting evolutionary histories, can be a strong indicator of past recombination events.

### Workflow and Tools Conceptual Overview
A typical viral phylogenetic analysis workflow involves several steps:

*   **1. Data Preparation:**
    *   **Sequence Quality and Metadata:** High-quality viral sequences are essential. Equally important is accurate metadata, including the time of sampling, geographic location, host species, etc. This metadata is crucial for interpreting the phylogenetic tree in an epidemiological or ecological context.
    *   EvolCat-Python scripts like `pylib/scripts/clean_fasta_name.py` and `pylib/scripts/merge_fastas.py` can aid in pre-processing.

*   **2. Multiple Sequence Alignment (MSA):**
    *   Align homologous sequences using tools like MAFFT, MUSCLE, or Clustal Omega.

*   **3. Alignment Curation:**
    *   Inspect and curate alignments using viewers like Jalview or AliView. `pylib/scripts/nogaps.py` can help remove entirely gappy columns.

*   **4. Phylogenetic Tree Inference:**
    *   Use methods like Maximum Likelihood (ML) or Bayesian Inference (BI) with software like IQ-TREE, RAxML, MrBayes, or BEAST. `pylib/scripts/fas2phy.py` can help with format conversion if needed.

*   **5. Tree Visualization and Interpretation:**
    *   Visualize trees with FigTree, iTOL, or Dendroscope. For Python-based tree manipulation and basic interpretation, refer to the [Interpreting Phylogenetic Trees with Python](./phylogenetic-tree-interpretation.html) guide.

### Considerations for Viral Phylogenetics

*   **High Mutation Rates:** Can lead to homoplasy; appropriate substitution models are key.
*   **Recombination:** Can violate the assumption of a single tree; may require specific detection methods.
*   **Substitution Models:** Crucial for accuracy; use model selection tools (e.g., ModelFinder in IQ-TREE).
*   **Sampling Bias:** Non-random sampling can affect inferences.

## Core Bioinformatic Tasks and Python Implementations

This section provides Python examples for several common tasks in viral genomics, primarily using Biopython.

### Data Acquisition from NCBI

```python
from Bio import Entrez, SeqIO

Entrez.email = "your_email@example.com" # ALWAYS provide your email

def fetch_viral_sequences_ncbi(search_term, database="nucleotide", retmax=5, output_filename="viral_sequences.fasta"):
    """
    Searches NCBI and fetches viral sequences.
    Example search_term: "Zaire ebolavirus"[Organism] AND complete genome
                         "SARS-CoV-2"[Organism] AND "Spike gene"[Gene Name]
    """
    print(f"[NCBI] Searching for: {search_term}")
    try:
        handle = Entrez.esearch(db=database, term=search_term, retmax=retmax, idtype="acc")
        record = Entrez.read(handle)
        handle.close()
        ids = record["IdList"]
        print(f"[NCBI] Found {len(ids)} IDs. Fetching up to {retmax}.")

        if not ids:
            print("[NCBI] No records found.")
            return None

        handle = Entrez.efetch(db=database, id=ids, rettype="fasta", retmode="text")
        fasta_sequences = handle.read()
        handle.close()

        with open(output_filename, "w") as f:
            f.write(fasta_sequences)
        print(f"[NCBI] Sequences saved to {output_filename}")
        return output_filename
    except Exception as e:
        print(f"[NCBI] Error: {e}")
        return None

# Example usage:
# fetch_viral_sequences_ncbi(
#    search_term='"Zaire ebolavirus"[Organism] AND complete genome AND 2014:2016[Publication Date]',
#    retmax=3,
#    output_filename="ebolavirus_genomes_2014-2016.fasta"
# )
```

### Sequence Quality Control QC and Assembly
*   **QC:** Tools like FastQC are standard. Python can parse FastQC reports.
*   **Assembly:** Typically done with specialized assemblers (e.g., SPAdes, Flye, or reference-based with Bowtie2/BWA then iVar/bcftools). Python is useful for scripting pipelines and parsing assembly statistics (e.g., from QUAST).

### Genome Annotation
*   Dedicated tools like Prokka, VAPiD, or RAST are commonly used.
*   Python can parse annotation files (GFF3, GenBank) and implement simple ORF finders (see below).

### Basic Sequence Analysis with Python

```python
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio import SeqIO

def analyze_viral_genome_details(fasta_file):
    """
    Performs basic analysis (length, GC, ORF finding) on viral sequences in a FASTA file.
    """
    print(f"\n--- Analyzing sequences in {fasta_file} ---")
    for record in SeqIO.parse(fasta_file, "fasta"):
        print(f"\n> ID: {record.id}")
        print(f"  Description: {record.description}")
        seq = record.seq
        print(f"  Length: {len(seq)} bp")
        print(f"  GC Content: {GC(seq):.2f}%")

        print("  Potential ORFs (minimal length 90aa, standard genetic code, forward strands):")
        min_protein_len = 90 

        for frame in range(3):
            trimmed_seq = seq[frame:]
            if len(trimmed_seq) % 3 != 0:
                trimmed_seq = trimmed_seq[:-(len(trimmed_seq)%3)]
            if not trimmed_seq: continue
            try:
                trans = trimmed_seq.translate(to_stop=False) 
                for protein_segment in trans.split('*'):
                    start_codon_pos = protein_segment.find('M')
                    if start_codon_pos != -1:
                        orf_protein = protein_segment[start_codon_pos:]
                        if len(orf_protein) >= min_protein_len:
                            # Approximate nucleotide start/end for illustration
                            # For precise ORF coordinates, more careful calculation relative to original sequence and frame is needed.
                            print(f"    - Frame {frame+1}: Approx. Length {len(orf_protein)} aa, Seq: {orf_protein[:20]}...")
            except Exception as e:
                print(f"    - Error during translation/ORF finding in frame {frame+1}: {e}")
```

### Variant Calling and Analysis
1.  **Alignment:** Align reads to a reference (Bowtie2, BWA, Minimap2) -> SAM/BAM.
2.  **Variant Calling:** Use bcftools, LoFreq, iVar -> VCF.
3.  **Python for VCF Analysis (PyVCF):**

```python
import vcf # Ensure: pip install pyvcf

def analyze_vcf_details(vcf_filepath, min_qual=30, min_dp=10):
    print(f"\n--- Analyzing VCF file: {vcf_filepath} ---")
    try:
        vcf_reader = vcf.Reader(open(vcf_filepath, 'r'))
        count = 0
        for record in vcf_reader:
            qual_pass = record.QUAL is not None and record.QUAL >= min_qual
            
            dp_pass = False
            depth_info = "N/A"
            if 'DP' in record.INFO and isinstance(record.INFO.get('DP'), (int, float)) and record.INFO['DP'] >= min_dp:
                dp_pass = True
                depth_info = f"INFO DP: {record.INFO['DP']}"
            elif hasattr(record, 'samples') and record.samples and 'DP' in record.samples[0].data._fields and \
                 isinstance(record.samples[0]['DP'], (int,float)) and record.samples[0]['DP'] >= min_dp:
                dp_pass = True
                depth_info = f"Sample {record.samples[0].sample} DP: {record.samples[0]['DP']}"
            
            if qual_pass and dp_pass:
                count += 1
                print(f"  Variant at {record.CHROM}:{record.POS} {record.REF} -> {record.ALT}")
                print(f"    Quality: {record.QUAL}, Depth: {depth_info}")
        print(f"\nFound {count} variants passing filters (QUAL>={min_qual}, DP>={min_dp}).")
    except Exception as e:
        print(f"Error processing VCF file {vcf_filepath}: {e}")
```

### Multiple Sequence Alignment MSA with Python Wrapper

```python
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
import os

def run_mafft_alignment_wrapper(input_fasta, output_alignment_clustal):
    if not os.path.exists(input_fasta):
        print(f"Error: Input FASTA file not found: {input_fasta}")
        return None
    mafft_cline = MafftCommandline(input=input_fasta) # Assumes MAFFT is in PATH
    print(f"[MAFFT] Running alignment for {input_fasta}...")
    try:
        stdout, stderr = mafft_cline()
        temp_fasta_aln = "temp_mafft_output.fasta"
        with open(temp_fasta_aln, "w") as f: f.write(stdout)
        AlignIO.convert(temp_fasta_aln, "fasta", output_alignment_clustal, "clustal")
        print(f"[MAFFT] Alignment saved to {output_alignment_clustal} (Clustal format)")
        os.remove(temp_fasta_aln)
        return output_alignment_clustal
    except Exception as e:
        print(f"[MAFFT] Error: {e}")
        return None
```

### Comparative Genomics
Involves comparing genome structures, gene content, and evolutionary dynamics. Often uses a combination of MSA, phylogenetics, and custom scripting to parse annotations and alignment results. Tools like Mauve can be useful for visualizing whole-genome alignments.

## Advanced Modeling Population Dynamics and Immunogenetics

Beyond basic sequence analysis and phylogenetics, advanced computational models are increasingly used to understand complex virus-host interactions, including population-level dynamics and the intricacies of immune recognition.

### Population Level Virus Host Interactions
As detailed by [Friedman, 2022](#references), the interaction between a pathogenic virus and its vertebrate host can be conceptualized at the population level using ecological models like the Lotka-Volterra predator-prey system. In this framework:
*   The virus acts as the "predator" and the susceptible host population as the "prey."
*   Population fluctuations are driven by rates of host birth/recovery (which can relate to decaying immunity or new susceptible individuals) and virus transmission/clearance.
*   Factors like host genetic heterogeneity and spatial distribution can influence system stability and viral persistence.
*   Evolutionary changes in the virus (e.g., to evade host immunity) or in the host (e.g., acquired immunity) are critical components that modulate these population dynamics.

### Predicting Immunogenic Peptides and TCR Interactions
A key aspect of the host immune response involves the presentation of viral peptides by Major Histocompatibility Complex (MHC) molecules to T-cell Receptors (TCRs). Predicting which peptides will be immunogenic is a major challenge.
*   **MHC-Peptide Binding:** While models exist to predict peptide binding to MHC molecules, true immunogenicity (i.e., eliciting a T-cell response) is harder to predict [Friedman, 2024](#references).
*   **TCR-pMHC Interaction:** The recognition of a peptide-MHC (pMHC) complex by a TCR is highly specific.
    *   Structural modeling of TCRs, like with **TCRBuilder2** (part of the ImmuneBuilder suite), aims to predict the 3D structure of TCRs, which is crucial for understanding their binding specificity [Friedman, 2024, Abanades et al., 2023](#references).
    *   Computational methods like **PanPep** attempt to predict TCR-antigen binding using meta-learning principles, even with limited prior data (zero-shot or few-shot learning) [Friedman, 2024, Gao et al., 2023](#references).
    *   Structural similarity metrics like **TM-score** and **RMSD** are used to compare predicted protein structures (e.g., of TCRs) against experimentally determined structures or other predictions [Friedman, 2024](#references).

### Deep Learning in Viral Evolution and Immunogenicity
Deep learning approaches are becoming increasingly important in virology and immunogenetics:
*   **Protein Structure Prediction:** Models like AlphaFold and specialized tools like ImmuneBuilder/TCRBuilder2 predict 3D structures of viral proteins and immune receptors [Friedman, 2024](#references).
*   **Modeling Viral Evolution:** Deep learning can be used to model the complex interplay of mutation, recombination, and selection that drives viral evolution in response to host immunity [Friedman, 2022](#references).
*   **Predicting Immunogenicity:** Transformer-based models and other deep learning architectures are being applied to predict peptide immunogenicity and TCR-pMHC interactions, aiming to capture higher-order features beyond simple sequence motifs [Friedman, 2024](#references).
*   **Generative Models:** Systems like ProtGPT2 can generate novel protein sequences, and with appropriate training data and prompting, could potentially explore viral protein sequence space or design immunogenic peptides [Friedman, 2022](#references).
*   **Reinforcement Learning:** This paradigm could be applied to simulate viral adaptation, where the virus (agent) makes "moves" (mutations) and receives "rewards" based on its ability to evade a modeled immune system [Friedman, 2024](#references).

## Tools and Resources

This section highlights specific EvolCat-Python scripts relevant to viral genomics and key external resources.

### EvolCat Python Scripts for Viral Genomics

*   **`calculate_site_specific_ds_dn.py`**:
    *   **Purpose:** A Python wrapper for the `codeml` program from the PAML package. It automates identifying natural selection at individual codon sites in aligned coding sequences.
    *   **dN/dS Ratio:**
        *   `dN/dS > 1`: Positive (Darwinian) selection.
        *   `dN/dS < 1`: Purifying (negative) selection.
        *   `dN/dS = 1`: Neutral evolution.
    *   **Inputs:** Coding sequence alignment (FASTA), phylogenetic tree (Newick), PAML model (e.g., M0, M1a, M2a, M8).
    *   **Output:** Generates several files, most importantly `<outfile_prefix>_site_analysis.tsv`, which includes site number, dN/dS, and Bayes Empirical Bayes (BEB) posterior probabilities for identifying sites under positive selection.
    *   **Dependency:** Requires PAML (`codeml` executable) to be installed and accessible. This tool is particularly relevant for understanding selective pressures discussed in [Friedman, 2024](#references) concerning pathogen adaptation.

*   **`pylib/scripts/viral_tools/calculate_nucleotide_diversity.py`**: Calculates nucleotide diversity (π) from a FASTA alignment.
*   **`pylib/scripts/calculate_dna_distances.py`**: Computes pairwise genetic distances between sequences using various substitution models.
*   **`pylib/scripts/calculate_k2p.py`**: Calculates pairwise Kimura 2-Parameter (K2P) distances.
*   **`pylib/scripts/clean_fasta_name.py`**: Standardizes FASTA sequence headers.
*   **`pylib/scripts/merge_fastas.py`**: Combines multiple FASTA files.
*   **`pylib/scripts/nogaps.py`**: Removes columns from an MSA that consist entirely or predominantly of gaps.
*   **`pylib/scripts/fas2phy.py`**: Converts sequence alignments from FASTA to PHYLIP format.
*   **(Conceptual) `pylib/scripts/extract_region.py`**: For extracting specific genomic regions (functionality may vary).

### Key External Databases and Resources Recap

*   NCBI Viral Genomes, GISAID, ViPR, Nextstrain.

### Key External Software Recap

Essential external software for tasks beyond simple Python scripting:
*   **Multiple Sequence Alignment (MSA):** MAFFT, MUSCLE, Clustal Omega.
*   **Phylogenetic Inference:** IQ-TREE, RAxML, MrBayes, BEAST.
*   **Advanced dN/dS Analysis:** PAML package (specifically `codeml`).
*   **Protein Structure:** AlphaFold, ImmuneBuilder/TCRBuilder2.
*   **Sequence QC & Assembly:** FastQC, Trimmomatic, SPAdes, Bowtie2, BWA, Samtools, iVar.
*   **Variant Calling:** bcftools, LoFreq.
*   **Annotation:** Prokka, RAST.
*   **Visualization:** FigTree, iTOL, Jalview, AliView, RasMol (for PDB viewing, relevant to [Friedman, 2024](#references)).

## Example Workflow Analyzing a Viral Dataset

This section outlines a hypothetical step-by-step workflow to illustrate how a researcher might analyze a small viral dataset.

### Assumptions:
*   You have a set of viral sequences in FASTA format.
*   You have access to a command-line environment where EvolCat-Python scripts and necessary external tools can be run.

---

### Step 1 Data Acquisition and Initial Assessment

*   **1.1. Obtain Viral Sequences:** From public databases or own sequencing (`raw_sequences.fasta`).
*   **1.2. Initial Quality Control (Conceptual):** If from raw reads, use tools like Trimmomatic, BWA, Samtools, iVar.
*   **1.3. Gather Metadata:** Sampling dates, locations, host species, etc.

---

### Step 2 Sequence Preparation and Cleaning

*   **2.1. Standardize FASTA Headers:**
    ```bash
    python3 pylib/scripts/clean_fasta_name.py raw_sequences.fasta > cleaned_sequences.fasta
    ```
*   **2.2. Merge Sequence Files (If Applicable):**
    ```bash
    # If sequences are in batch1.fasta, batch2.fasta
    # python3 pylib/scripts/merge_fastas.py batch1.fasta batch2.fasta > combined_sequences.fasta
    # For this workflow, let's assume 'cleaned_sequences.fasta' is our input, renamed to:
    # mv cleaned_sequences.fasta sequences_for_analysis.fasta
    ```
*   **2.3. Extract Region of Interest (Optional):**
    ```bash
    # Conceptual: python3 pylib/scripts/extract_region.py sequences_for_analysis.fasta --coords START-END > gene_sequences.fasta
    # Let's assume 'sequences_for_analysis.fasta' is the file to align.
    ```

---

### Step 3 Multiple Sequence Alignment MSA

*   Use an external MSA program like MAFFT:
    ```bash
    mafft --auto sequences_for_analysis.fasta > aligned_sequences.afa
    ```

---

### Step 4 Alignment Curation

*   **4.1. Remove Gappy Columns:**
    ```bash
    python3 pylib/scripts/nogaps.py aligned_sequences.afa > curated_alignment.afa
    ```
*   **4.2. Manual Inspection:** Use AliView, Jalview, or SeaView.

---

### Step 5 Calculate Genetic Diversity

*   Use `pylib/scripts/viral_tools/calculate_nucleotide_diversity.py`:
    ```bash
    python3 pylib/scripts/viral_tools/calculate_nucleotide_diversity.py curated_alignment.afa > nucleotide_diversity_report.txt
    ```
*   **Optional: Pairwise Distances:**
    ```bash
    # python3 pylib/scripts/calculate_dna_distances.py curated_alignment.afa > pairwise_distances.tsv
    # python3 pylib/scripts/calculate_k2p.py curated_alignment.afa > pairwise_k2p_distances.tsv
    ```

---

### Step 6 Phylogenetic Tree Construction

*   **6.1. Format Conversion (If Needed for specific tools):**
    ```bash
    # python3 pylib/scripts/fas2phy.py curated_alignment.afa > curated_alignment.phy
    ```
*   **6.2. Tree Inference (e.g., IQ-TREE):**
    ```bash
    iqtree -s curated_alignment.afa -m MFP -B 1000 -T AUTO --prefix viral_phylogeny
    ```
    *(This creates `viral_phylogeny.treefile`, among others)*
*   **6.3. Tree Visualization:** Use FigTree, iTOL. Refer to [Phylogenetic Tree Interpretation Guide](./phylogenetic-tree-interpretation.html).

---

### Step 7 Test for Selection using Site Specific dN dS Analysis

*   Use `pylib/scripts/viral_tools/calculate_site_specific_ds_dn.py` (requires PAML's `codeml`).
    ```bash
    # Ensure curated_alignment.afa contains in-frame coding sequences.
    python3 pylib/scripts/viral_tools/calculate_site_specific_ds_dn.py \
        --alignment curated_alignment.afa \
        --tree viral_phylogeny.treefile \
        --model M8 \
        --outfile_prefix viral_selection_m8 \
        --verbose 
    ```
    *(Outputs include `viral_selection_m8_site_analysis.tsv`)*

---

### Step 8 Interpretation and Reporting

*   Synthesize results from diversity, phylogeny, and selection analyses.
*   Draw conclusions about viral population structure, evolution, and adaptive pressures.
*   Document methods and interpretations.

---

## Concluding Thoughts
Viral genomics integrates diverse computational approaches, from sequence analysis to complex modeling of evolutionary and immunological processes. The interplay between viral evolution, host population dynamics [Friedman, 2022](#references), and the specifics of immune recognition [Friedman, 2024](#references) highlights the need for multi-scale models. Python and its ecosystem, complemented by specialized bioinformatics tools and emerging deep learning techniques, provide a powerful framework for advancing our understanding in this critical field.

Always remember to consult the documentation for specific tools and databases, as best practices and software versions evolve.

## References
1.  Friedman, R. (2022). A Hierarchy of Interactions between Pathogenic Virus and Vertebrate Host. *Symmetry*, 14(11), 2274. [https://doi.org/10.3390/sym14112274](https://doi.org/10.3390/sym14112274)
2.  Friedman, R. (2024). Techniques for Theoretical Prediction of Immunogenic Peptides. *Encyclopedia*, 4(1), 600-621. [https://doi.org/10.3390/encyclopedia4010038](https://doi.org/10.3390/encyclopedia4010038)
3.  Abanades, B., Wong, W.K., Boyles, F., Georges, G., Bujotzek, A., & Deane, C.M. (2023). ImmuneBuilder: Deep-Learning models for predicting the structures of immune proteins. *Communications Biology*, 6, 575. (Cited in Friedman, 2024)
4.  Gao, Y., Gao, Y., Fan, Y., Zhu, C., Wei, Z., Zhou, C., Chuai, G., Chen, Q., Zhang, H., & Liu, Q. (2023). Pan-Peptide Meta Learning for T-cell receptor-antigen binding recognition. *Nature Machine Intelligence*, 5, 236–249. (Cited in Friedman, 2024)

## Acknowledgements
This guide has been developed with significant assistance from an AI language model. I would like to acknowledge the contributions of **Gemini AI (v2.5 Pro)** in the creation of this document.

The AI's assistance was invaluable in:
*   **Structuring and Organizing Content:** Helping to outline the guide and create a navigable Table of Contents.
*   **Drafting and Composing Text:** Generating initial drafts for various sections, including explanations of complex topics, Python code examples, and lists of databases/tools.
*   **Markdown Formatting:** Assisting with the correct Markdown syntax for display on GitHub Pages, including the iterative troubleshooting of Table of Contents links.
*   **Integrating Research:** Incorporating details from the author's published academic papers into relevant sections of the guide.
*   **Iterative Refinement:** Responding to feedback and making revisions to improve clarity, accuracy, and completeness.

While the AI provided substantial support in these areas, the overall direction, core scientific insights, specific script development (for EvolCat-Python), and final editorial decisions for this guide remain the work of the human author. This collaborative effort aimed to produce a comprehensive and helpful resource for the community.
