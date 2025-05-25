# Guide to Virus Genomics, Diversity, and Analysis

## Introduction to Virus Genomics
Virus genomics is the field of virology dedicated to studying the complete genetic material (genome) of viruses. This comprehensive analysis involves sequencing viral genomes and subsequently examining their structure (organization and architecture), function (how genes are expressed and proteins operate), and evolution (how viruses change over time and adapt to their hosts and environments).

**Importance of Virus Genomics:**
The study of viral genomes is crucial for numerous reasons, impacting both basic research and applied sciences:
*   **Understanding Viral Evolution:** Genomics allows scientists to track how viruses change over time, infer their origins (e.g., zoonotic spillovers), and elucidate the evolutionary relationships between different viral strains or species. This is fundamental to understanding viral diversity and adaptation.
*   **Epidemiology and Public Health:** Viral genomics is a cornerstone of modern epidemiology. It enables the monitoring of viral outbreaks in real-time, identification of infection sources, and tracking of transmission pathways within populations (a field known as phylodynamics). This information is vital for informing public health responses, such as quarantine measures or targeted interventions.
*   **Pathogenesis:** By comparing genomes of virulent and attenuated strains, or by analyzing mutations that arise during infection, researchers can investigate how viruses cause disease. This includes identifying specific genes, genetic markers, or mutations associated with increased virulence, host tropism, or immune evasion.
*   **Antiviral Drug Development:** Genomics helps identify potential viral targets (e.g., essential enzymes or structural proteins) for the development of new antiviral drugs. Furthermore, it is used to monitor the emergence and spread of drug-resistant viral variants.
*   **Vaccine Development:** Designing effective vaccines relies heavily on understanding the viral antigens (proteins or parts of proteins) that elicit protective immune responses. Viral genomics helps identify these antigens, track their evolution (antigenic drift/shift), and evaluate the efficacy of vaccine candidates.
*   **Virus Discovery:** High-throughput sequencing and metagenomic approaches applied to various environmental (e.g., water, soil) or clinical (e.g., patient tissues) samples allow for the identification of novel viruses, expanding our knowledge of the virosphere and potentially uncovering new pathogens.

**Key Characteristics of Viral Genomes:**
Viral genomes are remarkably diverse and possess several unique characteristics that distinguish them from the genomes of cellular organisms:
*   **Nature of Genetic Material:** Viral genomes exhibit extraordinary diversity in their basic composition. They can be made of:
    *   Deoxyribonucleic acid (DNA) or Ribonucleic acid (RNA).
    *   Single-stranded (ss) or double-stranded (ds).
    *   For RNA viruses, the genome can be positive-sense (+ssRNA, directly translatable by host ribosomes), negative-sense (-ssRNA, must be transcribed into +ssRNA before translation), or ambisense (containing both positive and negative-sense regions).
*   **Genome Size:** Viral genomes are typically much smaller and more compact than those of cellular organisms. They range from just a few kilobases (kb) – for example, Parvoviruses (around 5 kb) or Picornaviruses (around 7-8 kb) – to several megabases (Mb) in the case of some "giant viruses" like Pandoraviruses (up to 2.5 Mb).
*   **Genome Structure:** The physical organization of viral genomes also varies:
    *   **Linear:** A single, linear piece of nucleic acid (e.g., Poxviruses, Herpesviruses, many RNA viruses like HIV).
    *   **Circular:** A circular molecule of nucleic acid (e.g., Papillomaviruses, Polyomaviruses, some bacteriophages).
    *   **Segmented:** The genome is divided into two or more physically separate nucleic acid molecules, all of which are required for a productive infection (e.g., Influenza virus has 8 RNA segments, Reoviruses have 10-12 dsRNA segments).
*   **Mutation Rates:** Viruses, particularly RNA viruses, generally have very high mutation rates. This is often because their polymerases (enzymes that copy their genomes) lack the proofreading mechanisms found in cellular organisms or DNA viruses. These high mutation rates lead to rapid evolution, the generation of extensive genetic diversity, and the existence of viral populations as "quasispecies" (complex mixtures of related but non-identical genomes).
*   **Genome Compactness:** To function with such small genomes, viruses have evolved remarkable strategies for genomic economy. These include:
    *   **Overlapping Open Reading Frames (ORFs):** Different proteins are encoded from the same stretch of DNA/RNA by using different reading frames.
    *   **Polycistronic mRNAs:** A single messenger RNA (mRNA) molecule can encode multiple distinct proteins (common in RNA viruses and bacteriophages).
    *   **Alternative Splicing:** (More common in DNA viruses and retroviruses) A single pre-mRNA transcript can be spliced in different ways to produce multiple distinct mRNAs, and thus multiple proteins.

## Measuring Viral Diversity
Measuring viral diversity is crucial for understanding the adaptability of viruses, their mechanisms of immune escape, the development of drug resistance, and for tracking their evolution and transmission patterns. Viral populations with high genetic diversity are more likely to possess variants that can adapt to new hosts, evade pre-existing or vaccine-induced immunity, or become resistant to antiviral therapies.

**Common Metrics for Viral Diversity:**

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

**Key Concepts in Viral Diversity:**

*   **Quasispecies:**
    *   **Definition:** Many viral populations, particularly RNA viruses with their high mutation rates and large population sizes, do not exist as a single, uniform genetic entity. Instead, they form a "quasispecies" – a complex, dynamic, and heterogeneous cloud of closely related but genetically distinct variants (haplotypes). These variants are genetically linked through mutation, competition, and selection, and the population is often centered around a consensus or master sequence.
    *   **Significance:** The quasispecies nature of viruses is a major driver of their adaptability. The presence of a diverse array of pre-existing variants allows the viral population to rapidly respond to changing environmental pressures, such as host immune responses or the introduction of antiviral drugs. A variant that can evade the immune system or resist a drug might already exist at low frequency and can then be selected for, becoming dominant.

*   **Viral Evolution Rates:**
    *   **Explanation:** Different viruses evolve at vastly different rates. This rate is influenced by multiple factors, including:
        *   The fidelity of their replication enzyme (polymerase): RNA polymerases often lack proofreading capabilities, leading to higher error rates than DNA polymerases.
        *   The type of genome (RNA viruses generally evolve faster than DNA viruses).
        *   The replication rate of the virus.
        *   The generation time.
        *   Host factors and selective pressures (e.g., immune system, antiviral drugs).
    *   **Impact on Diversity:** High rates of evolution, driven by frequent mutations, directly contribute to the generation of genetic diversity within viral populations.

**Methods for Assessing Diversity:**

*   **Sequence Alignment:**
    *   A critical prerequisite for most diversity analyses is an accurate multiple sequence alignment (MSA). The alignment places homologous residues (nucleotides or amino acids) from different sequences into the same columns, allowing for meaningful site-by-site comparisons. Errors in alignment can lead to grossly inaccurate diversity estimates.

*   **Software Tools:**
    *   While EvolCat-Python provides scripts for calculating specific diversity metrics (as mentioned above), more comprehensive population genetics and viral diversity analyses often require specialized software packages. Examples include:
        *   **DnaSP (DNA Sequence Polymorphism):** Widely used for analyzing nucleotide polymorphism from population sequence data, calculating various diversity indices, and performing neutrality tests.
        *   **Arlequin:** A comprehensive population genetics software package that can compute diversity indices, population structure, and perform tests of selection.
        *   **VAPiD (Viral Analysis Pipeline for Intrahost Diversity):** Specifically designed for analyzing viral intra-host population diversity from next-generation sequencing data.

## Phylogenetic Analysis of Viruses
Phylogenetic analysis is a powerful bioinformatic approach used to study the evolutionary history and relationships of organisms, including viruses. In virology, it involves comparing the genetic sequences of different viral isolates to infer an evolutionary tree (phylogeny) that depicts their relatedness. These trees show how viruses have diverged from common ancestors over time.

**Applications in Virology:**
Phylogenetic methods have a wide array of applications in understanding viral biology, evolution, and epidemiology:

*   **Tracking Outbreak Origins and Spread (Phylodynamics):** This is a major application, especially for rapidly evolving RNA viruses. By analyzing sequences collected during an epidemic, researchers can:
    *   Determine the likely geographic origin of an outbreak.
    *   Track how the virus spreads geographically and chronologically.
    *   Identify key transmission events and patterns (e.g., superspreading events).
    *   Estimate the effective reproductive number (R0/Re) of the virus.
    This field, often called phylodynamics, is crucial for informing public health interventions.

*   **Identifying Transmission Chains:** In specific outbreak investigations (e.g., within a hospital or community), phylogenetics can help understand who likely infected whom, providing valuable information for contact tracing and targeted public health measures.

*   **Understanding Evolutionary Relationships:** Phylogenetics clarifies the evolutionary relationships between different viral strains, species, or even higher taxonomic groups (families, orders). This aids in viral classification, understanding the origins of new viruses (e.g., through zoonotic spillover), and mapping out broader evolutionary patterns.

*   **Studying the Evolution of Virulence or Drug Resistance:** By mapping phenotypic traits (like disease severity, host range, or resistance to antiviral drugs) onto a phylogenetic tree, researchers can correlate these traits with specific genetic changes (mutations or acquisition of genes). This helps identify genetic determinants of virulence or resistance.

*   **Vaccine Strain Selection:** For viruses like influenza that undergo rapid antigenic evolution, phylogenetic analysis of currently circulating strains helps predict which variants are likely to become dominant in the near future. This information is critical for selecting appropriate strains to include in seasonal vaccine formulations.

*   **Recombination Detection:** Recombination (exchange of genetic material between different viral genomes co-infecting the same cell) can be a significant evolutionary mechanism for some viruses. Phylogenetic incongruence, where trees built from different parts of the viral genome show conflicting evolutionary histories, can be a strong indicator of past recombination events.

**Workflow and Tools (Conceptual Overview):**
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

**Considerations for Viral Phylogenetics:**

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

### Sequence Quality Control (QC) and Assembly
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

### Multiple Sequence Alignment (MSA) with Python Wrapper

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

## Tools and Resources

This section highlights specific EvolCat-Python scripts relevant to viral genomics and key external resources.

### EvolCat-Python Scripts for Viral Genomics

*   **`calculate_site_specific_ds_dn.py`**:
    *   **Purpose:** A Python wrapper for the `codeml` program from the PAML package. It automates identifying natural selection at individual codon sites in aligned coding sequences.
    *   **dN/dS Ratio:**
        *   `dN/dS > 1`: Positive (Darwinian) selection.
        *   `dN/dS < 1`: Purifying (negative) selection.
        *   `dN/dS = 1`: Neutral evolution.
    *   **Inputs:** Coding sequence alignment (FASTA), phylogenetic tree (Newick), PAML model (e.g., M0, M1a, M2a, M8).
    *   **Output:** Generates several files, most importantly `<outfile_prefix>_site_analysis.tsv`, which includes site number, dN/dS, and Bayes Empirical Bayes (BEB) posterior probabilities for identifying sites under positive selection.
    *   **Dependency:** Requires PAML (`codeml` executable) to be installed and accessible.

*   **`pylib/scripts/viral_tools/calculate_nucleotide_diversity.py`**: Calculates nucleotide diversity (π) from a FASTA alignment.
*   **`pylib/scripts/calculate_dna_distances.py`**: Computes pairwise genetic distances between sequences using various substitution models.
*   **`pylib/scripts/calculate_k2p.py`**: Calculates pairwise Kimura 2-Parameter (K2P) distances.
*   **`pylib/scripts/clean_fasta_name.py`**: Standardizes FASTA sequence headers.
*   **`pylib/scripts/merge_fastas.py`**: Combines multiple FASTA files.
*   **`pylib/scripts/nogaps.py`**: Removes columns from an MSA that consist entirely or predominantly of gaps.
*   **`pylib/scripts/fas2phy.py`**: Converts sequence alignments from FASTA to PHYLIP format.
*   **(Conceptual) `pylib/scripts/extract_region.py`**: For extracting specific genomic regions (functionality may vary).

### Key External Databases and Resources

*   **NCBI Viral Genomes** ([https://www.ncbi.nlm.nih.gov/genome/viruses/](https://www.ncbi.nlm.nih.gov/genome/viruses/)): Comprehensive resource for viral sequences, reference genomes, and metadata.
*   **GISAID (Global Initiative on Sharing All Influenza Data)** ([https://www.gisaid.org/](https://www.gisaid.org/)): Critical platform for rapid sharing of genomic data for epidemic-prone viruses like influenza and SARS-CoV-2.
*   **Virus Pathogen Resource (ViPR)** ([https://www.viprbrc.org/](https://www.viprbrc.org/)): Provides curated viral sequence data and analysis tools.
*   **Nextstrain** ([https://nextstrain.org/](https://nextstrain.org/)): Open-source project for interactive visualizations of pathogen evolution.

### Key External Software (Recap)

Essential external software for tasks beyond simple Python scripting:
*   **Multiple Sequence Alignment (MSA):** MAFFT, MUSCLE, Clustal Omega.
*   **Phylogenetic Inference:** IQ-TREE, RAxML, MrBayes, BEAST.
*   **Advanced dN/dS Analysis:** PAML package (specifically `codeml`).
*   **Sequence QC & Assembly:** FastQC, Trimmomatic, SPAdes, Bowtie2, BWA, Samtools, iVar.
*   **Variant Calling:** bcftools, LoFreq.
*   **Annotation:** Prokka, RAST.
*   **Visualization:** FigTree, iTOL, Jalview, AliView.

## Example Workflow: Analyzing a Viral Dataset

This section outlines a hypothetical step-by-step workflow to illustrate how a researcher might analyze a small viral dataset.

**Assumptions:**
*   You have a set of viral sequences in FASTA format.
*   You have access to a command-line environment where EvolCat-Python scripts and necessary external tools can be run.

---

**Step 1: Data Acquisition & Initial Assessment**

*   **1.1. Obtain Viral Sequences:** From public databases or own sequencing (`raw_sequences.fasta`).
*   **1.2. Initial Quality Control (Conceptual):** If from raw reads, use tools like Trimmomatic, BWA, Samtools, iVar.
*   **1.3. Gather Metadata:** Sampling dates, locations, host species, etc.

---

**Step 2: Sequence Preparation and Cleaning**

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

**Step 3: Multiple Sequence Alignment (MSA)**

*   Use an external MSA program like MAFFT:
    ```bash
    mafft --auto sequences_for_analysis.fasta > aligned_sequences.afa
    ```

---

**Step 4: Alignment Curation**

*   **4.1. Remove Gappy Columns:**
    ```bash
    python3 pylib/scripts/nogaps.py aligned_sequences.afa > curated_alignment.afa
    ```
*   **4.2. Manual Inspection:** Use AliView, Jalview, or SeaView.

---

**Step 5: Calculate Genetic Diversity**

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

**Step 6: Phylogenetic Tree Construction**

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

**Step 7: Test for Selection (Site-Specific dN/dS Analysis)**

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

**Step 8: Interpretation and Reporting**

*   Synthesize results from diversity, phylogeny, and selection analyses.
*   Draw conclusions about viral population structure, evolution, and adaptive pressures.
*   Document methods and interpretations.

---

## Concluding Thoughts

Viral genomics is a dynamic field. This guide provides a starting point for using Python and associated tools in your viral genomics workflows. Many tasks involve integrating Python scripts with specialized command-line bioinformatics software. The ability to parse various file formats, automate repetitive tasks, and analyze large datasets makes Python an indispensable skill for viral genomicists.

Always remember to consult the documentation for specific tools and databases, as best practices and software versions evolve.
