# SARS-CoV-2 Lineage Classification Pipeline

## Overview
This document describes an end-to-end pipeline for processing SARS-CoV-2 mutation data, classifying fine-grained lineages using a neural network, and interpreting the model's decisions through saliency analysis. The goal is to identify key mutations that define these lineages, potentially corresponding to known biological markers.

The pipeline involves several stages:
1.  **Setup and Data Acquisition**: Preparing the environment and downloading necessary data, including the UShER Mutation-Annotated Tree (MAT).
2.  **Clade-Specific Data Extraction**: Isolating data for specific clades of interest into VCF format.
3.  **Sequence and Variant Processing**: Generating consensus sequences and programmatically accessing VCF data.
4.  **Feature Engineering**: Transforming raw variant data into a numerical matrix suitable for machine learning, with significant attention to memory optimization and data integrity.
5.  **Modeling**: Building, training, and evaluating a neural network classifier.
6.  **Interpretation**: Using saliency mapping to understand which mutations are most influential in the model's classification decisions.

This pipeline is designed to be run in an environment like Google Colab, leveraging its capabilities for package installation and distributed computing, though it can be adapted for other environments.

## Python Script: `sars_cov2_lineage_classifier.py`
[Script details will be added here]

### Dependencies and Execution Notes
*   **Environment**: This script is intended to be run in a Google Colab environment due to its use of Colab-specific commands (`drive.mount`, `!pip install`, `!conda install`).
*   **Google Drive**: The script assumes that input data (like `public-latest.all.masked.pb` and `public-latest.metadata.tsv`) and output files will be stored in Google MyDrive. Paths might need adjustment based on your Drive structure.
*   **Conda and UShER**: Conda is installed within the Colab environment to manage the installation of the UShER toolkit.
*   **Python Libraries**:
    *   `google.colab`
    *   `os`
    *   `condacolab`
    *   `pysam`
    *   `pandas`
    *   `gc` (garbage collector)
    *   `tensorflow`
    *   `sklearn`
    *   `numpy`
    *   `matplotlib`
    *   `seaborn`
    These are installed via `pip` or `conda` within the script.
*   **Input Files**:
    *   `public-latest.all.masked.pb.gz`: UShER MAT data (downloaded).
    *   `NC_045512v2.fa.gz`: SARS-CoV-2 reference sequence (downloaded).
    *   `public-latest.metadata.tsv`: Metadata file with sample IDs and lineages (expected to be in Google Drive).
*   **Output Files**:
    *   `my_clade.vcf`, `my_clade.vcf.gz`, `my_clade.vcf.gz.tbi`: VCF files for the clade of interest.
    *   `consensus.fa`: Consensus FASTA sequence.
    *   `confusion_matrix_top_N.png`: Image of the confusion matrix.
    *   `saliency_report_top10.csv`: CSV file with saliency analysis results.

## Utility and Critique of the Pipeline
### 1. Analysis of Methodology

This Python pipeline is robust, logical, and follows best practices.

*   **Part 1-4: Setup & Data Acquisition:** Standard and correct. It uses `conda` for environment management and `wget` to pull the latest, most comprehensive Mutation-Annotated Tree (MAT) from UShER. This is the right foundation.
*   **Part 5: Clade-Specific Data Extraction:** This is a critical step for making the problem manageable. Using `matUtils extract` with the `--clade` and `-v` (VCF output) flags is the perfect way to isolate the data for your clades of interest (e.g., `20H`, `20J`, `23E`). VCF is an excellent, standardized format for this.
*   **Part 6-7: Sequence & Variant Processing:** The use of `bcftools` to create a consensus sequence and `pysam` to programmatically access the VCF data is spot-on. Using `vcf.fetch()` on an indexed VCF (`my_clade.vcf.gz`) is highly efficient and the professionally recommended way to handle these files.
*   **Part 8: Feature Engineering:** This is arguably the most important part of the pipeline, and it's executed beautifully.
    *   **Memory Optimization:** Filtering the massive `metadata.tsv` *before* processing by keeping only the strains present in your VCF is a crucial optimization that prevents memory crashes. Using `gc.collect()` is also good practice.
    *   **Feature Matrix `X`:** Correctly identified that each unique mutation is a feature. Building the `variant_to_samples_map` and then creating the `(samples x mutations)` matrix is the right approach. Using `dtype='int8'` is a very smart memory-saving technique.
    *   **Singleton Removal:** **This is a critical and sophisticated step.** Removing lineages with only one sample is essential for using stratified splitting (`train_test_split`) and for building a generalizable model. A model can't learn to classify a group from a single example.
*   **Part 9-11: Modeling, Evaluation & Saliency:**
    *   **Data Splitting:** Using `train_test_split` with `stratify=y_np` is absolutely essential for this kind of imbalanced dataset. It ensures that the distribution of lineages is the same in your training and test sets, leading to a much more reliable evaluation.
    *   **Model Architecture:** A simple, dense neural network is a great starting point for this type of tabular data. The architecture (128 -> 64 -> num_classes) with ReLU and Softmax activations is a classic and effective choice.
    *   **Saliency Analysis:** This is the most advanced part. It correctly uses `tf.GradientTape` to calculate the gradient of the correct class output with respect to the input features. This directly answers the question: *"To correctly identify this sample, which mutations did the model 'pay attention' to the most?"* Generating a clean `.csv` report is far more useful for analysis than just plots.

---

### 2. Interpretation of Saliency Reports

The saliency reports are the payoff. They let us peer inside the "black box" of the neural network.

A key concept here is the **"importance" score**. A high score for a mutation means that changing this feature would have a large impact on the model's output for a specific lineage. In essence, it's a mutation the model has learned is a **strong and reliable identifier** for that lineage.

#### Example Analysis of `saliency_report_top10_23e.csv` (Clade XBB.2.3)

This report (from the sample run data) would show varying levels of importance for mutations across different sublineages of XBB.2.3.

*   **Highly Distinct Lineages:** For example, `XBB.2.3.13` (importance up to **0.8**), `HH.1` (**0.54**), and `XCQ` (**0.38**) are highlighted as extremely clear to the model in the provided analysis. Mutations like `5730_c` (XBB.2.3.13) and `20697_t` (HH.1) are identified as powerful identifiers.
*   **Broad Range:** The importance scores cover the full spectrum, from very high scores down to very subtle ones for lineages like `XBB.2.3.3` (`~1e-08` in the example).
*   **Lineage `HH.1`:** The top mutation, `20697_t`, has an importance of `0.54`. This mutation (T690I in the Spike protein, though the example later clarifies it as nsp15:T690I) is a known marker, and the model correctly identified it as a key feature. This is a great validation of the approach.

---

## Interpretation of Saliency Reports
[Interpretation of Saliency Reports will be added here]

## Overall Conclusions & Synthesis
[Overall Conclusions & Synthesis will be added here]

## Potential Next Steps
[Potential Next Steps will be added here]

## Example Run: Analysis of Clade 23E (XBB.2.3) Results
[Example Run details will be added here]

## Onto the Next Step: Biological Analysis
[Biological Analysis details will be added here]

## On Model Size and Further Experiments
[Model Size and Further Experiments details will be added here]

---
**The core goal of this project is:** "to develop an end-to-end pipeline that can accurately classify fine-grained SARS-CoV-2 lineages and, using saliency analysis, we can interpret the model's decisions to identify the key mutations that define these lineages, which correspond to known biological markers."
---
