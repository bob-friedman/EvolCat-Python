# SARS-CoV-2 Lineage Classification Pipeline

This repository contains an end-to-end pipeline for processing SARS-CoV-2 mutation data, classifying fine-grained lineages with a high-accuracy neural network, and interpreting the model's decisions using saliency mapping to identify key biological markers.

This project demonstrates a complete workflow from raw genomic data to interpretable machine learning insights, achieving **over 97% accuracy** on a challenging 95-class classification problem.

---

## Key Features

*   **End-to-End Workflow:** A complete, runnable pipeline from data acquisition to model interpretation.
*   **Scalable Data Processing:** Efficiently handles massive datasets using `matUtils`, `pysam`, and memory-optimized `pandas` workflows.
*   **High-Accuracy Deep Learning Model:** A Keras/TensorFlow neural network designed for classifying fine-grained lineages.
*   **Model Interpretability:** Uses saliency mapping to identify the specific mutations the model finds most important for its classifications, providing a bridge between machine learning and virology.
*   **Colab-Ready:** Designed and documented to run seamlessly in a Google Colab environment.

---

## Pipeline Stages

The pipeline is organized into a logical sequence of stages:

1.  **Setup and Data Acquisition**: Prepares the Conda environment and downloads the latest UShER Mutation-Annotated Tree (MAT).
2.  **Clade-Specific Data Extraction**: Isolates variant data for specific clades of interest (e.g., XBB.2.3) into VCF format using `matUtils`.
3.  **Sequence and Variant Processing**: Generates consensus sequences with `bcftools` and programmatically accesses VCF data using `pysam`.
4.  **Feature Engineering**: Transforms raw variant data into a numerical `(samples x mutations)` matrix suitable for machine learning, with critical steps for memory optimization and data integrity.
5.  **Modeling**: Builds, trains, and evaluates a neural network classifier using Keras/TensorFlow.
6.  **Interpretation**: Generates saliency reports that highlight the most influential mutations for each lineage's classification.

---

## How to Run the Pipeline

This project is designed to be run as a single, comprehensive script within Google Colab.

**Prerequisites:**
*   A Google Account with access to Google Colab and Google Drive.
*   The primary metadata file (`public-latest.metadata.tsv`) should be present in your Google Drive root directory.

**Running the Pipeline:**
1.  Clone this repository or download the `sars_cov2_lineage_classifier.py` script.
2.  Upload the script to your Google Colab environment.
3.  Open and run the script from top to bottom. The script will handle all dependencies and data downloads automatically.

> **Note:** The entire script is a single file. All the numbered sections (1-11) from your original document are contained within it and will execute in order.

---

## Project Structure

*   `sars_cov_2_lineage_classifier.py`: The main Python script containing the entire end-to-end pipeline.
*   `matutils_implementation.md`: A deep dive into advanced `matUtils` usage and strategies for full-scale data extraction.
*   `proof_of_concept/`: Directory containing experimental Transformer models and related scripts.
*   `README.md`: This file.
*   `PAPER.md` : Technical report on this study.

---

## Results and Interpretation

This pipeline has been successfully tested on the 23E (XBB.2.3) clade, demonstrating its effectiveness on a complex, real-world dataset.

*   **Performance:** The model achieved **97.98% accuracy** on the test set, successfully distinguishing between 95 different sublineages.
*   **Insights:** The generated saliency reports successfully identified key mutations that define specific lineages. These computationally-derived markers often correspond to known biological markers documented in virology literature, validating the entire approach.

For a detailed analysis of the results, including the confusion matrix and a deep dive into the saliency report for lineage `HH.1`, please see the **[Analysis and Interpretation of the 23E Clade Run](./CODE_RESULTS.md)**.

---

## Credits

This work is made possible by the collaborative efforts of Jules and Gemini.
