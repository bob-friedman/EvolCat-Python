# Guides in the EvolCat-Python Project

> Welcome to the Guides page in the EvolCat-Python project. This section provides a collection of guides, reports, and tutorials covering everything from core virology concepts to advanced bioinformatic pipelines.

---

## 📖 Main Guides

These guides cover the core functionality and foundational concepts of the project.

*   [**Usage and Command-Line Options**](../docs/USAGE.md)
    *   *A comprehensive reference for all script arguments and options.*
*   [**Accessing MHC Sequence Databases**](./mhc-database-guide.md)
    *   *Step-by-step instructions for downloading and preparing MHC data.*
*   [**Interpreting Phylogenetic Trees with Python**](./phylogenetic-tree-interpretation.md)
    *   *A guide to using libraries like ETE3 and Biopython for tree analysis.*
*   [**Understanding Transformer Core Concepts**](./transformer_core_concepts.md)
    *   *An introduction to the attention mechanism and architecture used in the models.*
*   [**PyTorch vs. Keras for Transformers: A Comparison**](./pytorch_keras_transformer_comparison.md)
    *   *A discussion of the pros and cons of each framework for this project's goals.*
 
---

## 📖 Virus Biology Guides

This section contains conceptual guides and reports on key topics in virology and bioinformatics.

*   [**Virus Genomics, Diversity, and Analysis Guide**](./virus_genomics_guide.md)
    *   *A comprehensive overview of viral genomes, diversity metrics, phylogenetic analysis, and core bioinformatic tasks.*
*   [**Condensed Report on Virus Evolution and Biotechnological Innovation**](./condensed_virus_evolution_report.md)
    *   *A summary of key advancements in virology, with a focus on recombination and biotechnology.*
*   [**Guide to Key Biotechnologies in Virology Research**](./biotechnologies_in_virology_guide.md)
    *   *An overview of the foundational and cutting-edge technologies used to study viruses.*
*   [**Guide to Estimating Mutational Fitness Effects**](./estimating_mutation_fitness_effects_guide.md)
    *   *A methodology for estimating mutation fitness effects from large-scale sequence data.*
*   [**UShER Toolkit and matUtils: A Condensed Guide**](./usher_toolkit_report.md)
    *   *A practical guide to using the UShER toolkit for analyzing massive viral phylogenies like SARS-CoV-2.*

---

## ⚙️ Practical Pipelines and Tutorials

This section contains step-by-step tutorials and end-to-end pipelines for specific analyses.

*   [**Ancestral Reconstruction Tutorial**](./ancestral_reconstruction_tutorial.md)
    *   *A tutorial on performing ancestral sequence reconstruction using tools like TreeTime and DendroPy.*
*   [**SARS-CoV-2 Transformer Modeling Pipeline**](../pipelines/sars_cov2_transformer_pipeline.md)
    *   *A tutorial on processing SARS-CoV-2 data for use in Transformer-based generative models.*
*   [**SARS-CoV-2 Lineage Classification Pipeline**](../pipelines/sars_cov2_lineage_classification/README.md)
    *   *Documentation for the end-to-end pipeline that classifies SARS-CoV-2 lineages with a neural network and interprets the results.*
*   [**Guided Pipeline #1: Ancestral Sequence-Based Viral Evolution Modeling**](./PIPELINE.md)
    *   *An end-to-end guided pipeline for training a Transformer-based sequence-to-sequence model on viral evolution data. It demonstrates data acquisition from UShER and NCBI, ancestral state reconstruction (ASR) with IQ-TREE, data preparation, model training with TensorFlow/Keras, and inference. The guide also discusses the computational and biological limitations encountered on standard cloud-based hardware and serves as a documented baseline for future research. (Credits: Jules AI for presentation of this work and Gemini Pro for assistance in pipeline design, code development and debugging, scientific explanation, and the drafting of this technical guide.*
*   [**Guided Pipeline #2: Interpretable Viral Evolution Modeling**](./PIPELINE_2.md)
    *   *An end-to-end pipeline demonstrating how a Transformer model, trained on phylogenetically-derived ancestor-descendant SARS-CoV-2 Spike protein sequences, can autonomously learn to identify key functional regions (RBD and NTD) through its attention mechanism. This guide covers data acquisition, ancestral state reconstruction, model training, and attention-based interpretability. (Credits: Jules AI for presentation of this work and Gemini Pro for assistance in pipeline design, code development and debugging, scientific explanation, and the drafting of this technical guide.*

---

Navigate back to the [**main documentation page**](../README.md).
