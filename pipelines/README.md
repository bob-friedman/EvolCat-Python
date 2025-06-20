# Virus Biology and Analysis Scripts

This directory contains Python scripts for processing viral mutation data.

## `proof_of_concept/` Directory

This subdirectory contains an experimental implementation of a Transformer-based model for predicting viral mutations.

> **⚠️ Experimental Status**
> These scripts are a **proof of concept** and are not fully operational for training. They are provided for architectural reference and to showcase the primary, stable functionality of this suite: FASTA sequence generation.

### File Descriptions

*   `prepare_transformer_data.py`: **(Stable)** Preprocesses FASTA files into a format suitable for the Transformer model. This is the primary stable script in this directory.
*   `transformer_model.py`: Defines the Transformer model architecture using Keras.
*   `train_transformer.py`: Script for training the Transformer model.
*   `predict_mutations.py`: Script for using a trained model to make predictions.

### Known Issues & Areas for Development

*   **Training Halts:** The training process in `train_transformer.py` currently fails with a `ValueError` related to Keras argument passing.
*   **Contributions Welcome:** The model requires further debugging and development to be fully operational. We welcome contributions to resolve the training issues and improve the prediction pipeline.
