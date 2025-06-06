# Virus Biology and Analysis Scripts

This directory contains Python scripts for processing viral mutation data.

## `proof_of_concept/` Directory

The `proof_of_concept/` subdirectory contains Python scripts related to a Transformer-based model for predicting viral mutations:

*   `transformer_model.py`: Defines the Transformer model architecture.
*   `train_transformer.py`: Script for training the Transformer model.
*   `predict_mutations.py`: Script for using a trained model to make predictions.

**Note:** These Transformer-related scripts were part of an experimental implementation. While the model architecture is defined, there were unresolved issues during the training phase (`ValueError` related to Keras argument passing). They are provided as a proof of concept and may require further debugging and development to be fully operational for training and reliable prediction. The primary, stable functionality offered by this suite is the FASTA sequence generation via `prepare_transformer_data.py`.
