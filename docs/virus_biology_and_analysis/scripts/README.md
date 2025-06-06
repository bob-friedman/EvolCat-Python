# Virus Biology and Analysis Scripts

This directory contains Python scripts for processing viral mutation data.

## `prepare_transformer_data.py` - FASTA Sequence Generation

This script processes a reference FASTA file and a list of mutations to generate descendant sequences. It can output these sequences into a single multi-FASTA file.

**Primary Functionality: Multi-FASTA File Generation**

This script's main use is to reconstruct potential viral sequences based on a reference and observed mutations.

**Command-Line Usage:**

To generate a multi-FASTA file:

```bash
python prepare_transformer_data.py \
    --reference_fasta path/to/your/reference.fasta \
    --mutations_file path/to/your/mutations_for_clade.txt \
    --output_multifasta_file path/to/output/directory/descendant_sequences.fasta \
    --mutations_delimiter ":"
```

**Arguments:**

*   `--reference_fasta`: Path to the input reference genome FASTA file. (Default: `dummy_reference.fasta`)
*   `--mutations_file`: Path to the input file containing mutations. Each line should be `identifier<delimiter>mutation_string` (e.g., `node_X:A123G,C456T`). (Default: `dummy_mutations.txt`)
*   `--mutations_delimiter`: The delimiter character or string used between the identifier and the mutation string in the mutations file. (Default: `:`)
*   `--output_multifasta_file`: Path where the output multi-FASTA file containing descendant sequences will be saved. (Default: `output_fasta_sequences/descendant_sequences.fasta`)

The script will create the output directory if it doesn't exist. The output FASTA file will include the original reference sequence (with a header like `filename_reference`) and then each descendant sequence derived from the input mutations, using the identifier from the mutations file as the FASTA header.

**Secondary Functionality: Transformer Data Preparation**

The script also retains its original functionality to prepare tokenized and padded data suitable for training a Transformer model. This includes generating (ancestral, descendant) pairs, tokenizing them using a defined nucleotide vocabulary, padding them to a fixed length, and splitting them into training, validation, and test sets. This data is printed to standard output or can be further processed if the script is imported as a module.

## `proof_of_concept/` Directory

The `proof_of_concept/` subdirectory contains Python scripts related to a Transformer-based model for predicting viral mutations:

*   `transformer_model.py`: Defines the Transformer model architecture.
*   `train_transformer.py`: Script for training the Transformer model.
*   `predict_mutations.py`: Script for using a trained model to make predictions.

**Note:** These Transformer-related scripts were part of an experimental implementation. While the model architecture is defined, there were unresolved issues during the training phase (`ValueError` related to Keras argument passing). They are provided as a proof of concept and may require further debugging and development to be fully operational for training and reliable prediction. The primary, stable functionality offered by this suite is the FASTA sequence generation via `prepare_transformer_data.py`.
