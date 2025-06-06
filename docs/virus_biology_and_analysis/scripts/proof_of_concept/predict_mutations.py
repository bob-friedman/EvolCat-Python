import tensorflow as tf
import numpy as np
import argparse # For command-line arguments
import os

# Import from your other scripts
from prepare_transformer_data import tokenize_sequence, NUCLEOTIDE_VOCAB, REV_NUCLEOTIDE_VOCAB, parse_fasta_file, pad_sequences
from transformer_model import Transformer

# --- 1. Global Hyperparameters/Constants (MUST MATCH TRAINED MODEL) ---
# These should ideally be loaded from a config file or saved with the model
# For now, define them here. Ensure they match the 'train_transformer.py' settings used for the saved model.

# Model parameters (example values, replace with actual trained model's params)
NUM_LAYERS = 2 # Must match trained model
D_MODEL = 64   # Must match trained model
NUM_HEADS = 4  # Must match trained model
DFF = 128      # Must match trained model
DROPOUT_RATE = 0.1 # Only relevant if model uses it during inference, but good to define

# Vocabulary sizes (should be consistent)
INPUT_VOCAB_SIZE = len(NUCLEOTIDE_VOCAB)
TARGET_VOCAB_SIZE = len(NUCLEOTIDE_VOCAB)

# Max positional encoding lengths used during training
# This also defines the fixed sequence length the model expects after padding.
MAX_SEQUENCE_LENGTH = 60 # Example: Must match pe_input/pe_target from training, and padding length

# Path to the saved model weights (HDF5 file)
DEFAULT_MODEL_WEIGHTS_PATH = "./saved_model_train_script/final_transformer_weights.h5"

# --- 2. Load Trained Model ---
def load_trained_model(num_layers, d_model, num_heads, dff,
                       input_vocab_size, target_vocab_size,
                       max_pos_enc_input, max_pos_enc_target,
                       dropout_rate, model_weights_path):
    """
    Instantiates the Transformer model and loads trained weights.
    """
    print(f"Loading model with the following parameters:")
    print(f"  NUM_LAYERS: {num_layers}, D_MODEL: {d_model}, NUM_HEADS: {num_heads}, DFF: {dff}")
    print(f"  INPUT_VOCAB_SIZE: {input_vocab_size}, TARGET_VOCAB_SIZE: {target_vocab_size}")
    print(f"  MAX_POS_ENC_INPUT: {max_pos_enc_input}, MAX_POS_ENC_TARGET: {max_pos_enc_target}")
    print(f"  DROPOUT_RATE: {dropout_rate}")

    transformer = Transformer(
        num_layers=num_layers,
        d_model=d_model,
        num_heads=num_heads,
        dff=dff,
        input_vocab_size=input_vocab_size,
        target_vocab_size=target_vocab_size,
        pe_input=max_pos_enc_input,
        pe_target=max_pos_enc_target,
        rate=dropout_rate  # Dropout rate for consistency, though usually not active in inference
    )

    # Build the model by calling it with dummy data that has the expected shape
    # This is often necessary before loading weights for subclassed models.
    # Input shape: (batch_size, seq_len)
    dummy_encoder_input = tf.zeros((1, max_pos_enc_input), dtype=tf.int32)
    dummy_decoder_input = tf.zeros((1, 1), dtype=tf.int32) # Start with a single token for decoder

    try:
        _ = transformer((dummy_encoder_input, dummy_decoder_input), training=False)
        print("Model built successfully with dummy data.")
    except Exception as e:
        print(f"Error building model with dummy data: {e}")
        # Depending on the error, you might still try to load weights or exit

    print(f"Loading weights from: {model_weights_path}")
    transformer.load_weights(model_weights_path)
    print("Model weights loaded successfully.")
    return transformer

# --- 3. Predict Sequence (Greedy Decoding) ---
def predict_sequence(model, input_sequence_str, max_output_length=MAX_SEQUENCE_LENGTH):
    """
    Predicts the descendant sequence from an ancestral sequence string.
    """
    # Tokenize and pad the input sequence
    input_tokens = tokenize_sequence(input_sequence_str, NUCLEOTIDE_VOCAB)
    # Pad to MAX_SEQUENCE_LENGTH (which is pe_input for the encoder)
    # The pad_sequences function expects a list of token lists (pairs)
    # We adapt it here for a single sequence.
    # pad_sequences returns list of (padded_anc, padded_desc, mut_str), we only need padded_anc

    # Simpler padding for a single sequence:
    padded_input_tokens = input_tokens[:max_output_length] + \
                          [NUCLEOTIDE_VOCAB['<PAD>']] * max(0, max_output_length - len(input_tokens))

    encoder_input = tf.expand_dims(padded_input_tokens, axis=0) # Shape: (1, max_output_length)

    # Decoder input starts with a "start token".
    # Using <PAD> as a proxy start token. In a more refined setup, a dedicated <START> token is better.
    # The target sequence for the decoder during training starts with this token.
    decoder_input_tokens = [NUCLEOTIDE_VOCAB['<PAD>']]
    decoder_input = tf.expand_dims(decoder_input_tokens, 0) # Shape: (1, 1)

    print(f"Starting prediction loop for max output length: {max_output_length}")
    for i in range(max_output_length -1): # -1 because we already have one start token
        # predictions shape: (batch_size, tar_seq_len, target_vocab_size)
        predictions, attention_weights = model((encoder_input, decoder_input), training=False)

        # Select the last token from the seq_len dimension
        predicted_token_logits = predictions[:, -1:, :]  # Shape: (batch_size, 1, target_vocab_size)

        # Apply argmax to get the most likely next token ID
        predicted_id = tf.argmax(predicted_token_logits, axis=-1, output_type=tf.int32) # Shape: (batch_size, 1)

        # Concatenate the predicted_id to the decoder input (which is used as next input)
        decoder_input = tf.concat([decoder_input, predicted_id], axis=-1)

        # Optional: if predicted_id is an end token, break
        # if predicted_id.numpy()[0][0] == NUCLEOTIDE_VOCAB.get('<END>', -1): # Assuming <END> token
        #     print(f"End token predicted at step {i+1}.")
        #     break

    # Remove the initial start token from the result
    predicted_tokens = decoder_input.numpy()[0][1:]

    # Convert token IDs back to nucleotide characters
    predicted_nucleotides = [REV_NUCLEOTIDE_VOCAB.get(token_id, '?') for token_id in predicted_tokens]

    return "".join(predicted_nucleotides)

# --- 4. Identify Mutations ---
def identify_mutations(original_seq_str, predicted_seq_str):
    """
    Compares original and predicted sequences to identify mutations.
    Returns a list of mutation strings (e.g., "A1G").
    """
    mutations = []
    min_len = min(len(original_seq_str), len(predicted_seq_str))

    for i in range(min_len):
        if original_seq_str[i] != predicted_seq_str[i]:
            mutation = f"{original_seq_str[i]}{i+1}{predicted_seq_str[i]}"
            mutations.append(mutation)

    # Handle cases where predicted sequence is shorter/longer (insertions/deletions)
    # For this model, output length is fixed, so this part might not be triggered
    # unless the <PAD> token is interpreted as an early stop.
    if len(predicted_seq_str) > min_len:
        for i in range(min_len, len(predicted_seq_str)):
            mutations.append(f"_{i+1}{predicted_seq_str[i]}") # Insertion relative to original
    elif len(original_seq_str) > min_len:
         for i in range(min_len, len(original_seq_str)):
            mutations.append(f"{original_seq_str[i]}{i+1}_") # Deletion relative to original

    return mutations

# --- 5. Main Function ---
def main():
    parser = argparse.ArgumentParser(description="Predict viral mutations using a Transformer model.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--input_sequence", type=str, help="Ancestral nucleotide sequence string.")
    group.add_argument("--input_fasta", type=str, help="Path to a FASTA file containing the ancestral sequence (first entry used).")

    parser.add_argument("--weights_path", type=str, default=DEFAULT_MODEL_WEIGHTS_PATH,
                        help=f"Path to the trained model weights HDF5 file (default: {DEFAULT_MODEL_WEIGHTS_PATH}).")
    # Add arguments for other hyperparameters if they need to be highly configurable at runtime
    # For now, we use the globally defined ones.
    # parser.add_argument("--d_model", type=int, default=D_MODEL, help="D_MODEL for the transformer.")
    # ... other params

    args = parser.parse_args()

    print("--- Prediction Script Initialized ---")
    print(f"Using NUCLEOTIDE_VOCAB with {len(NUCLEOTIDE_VOCAB)} entries.")

    input_sequence_str = ""
    if args.input_sequence:
        input_sequence_str = args.input_sequence.upper()
    elif args.input_fasta:
        print(f"Loading sequence from FASTA: {args.input_fasta}")
        # Assuming parse_fasta_file returns a single string (the first sequence)
        # and converts to uppercase.
        input_sequence_str = parse_fasta_file(args.input_fasta)
        if not input_sequence_str:
            print(f"Error: Could not read sequence from FASTA file {args.input_fasta}. Exiting.")
            return

    if not input_sequence_str:
        print("Error: No input sequence provided. Exiting.")
        return

    # Ensure model weights path exists
    if not os.path.exists(args.weights_path):
        print(f"Error: Model weights file not found at {args.weights_path}.")
        print("Please ensure the path is correct and the model has been trained and saved.")
        # Create a dummy weights file for the script to run without error for structure check
        # This is NOT a valid model, just for testing the script's flow.
        print(f"Attempting to create a placeholder dummy weights file at {args.weights_path} for structure testing...")
        try:
            # Need to instantiate a model to save dummy weights
            dummy_model_for_weights = Transformer(
                NUM_LAYERS, D_MODEL, NUM_HEADS, DFF,
                INPUT_VOCAB_SIZE, TARGET_VOCAB_SIZE,
                MAX_SEQUENCE_LENGTH, MAX_SEQUENCE_LENGTH, DROPOUT_RATE
            )
            # Build it with a dummy call
            dummy_enc_inp = tf.zeros((1, MAX_SEQUENCE_LENGTH), dtype=tf.int32)
            dummy_dec_inp = tf.zeros((1, 1), dtype=tf.int32)
            _ = dummy_model_for_weights((dummy_enc_inp, dummy_dec_inp), training=False)

            os.makedirs(os.path.dirname(args.weights_path), exist_ok=True)
            dummy_model_for_weights.save_weights(args.weights_path)
            print(f"Placeholder dummy weights file created at {args.weights_path}.")
            print("WARNING: This is NOT a trained model. Predictions will be random.")
        except Exception as e:
            print(f"Could not create placeholder weights file: {e}. Please provide valid weights.")
            return


    # Load the model
    print("\n--- Loading Trained Model ---")
    try:
        model = load_trained_model(
            num_layers=NUM_LAYERS, d_model=D_MODEL, num_heads=NUM_HEADS, dff=DFF,
            input_vocab_size=INPUT_VOCAB_SIZE, target_vocab_size=TARGET_VOCAB_SIZE,
            max_pos_enc_input=MAX_SEQUENCE_LENGTH, # Consistent with padding
            max_pos_enc_target=MAX_SEQUENCE_LENGTH, # Consistent with padding
            dropout_rate=DROPOUT_RATE,
            model_weights_path=args.weights_path
        )
    except Exception as e:
        print(f"Error loading the model: {e}")
        return

    # Perform prediction
    print("\n--- Performing Prediction ---")
    print(f"Input Sequence (first 100 chars): {input_sequence_str[:100]}")
    predicted_sequence_str = predict_sequence(model, input_sequence_str, max_output_length=MAX_SEQUENCE_LENGTH)
    print(f"Predicted Sequence (first 100 chars): {predicted_sequence_str[:100]}")

    # Identify and print mutations
    print("\n--- Identifying Mutations ---")
    mutations = identify_mutations(input_sequence_str, predicted_sequence_str)
    if mutations:
        print("Identified mutations:")
        for mut in mutations:
            print(f"  {mut}")
    else:
        print("No mutations identified (or predicted sequence is identical up to min length).")

    print("\n--- Prediction Script Finished ---")

if __name__ == "__main__":
    print("predict_mutations.py script structure defined.")
    # To run from command line:
    # python predict_mutations.py --input_sequence "AGCT..."
    # or
    # python predict_mutations.py --input_fasta "path/to/your/sequence.fasta"
    # (Potentially add --weights_path if not using default)
    main()

```
