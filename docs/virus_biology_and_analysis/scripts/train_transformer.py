import tensorflow as tf
import numpy as np
import time
import os # For saving models

# Import from your other scripts
from prepare_transformer_data import (
    parse_fasta_file,
    parse_mutation_file,
    create_input_output_pairs,
    NUCLEOTIDE_VOCAB,
    pad_sequences,
    split_data
)
from transformer_model import Transformer

# --- 1. Hyperparameters ---
# Model parameters
NUM_LAYERS = 2 # Reduced for faster dummy training
D_MODEL = 64   # Reduced for faster dummy training
NUM_HEADS = 4  # Reduced
DFF = 128      # Reduced
INPUT_VOCAB_SIZE = len(NUCLEOTIDE_VOCAB)
TARGET_VOCAB_SIZE = len(NUCLEOTIDE_VOCAB)
MAX_POSITIONAL_ENCODING_INPUT = 60 # Adjusted to dummy data length
MAX_POSITIONAL_ENCODING_TARGET = 60 # Adjusted to dummy data length
DROPOUT_RATE = 0.1

# Training parameters
EPOCHS = 3 # Reduced for faster dummy training
BATCH_SIZE = 2 # Reduced for dummy data
WARMUP_STEPS = 100 # Reduced for dummy data

# Data and Model Paths
REFERENCE_FASTA_PATH = "reference_genome_train.fasta" # Unique name for this script's dummy
MUTATIONS_FILE_PATH = "mutations_file_train.txt"   # Unique name
CHECKPOINT_PATH = "./checkpoints_train_script/train" # Checkpoint directory
MODEL_SAVE_PATH = "./saved_model_train_script/final_transformer" # Final model save path

# --- Custom Learning Rate Schedule ---
class CustomSchedule(tf.keras.optimizers.schedules.LearningRateSchedule):
    def __init__(self, d_model, warmup_steps=4000):
        super(CustomSchedule, self).__init__()
        self.d_model = d_model
        self.d_model = tf.cast(self.d_model, tf.float32)
        self.warmup_steps = warmup_steps

    def __call__(self, step):
        step = tf.cast(step, tf.float32)
        arg1 = tf.math.rsqrt(step)
        arg2 = step * (self.warmup_steps**-1.5)
        # Add small epsilon to prevent division by zero if step is 0
        return tf.math.rsqrt(self.d_model) * tf.math.minimum(arg1, arg2 + 1e-9)


def run_training():
    print("Starting training process...")
    os.makedirs(os.path.dirname(CHECKPOINT_PATH), exist_ok=True)
    os.makedirs(os.path.dirname(MODEL_SAVE_PATH), exist_ok=True)

    # --- 2. Data Loading and Preparation ---
    print("\n--- Loading and Preparing Data ---")
    seq_len_dummy = MAX_POSITIONAL_ENCODING_INPUT # Use consistent length
    if not os.path.exists(REFERENCE_FASTA_PATH) or os.path.getsize(REFERENCE_FASTA_PATH) == 0:
        print(f"Warning: Dummy {REFERENCE_FASTA_PATH} created.")
        with open(REFERENCE_FASTA_PATH, "w") as f:
            f.write(f">dummy_ref\n{'A'*seq_len_dummy}\n")
    if not os.path.exists(MUTATIONS_FILE_PATH) or os.path.getsize(MUTATIONS_FILE_PATH) == 0:
        print(f"Warning: Dummy {MUTATIONS_FILE_PATH} created.")
        with open(MUTATIONS_FILE_PATH, "w") as f:
            # Create 10 mutations for 8/1/1 split
            for i in range(1, 11):
                f.write(f"node{i}: A{i%seq_len_dummy+1}T,C{(i+5)%seq_len_dummy+1}G\n")

    reference_genome = parse_fasta_file(REFERENCE_FASTA_PATH)
    if not reference_genome:
        print(f"Error: Could not load reference genome. Exiting.")
        return

    mutation_strings = parse_mutation_file(MUTATIONS_FILE_PATH)
    if not mutation_strings:
        mutation_strings = [""] # Ensure pipeline runs with at least one "no mutation" data point
        print(f"Warning: No valid mutations loaded. Using reference as only sequence.")

    all_tokenized_pairs = create_input_output_pairs(reference_genome, mutation_strings, NUCLEOTIDE_VOCAB)

    max_len_for_padding = MAX_POSITIONAL_ENCODING_INPUT
    all_padded_pairs = pad_sequences(all_tokenized_pairs, max_len_for_padding, NUCLEOTIDE_VOCAB['<PAD>'])

    if not all_padded_pairs:
        print("No data after padding. Exiting.")
        return

    train_pairs, val_pairs, test_pairs = split_data(all_padded_pairs, train_ratio=0.8, val_ratio=0.1, random_seed=42)
    print(f"Data split: {len(train_pairs)} train, {len(val_pairs)} validation, {len(test_pairs)} test pairs.")

    def to_tf_dataset(pairs_list, batch_size, max_len):
        if not pairs_list:
            return tf.data.Dataset.from_tensor_slices(
                ((tf.zeros((0, max_len), dtype=tf.int32), tf.zeros((0, max_len - 1), dtype=tf.int32)),
                 tf.zeros((0, max_len - 1), dtype=tf.int32))
            ).batch(batch_size).prefetch(tf.data.AUTOTUNE)

        inputs_enc = np.array([p[0] for p in pairs_list], dtype=np.int32)
        inputs_dec = np.array([p[1][:-1] for p in pairs_list], dtype=np.int32)
        targets = np.array([p[1][1:] for p in pairs_list], dtype=np.int32)

        dataset = tf.data.Dataset.from_tensor_slices(((inputs_enc, inputs_dec), targets))
        dataset = dataset.shuffle(len(pairs_list), reshuffle_each_iteration=True)
        dataset = dataset.batch(batch_size)
        dataset = dataset.prefetch(tf.data.AUTOTUNE)
        return dataset

    train_dataset = to_tf_dataset(train_pairs, BATCH_SIZE, max_len_for_padding)
    val_dataset = to_tf_dataset(val_pairs, BATCH_SIZE, max_len_for_padding)
    print("TensorFlow Datasets created.")

    # --- 3. Model Instantiation ---
    print("\n--- Initializing Model ---")
    transformer = Transformer(
        num_layers=NUM_LAYERS, d_model=D_MODEL, num_heads=NUM_HEADS, dff=DFF,
        input_vocab_size=INPUT_VOCAB_SIZE, target_vocab_size=TARGET_VOCAB_SIZE,
        pe_input=max_len_for_padding, pe_target=max_len_for_padding, rate=DROPOUT_RATE)
    print("Transformer model instantiated.")

    # --- 4. Loss, Optimizer, Metrics, Checkpointing ---
    print("\n--- Setting up Training Components ---")
    loss_object = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True, reduction='none')

    def loss_function(real, pred):
        mask = tf.math.logical_not(tf.math.equal(real, NUCLEOTIDE_VOCAB['<PAD>']))
        loss_ = loss_object(real, pred)
        mask = tf.cast(mask, dtype=loss_.dtype)
        loss_ *= mask
        return tf.reduce_sum(loss_) / tf.reduce_sum(mask)

    learning_rate_custom = CustomSchedule(D_MODEL, warmup_steps=WARMUP_STEPS)
    optimizer = tf.keras.optimizers.Adam(learning_rate_custom, beta_1=0.9, beta_2=0.98, epsilon=1e-9)

    def accuracy_function(real, pred):
        mask = tf.math.logical_not(tf.math.equal(real, NUCLEOTIDE_VOCAB['<PAD>']))
        match = tf.equal(real, tf.cast(tf.argmax(pred, axis=2), dtype=real.dtype))
        match = tf.math.logical_and(match, mask)
        match_cast = tf.cast(match, dtype=tf.float32)
        mask_cast = tf.cast(mask, dtype=tf.float32)
        sum_match = tf.reduce_sum(match_cast)
        sum_mask = tf.reduce_sum(mask_cast)
        return tf.where(sum_mask == 0, 0.0, sum_match / sum_mask) # Avoid division by zero

    epoch_loss_avg = tf.keras.metrics.Mean(name='train_loss')
    epoch_accuracy_avg = tf.keras.metrics.Mean(name='train_accuracy')
    epoch_val_loss_avg = tf.keras.metrics.Mean(name='val_loss')
    epoch_val_accuracy_avg = tf.keras.metrics.Mean(name='val_accuracy')

    ckpt = tf.train.Checkpoint(transformer=transformer, optimizer=optimizer)
    ckpt_manager = tf.train.CheckpointManager(ckpt, CHECKPOINT_PATH, max_to_keep=5)

    if ckpt_manager.latest_checkpoint:
        ckpt.restore(ckpt_manager.latest_checkpoint).expect_partial() # expect_partial for flexibility
        print(f'Latest checkpoint restored from {ckpt_manager.latest_checkpoint}!')
    else:
        print('No checkpoint found. Initializing from scratch.')
    print("Training components (loss, optimizer, metrics, checkpointing) defined.")

    # --- 5. Training & Validation Steps ---
    @tf.function
    def train_step(inputs_tuple, tar_real):
        inp_enc, inp_dec = inputs_tuple
        with tf.GradientTape() as tape:
            predictions, _ = transformer((inp_enc, inp_dec), training=True)
            loss = loss_function(tar_real, predictions)
        gradients = tape.gradient(loss, transformer.trainable_variables)
        optimizer.apply_gradients(zip(gradients, transformer.trainable_variables))
        acc = accuracy_function(tar_real, predictions)
        return loss, acc

    @tf.function
    def validation_step(inputs_tuple, tar_real):
        inp_enc, inp_dec = inputs_tuple
        predictions, _ = transformer((inp_enc, inp_dec), training=False)
        loss = loss_function(tar_real, predictions)
        acc = accuracy_function(tar_real, predictions)
        return loss, acc
    print("Train and validation step functions defined.")

    # --- 6. Training Loop ---
    print("\n--- Starting Training ---")
    if len(train_pairs) == 0:
        print("No training data. Skipping training loop.")
    else:
        for epoch in range(EPOCHS):
            start_time = time.time()
            epoch_loss_avg.reset_states()
            epoch_accuracy_avg.reset_states()
            epoch_val_loss_avg.reset_states()
            epoch_val_accuracy_avg.reset_states()

            print(f"--- Epoch {epoch + 1}/{EPOCHS} ---")
            # Training
            for (batch_num, (inputs_for_model, target_for_loss)) in enumerate(train_dataset):
                batch_loss, batch_accuracy = train_step(inputs_for_model, target_for_loss)
                epoch_loss_avg.update_state(batch_loss)
                epoch_accuracy_avg.update_state(batch_accuracy)
                if batch_num % 1 == 0: # Print more frequently for small dummy data
                    print(f'  Training: Batch {batch_num} Loss {batch_loss.numpy():.4f} Accuracy {batch_accuracy.numpy():.4f}')
            print(f'Epoch {epoch + 1} Avg Training Loss: {epoch_loss_avg.result():.4f}, Avg Training Accuracy: {epoch_accuracy_avg.result():.4f}')

            # Validation
            if len(val_pairs) > 0:
                print("\n  Validation phase...")
                for (batch_num, (inputs_for_model, target_for_loss)) in enumerate(val_dataset):
                    batch_val_loss, batch_val_accuracy = validation_step(inputs_for_model, target_for_loss)
                    epoch_val_loss_avg.update_state(batch_val_loss)
                    epoch_val_accuracy_avg.update_state(batch_val_accuracy)
                print(f'Epoch {epoch + 1} Validation Loss: {epoch_val_loss_avg.result():.4f}, Validation Accuracy: {epoch_val_accuracy_avg.result():.4f}')
            else:
                print("No validation data to evaluate.")

            ckpt_save_path = ckpt_manager.save()
            print(f'Saving checkpoint for epoch {epoch+1} at {ckpt_save_path}')
            print(f'Time taken for epoch: {time.time() - start_time:.2f} secs\n')
    print("Training loop finished.")

    # --- 7. Final Model Saving ---
    print("\n--- Saving Final Model Weights ---")
    # transformer.save_weights(MODEL_SAVE_PATH) # For saving only weights
    # To save the entire model (if it's a subclassed Model and built):
    # try:
    #     # Build the model with a sample input if not already built
    #     if not transformer.built:
    #          # Create a sample input batch (use specs from one of the datasets)
    #         sample_input_batch, _ = next(iter(train_dataset)) # Get one batch
    #         transformer(sample_input_batch, training=False) # Call to build
    #     transformer.save(MODEL_SAVE_PATH) # Saves in SavedModel format
    #     print(f"Full model saved to {MODEL_SAVE_PATH}")
    # except Exception as e:
    #     print(f"Could not save full model: {e}. Saving weights instead.")
    transformer.save_weights(MODEL_SAVE_PATH + "_weights.h5")
    print(f"Model weights saved to {MODEL_SAVE_PATH}_weights.h5")

if __name__ == "__main__":
    print("train_transformer.py script starting...")
    run_training()
    print("\nTraining script execution finished.")
```
