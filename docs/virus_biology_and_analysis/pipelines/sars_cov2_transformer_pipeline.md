# SARS-CoV-2 Genomic Data Processing and Transformer-Based Modeling Pipeline

## Introduction

This tutorial provides a comprehensive guide to processing SARS-CoV-2 genomic data and implementing a Transformer-based model for phylogenetic analysis. We will cover data collection, feature engineering, the conceptual framework behind the model, and a Python implementation using TensorFlow. This pipeline aims to predict ancestral sequences from descendant sequences, leveraging the power of Transformer networks to capture complex patterns in genomic data.

## Prerequisites

### Tools:
*   **Conda:** For environment management.
*   **UShER:** For working with large phylogenetic trees and VCF data. (Installation: `conda install -c bioconda usher`)
*   **NCBI Datasets CLI:** For downloading reference genomes. (Installation: [NCBI Datasets CLI Installation Guide](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/))
*   **Python 3.8+**
*   **Pandas:** For data manipulation. (`pip install pandas`)
*   **NumPy:** For numerical operations. (`pip install numpy`)
*   **TensorFlow 2.x:** For building and training the Transformer model. (`pip install tensorflow`)
*   **BioPython:** For sequence manipulation. (`pip install biopython`)
*   **PySam:** For reading VCF files. (`pip install pysam`)
*   **DendroPy:** For phylogenetic computations, tree manipulation, and ASR. (`pip install dendropy`)
*   **(Optional) Nextflow:** For pipeline orchestration.

### Conceptual Knowledge:
*   Basic understanding of genomic data (FASTA, VCF formats).
*   Familiarity with phylogenetic concepts (ancestral sequence reconstruction, phylogenetic trees, mutations).
*   Basic knowledge of machine learning concepts (Transformers, sequence-to-sequence models, embeddings, attention mechanisms).
*   Understanding of basic Linux/command-line operations.

## Expected Outcomes

By completing this tutorial, you will:
*   Understand the workflow for collecting and processing SARS-CoV-2 genomic data, specifically focusing on creating ancestor-descendant sequence pairs.
*   Learn how to perform feature engineering using Pandas to prepare data for a machine learning model.
*   Gain insights into the conceptual framework of using Transformer models for predicting ancestral genomic sequences.
*   Have a trained Transformer model in TensorFlow capable of predicting ancestral sequences from descendant sequences.
*   Be able to adapt this pipeline for your own research questions in viral genomics or similar sequence-based tasks.
*   Understand how to prepare and preprocess data for a nucleotide-based Transformer model.

## DATA COLLECTION

This section outlines the steps to collect and preprocess SARS-CoV-2 data to create a dataset of (ancestor, descendant) sequence pairs.

### 1. Database Construction with UShER

The initial step involves constructing or updating a comprehensive SARS-CoV-2 database using UShER. This database typically includes a global phylogenetic tree (`public-latest.all.masked.pb.gz`) and corresponding mutation-annotated VCF files.

```bash
# Download the latest UShER tree and VCF (example paths)
# These files are usually obtained from sources like https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/
wget https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.all.masked.pb.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.all.masked.vcf.gz

# It's good practice to decompress the VCF for easier processing by some tools,
# though UShER itself can handle .vcf.gz
gunzip public-latest.all.masked.vcf.gz
# Now you have public-latest.all.masked.vcf
```
*Note: Ensure UShER is properly installed and configured. The paths to the database files should be specified.*

### 2. Downloading Reference Genome

We need a reference genome for alignment and variant calling. The Wuhan-Hu-1 genome (NC_045512.2) is commonly used for SARS-CoV-2.

```bash
# Download reference genome and annotation using NCBI Datasets CLI
datasets download genome accession NC_045512.2 --include genome --filename SC2_ref.zip
unzip SC2_ref.zip
mv ncbi_dataset/data/NC_045512.2/NC_045512.2.fna ./SC2_ref.fasta
rm -rf ncbi_dataset SC2_ref.zip README.md # Clean up
```

### 3. Extracting Ancestor-Descendant Sequence Pairs using UShER and `matUtils`

`matUtils extract` (a tool packaged with UShER) is crucial for generating the data needed. We will use it to extract information about ancestor-descendant pairs along with their mutations.

```bash
# Ensure UShER and matUtils are in your PATH
# Path to your UShER protobuf tree
TREE_PATH="./public-latest.all.masked.pb.gz"
# Path to your VCF file (unzipped)
VCF_PATH="./public-latest.all.masked.vcf"
# Output directory for extracted data
OUTPUT_DIR="./usher_extracted_data"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Use matUtils to extract ancestor-descendant information
# -t: input tree
# -v: input VCF
# -d: output directory
# -A: output ancestor sequences in FASTA format (ancestors.fa)
# -D: output descendant sequences in FASTA format (descendants.fa)
# -M: output mutations for each descendant relative to its ancestor (mutation_paths.txt)
# -N: output node names for descendants (descendant_nodes.txt)
# -P: output node names for ancestors (ancestor_nodes.txt)
matUtils extract -t ${TREE_PATH} -v ${VCF_PATH} -d ${OUTPUT_DIR} -A -D -M -N -P
```

This command will generate several files in `${OUTPUT_DIR}`:
*   `ancestors.fa`: FASTA file of ancestral sequences.
*   `descendants.fa`: FASTA file of descendant sequences.
*   `mutation_paths.txt`: Text file detailing mutations between each ancestor-descendant pair.
*   `ancestor_nodes.txt`: List of ancestor node identifiers.
*   `descendant_nodes.txt`: List of descendant node identifiers.

Each line in `ancestor_nodes.txt`, `descendant_nodes.txt`, `ancestors.fa` (sequence headers), and `descendants.fa` (sequence headers) corresponds to an ancestor-descendant pair. `mutation_paths.txt` also aligns with these pairs.

## Feature Engineering with Pandas

Now, we'll process the files generated by `matUtils` into a structured Pandas DataFrame. This DataFrame will be the input for our Transformer model.

### Gemini: Script to Create DataFrame for Ancestor-Descendant Pairs

This Python script reads the FASTA files and creates a Pandas DataFrame.

```python
from Bio import SeqIO
import pandas as pd
import numpy as np

def load_sequences_to_df(ancestor_fasta_path, descendant_fasta_path,
                         ancestor_nodes_path, descendant_nodes_path):
    """
    Loads ancestor and descendant sequences and their node IDs into a Pandas DataFrame.

    Args:
        ancestor_fasta_path (str): Path to the FASTA file of ancestral sequences.
        descendant_fasta_path (str): Path to the FASTA file of descendant sequences.
        ancestor_nodes_path (str): Path to the text file with ancestor node IDs.
        descendant_nodes_path (str): Path to the text file with descendant node IDs.

    Returns:
        pandas.DataFrame: DataFrame with columns ['ancestor_id', 'descendant_id',
                                                 'ancestor_seq', 'descendant_seq'].
    """
    ancestor_seqs = {record.id: str(record.seq) for record in SeqIO.parse(ancestor_fasta_path, "fasta")}
    descendant_seqs = {record.id: str(record.seq) for record in SeqIO.parse(descendant_fasta_path, "fasta")}

    with open(ancestor_nodes_path, 'r') as f_anc_nodes, open(descendant_nodes_path, 'r') as f_desc_nodes:
        ancestor_node_ids = [line.strip() for line in f_anc_nodes]
        descendant_node_ids = [line.strip() for line in f_desc_nodes]

    if len(ancestor_node_ids) != len(descendant_node_ids):
        raise ValueError("Mismatch in the number of ancestor and descendant node IDs.")
    if len(ancestor_seqs) != len(descendant_seqs) or len(ancestor_seqs) != len(ancestor_node_ids):
        print(f"Warning: Mismatch in lengths - anc_seqs: {len(ancestor_seqs)}, "
              f"desc_seqs: {len(descendant_seqs)}, node_ids: {len(ancestor_node_ids)}")
        # This might happen if matUtils outputs fewer sequences than node IDs due to some filtering.
        # We will proceed by pairing based on the provided node ID lists and looking up sequences.
        # Sequences not found will result in NaN and can be filtered later.

    data = []
    for i in range(len(ancestor_node_ids)):
        anc_id = ancestor_node_ids[i]
        desc_id = descendant_node_ids[i]

        data.append({
            'ancestor_id': anc_id,
            'descendant_id': desc_id,
            'ancestor_seq': ancestor_seqs.get(anc_id),
            'descendant_seq': descendant_seqs.get(desc_id)
        })

    df = pd.DataFrame(data)

    # Drop rows where sequences might be missing if IDs didn't match FASTA headers
    df.dropna(subset=['ancestor_seq', 'descendant_seq'], inplace=True)

    return df

# --- Main script execution ---
ancestor_fasta_file = "./usher_extracted_data/ancestors.fa"
descendant_fasta_file = "./usher_extracted_data/descendants.fa"
ancestor_nodes_file = "./usher_extracted_data/ancestor_nodes.txt"
descendant_nodes_file = "./usher_extracted_data/descendant_nodes.txt"

df_sequences = load_sequences_to_df(ancestor_fasta_file, descendant_fasta_file,
                                    ancestor_nodes_file, descendant_nodes_file)

print("Shape of the DataFrame:", df_sequences.shape)
print(df_sequences.head())

# Save the DataFrame (optional)
df_sequences.to_csv("ancestor_descendant_sequences.csv", index=False)
```

### Note: Identical Ancestor-Descendant Pairs

It's possible that `matUtils extract` outputs pairs where the ancestor and descendant sequences are identical (zero mutations between them). This can happen for various reasons related to tree structure and mutation imputation. For a model learning to predict changes, these identical pairs might not be informative or could even bias the training if they are overly represented.

It's advisable to check for and potentially filter out these identical pairs.

```python
# Assuming df_sequences is your DataFrame from the previous step
print(f"Original DataFrame shape: {df_sequences.shape}")

# Filter out identical sequences
df_sequences_filtered = df_sequences[df_sequences['ancestor_seq'] != df_sequences['descendant_seq']].copy()

print(f"Filtered DataFrame shape (non-identical pairs): {df_sequences_filtered.shape}")
num_filtered = len(df_sequences) - len(df_sequences_filtered)
print(f"Number of identical ancestor-descendant pairs removed: {num_filtered}")

# Use the filtered DataFrame for subsequent steps
df_model_input = df_sequences_filtered
# df_model_input.to_csv("ancestor_descendant_sequences_filtered.csv", index=False) # Optional save
```
*Gemini: The decision to remove identical pairs depends on the specific goals. If the model is meant to predict *any* ancestor, including those identical to the input, they could be kept. However, for focusing on mutational changes, filtering is common.*

## Conceptual framework for remaining steps

The core idea is to train a sequence-to-sequence Transformer model.
*   **Input (Encoder):** The descendant sequence.
*   **Output (Decoder):** The ancestor sequence.

The model learns to "translate" a descendant sequence back to its immediate ancestor.

### Data Preprocessing for Transformer
1.  **Tokenization:** Convert nucleotide sequences (A, C, G, T, N, and potentially other characters like '-') into integer tokens.
2.  **Padding:** Ensure all sequences in a batch have the same length by padding shorter sequences.
3.  **Special Tokens:** Add `[START]` and `[END]` tokens to target sequences (ancestor sequences) for the decoder.

## Implement the Transformer based method

This section provides a Python/TensorFlow implementation of the Transformer model.
*Gemini: This implementation is based on the TensorFlow tutorial "Transformer model for language understanding" ([https://www.tensorflow.org/text/tutorials/transformer](https://www.tensorflow.org/text/tutorials/transformer)), adapted for nucleotide sequences.*

### Setup

```python
import tensorflow as tf
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

# Assuming df_model_input is loaded and contains 'descendant_seq' and 'ancestor_seq'
# For demonstration, let's create a dummy DataFrame if not loaded
if 'df_model_input' not in locals():
    print("Creating dummy df_model_input for demonstration.")
    data = {
        'descendant_seq': ['AGCTNNA', 'CGTANCG', 'GTNCA--', 'NNAACGT'],
        'ancestor_seq': ['AGCCTNA', 'CGTACCG', 'GTNCAGG', 'NNATCGT']
    }
    df_model_input = pd.DataFrame(data)
    # Ensure sequences are of reasonable length for a real scenario
    # This dummy data is very short.
    # Example of making them longer:
    # df_model_input['descendant_seq'] = df_model_input['descendant_seq'] * 500
    # df_model_input['ancestor_seq'] = df_model_input['ancestor_seq'] * 500


# --- Parameters ---
VOCAB = sorted(list(set("".join(df_model_input['descendant_seq']) + "".join(df_model_input['ancestor_seq']))))
# Ensure padding token '-' is present, and special tokens for start/end
if '-' not in VOCAB: VOCAB.append('-') # Padding
# VOCAB should ideally be fixed based on expected characters e.g. ['A', 'C', 'G', 'T', 'N', '-']
# For this example, we derive it from data, but add common ones.
standard_chars = ['A', 'C', 'G', 'T', 'N', '-']
for char in standard_chars:
    if char not in VOCAB: VOCAB.append(char)
VOCAB = sorted(list(set(VOCAB)))


# Add START and END tokens to vocabulary
START_TOKEN = '[START]'
END_TOKEN = '[END]'
if START_TOKEN not in VOCAB: VOCAB.append(START_TOKEN)
if END_TOKEN not in VOCAB: VOCAB.append(END_TOKEN)
VOCAB_SIZE = len(VOCAB)

# Create tokenizers (char to int and int to char)
token_to_id = {char: i for i, char in enumerate(VOCAB)}
id_to_token = {i: char for i, char in enumerate(VOCAB)}

# Determine max sequence length from the data
# This should be done carefully on the actual dataset.
# For SARS-CoV-2, this is typically around 29903, but matUtils output might be shorter if only core genome used.
# Add 2 for [START] and [END] tokens for the target sequence.
MAX_LENGTH = max(max(df_model_input['descendant_seq'].apply(len)),
                 max(df_model_input['ancestor_seq'].apply(len))) + 2

print(f"Vocabulary: {VOCAB}")
print(f"Vocabulary Size: {VOCAB_SIZE}")
print(f"Max sequence length (with START/END): {MAX_LENGTH}")

# --- Tokenization and Padding Functions ---
def tokenize_sequence(sequence, tokenizer, max_length, add_start_end_tokens=False):
    tokens = [tokenizer[START_TOKEN]] if add_start_end_tokens else []
    tokens.extend([tokenizer.get(char, tokenizer['-']) for char in sequence]) # Default to padding char if unknown
    if add_start_end_tokens:
        tokens.append(tokenizer[END_TOKEN])

    # Pad
    padding_needed = max_length - len(tokens)
    tokens.extend([tokenizer['-']] * padding_needed)
    return tokens

def detokenize_sequence(token_ids, id_tokenizer):
    return "".join([id_tokenizer.get(tid, '') for tid in token_ids if tid != token_to_id[START_TOKEN] and tid != token_to_id[END_TOKEN] and tid != token_to_id['-']])

# --- Prepare data ---
# Input features (descendant sequences)
input_data = np.array([tokenize_sequence(seq, token_to_id, MAX_LENGTH, add_start_end_tokens=False)
                       for seq in df_model_input['descendant_seq']])

# Target features (ancestor sequences) - with START and END tokens
target_data = np.array([tokenize_sequence(seq, token_to_id, MAX_LENGTH, add_start_end_tokens=True)
                        for seq in df_model_input['ancestor_seq']])


# Split data
# Using a smaller test size for quick demo; adjust as needed
X_train, X_val, y_train, y_val = train_test_split(input_data, target_data, test_size=0.2, random_state=42)

print(f"X_train shape: {X_train.shape}, y_train shape: {y_train.shape}")
print(f"X_val shape: {X_val.shape}, y_val shape: {y_val.shape}")

# Create TensorFlow Datasets
BUFFER_SIZE = 10000 # Adjust based on dataset size
BATCH_SIZE = 64   # Adjust based on memory

train_dataset = tf.data.Dataset.from_tensor_slices((X_train, y_train))
train_dataset = train_dataset.shuffle(BUFFER_SIZE).batch(BATCH_SIZE).prefetch(tf.data.AUTOTUNE)

val_dataset = tf.data.Dataset.from_tensor_slices((X_val, y_val))
val_dataset = val_dataset.batch(BATCH_SIZE).prefetch(tf.data.AUTOTUNE)
```

### Transformer Model Components (adapted from TensorFlow tutorial)

```python
# --- Positional Encoding ---
def positional_encoding(length, depth):
    depth = depth / 2
    positions = np.arange(length)[:, np.newaxis]  # (seq, 1)
    depths = np.arange(depth)[np.newaxis, :] / depth  # (1, depth)
    angle_rates = 1 / (10000**depths)  # (1, depth)
    angle_rads = positions * angle_rates  # (pos, depth)
    pos_encoding = np.concatenate([np.sin(angle_rads), np.cos(angle_rads)], axis=-1)
    return tf.cast(pos_encoding, dtype=tf.float32)

class PositionalEmbedding(tf.keras.layers.Layer):
    def __init__(self, vocab_size, d_model, max_length):
        super().__init__()
        self.d_model = d_model
        self.embedding = tf.keras.layers.Embedding(vocab_size, d_model, mask_zero=True)
        self.pos_encoding = positional_encoding(length=max_length, depth=d_model)

    def compute_mask(self, *args, **kwargs):
        return self.embedding.compute_mask(*args, **kwargs)

    def call(self, x):
        length = tf.shape(x)[1]
        x = self.embedding(x)
        # This factor sets the relative scale of the embedding and positonal_encoding.
        x *= tf.math.sqrt(tf.cast(self.d_model, tf.float32))
        x = x + self.pos_encoding[tf.newaxis, :length, :]
        return x

# --- Base Attention ---
class BaseAttention(tf.keras.layers.Layer):
    def __init__(self, **kwargs):
        super().__init__()
        self.mha = tf.keras.layers.MultiHeadAttention(**kwargs)
        self.layernorm = tf.keras.layers.LayerNormalization()
        self.add = tf.keras.layers.Add()

# --- Cross Attention (for Decoder) ---
class CrossAttention(BaseAttention):
    def call(self, x, context):
        attn_output, attn_scores = self.mha(
            query=x,
            key=context,
            value=context,
            return_attention_scores=True)
        self.last_attn_scores = attn_scores # Cache attention scores
        x = self.add([x, attn_output])
        x = self.layernorm(x)
        return x

# --- Global Self Attention (for Encoder) ---
class GlobalSelfAttention(BaseAttention):
    def call(self, x):
        attn_output = self.mha(query=x, value=x, key=x)
        x = self.add([x, attn_output])
        x = self.layernorm(x)
        return x

# --- Causal Self Attention (for Decoder) ---
class CausalSelfAttention(BaseAttention):
    def call(self, x):
        attn_output = self.mha(query=x, value=x, key=x, use_causal_mask=True)
        x = self.add([x, attn_output])
        x = self.layernorm(x)
        return x

# --- Feed Forward Network ---
class FeedForward(tf.keras.layers.Layer):
    def __init__(self, d_model, dff, dropout_rate=0.1):
        super().__init__()
        self.seq = tf.keras.Sequential([
            tf.keras.layers.Dense(dff, activation='relu'),
            tf.keras.layers.Dense(d_model),
            tf.keras.layers.Dropout(dropout_rate)
        ])
        self.add = tf.keras.layers.Add()
        self.layernorm = tf.keras.layers.LayerNormalization()

    def call(self, x):
        seq_out = self.seq(x)
        x = self.add([x, seq_out])
        x = self.layernorm(x)
        return x

# --- Encoder Layer ---
class EncoderLayer(tf.keras.layers.Layer):
    def __init__(self,*, d_model, num_heads, dff, dropout_rate=0.1):
        super().__init__()
        self.self_attention = GlobalSelfAttention(num_heads=num_heads, key_dim=d_model, dropout=dropout_rate)
        self.ffn = FeedForward(d_model, dff)

    def call(self, x):
        x = self.self_attention(x)
        x = self.ffn(x)
        return x

# --- Encoder ---
class Encoder(tf.keras.layers.Layer):
    def __init__(self, *, num_layers, d_model, num_heads, dff, vocab_size, max_length, dropout_rate=0.1):
        super().__init__()
        self.d_model = d_model
        self.num_layers = num_layers
        self.pos_embedding = PositionalEmbedding(vocab_size=vocab_size, d_model=d_model, max_length=max_length)
        self.enc_layers = [EncoderLayer(d_model=d_model, num_heads=num_heads, dff=dff, dropout_rate=dropout_rate)
                           for _ in range(num_layers)]
        self.dropout = tf.keras.layers.Dropout(dropout_rate)

    def call(self, x):
        x = self.pos_embedding(x) # (batch, seq_len, d_model)
        x = self.dropout(x)
        for i in range(self.num_layers):
            x = self.enc_layers[i](x)
        return x # (batch, seq_len, d_model)

# --- Decoder Layer ---
class DecoderLayer(tf.keras.layers.Layer):
    def __init__(self, *, d_model, num_heads, dff, dropout_rate=0.1):
        super(DecoderLayer, self).__init__()
        self.causal_self_attention = CausalSelfAttention(num_heads=num_heads, key_dim=d_model, dropout=dropout_rate)
        self.cross_attention = CrossAttention(num_heads=num_heads, key_dim=d_model, dropout=dropout_rate)
        self.ffn = FeedForward(d_model, dff)

    def call(self, x, context):
        x = self.causal_self_attention(x=x)
        x = self.cross_attention(x=x, context=context)
        self.last_attn_scores = self.cross_attention.last_attn_scores # Cache attn scores
        x = self.ffn(x) # Shape `(batch_size, seq_len, d_model)`.
        return x

# --- Decoder ---
class Decoder(tf.keras.layers.Layer):
    def __init__(self, *, num_layers, d_model, num_heads, dff, vocab_size, max_length, dropout_rate=0.1):
        super(Decoder, self).__init__()
        self.d_model = d_model
        self.num_layers = num_layers
        self.pos_embedding = PositionalEmbedding(vocab_size=vocab_size, d_model=d_model, max_length=max_length)
        self.dropout = tf.keras.layers.Dropout(dropout_rate)
        self.dec_layers = [DecoderLayer(d_model=d_model, num_heads=num_heads, dff=dff, dropout_rate=dropout_rate)
                           for _ in range(num_layers)]
        self.last_attn_scores = None

    def call(self, x, context):
        x = self.pos_embedding(x) # (batch, seq_len, d_model)
        x = self.dropout(x)
        for i in range(self.num_layers):
            x = self.dec_layers[i](x, context)
        self.last_attn_scores = self.dec_layers[-1].last_attn_scores
        return x # (batch, seq_len, d_model)

# --- Transformer Model ---
class Transformer(tf.keras.Model):
    def __init__(self, *, num_layers, d_model, num_heads, dff,
                 input_vocab_size, target_vocab_size, max_length, dropout_rate=0.1):
        super().__init__()
        self.encoder = Encoder(num_layers=num_layers, d_model=d_model,
                               num_heads=num_heads, dff=dff,
                               vocab_size=input_vocab_size, max_length=max_length, dropout_rate=dropout_rate)

        self.decoder = Decoder(num_layers=num_layers, d_model=d_model,
                               num_heads=num_heads, dff=dff,
                               vocab_size=target_vocab_size, max_length=max_length, dropout_rate=dropout_rate)

        self.final_layer = tf.keras.layers.Dense(target_vocab_size)

    def call(self, inputs):
        context, x = inputs # context is descendant_seq (encoder input), x is ancestor_seq (decoder input)
        context = self.encoder(context)  # (batch_size, context_len, d_model)
        x = self.decoder(x, context)  # (batch_size, target_len, d_model)
        logits = self.final_layer(x)  # (batch_size, target_len, target_vocab_size)
        try:
            del logits._keras_mask # Remove mask before softmax
        except AttributeError:
            pass
        return logits
```

### Build and Train the Model

```python
# --- Hyperparameters ---
# These are small for demonstration. Increase for real training.
NUM_LAYERS = 2 # Originally 6 in many Transformer papers
D_MODEL = 128   # Originally 512
DFF = 256       # Originally 2048
NUM_HEADS = 4   # Originally 8
DROPOUT_RATE = 0.1

# Instantiate the Transformer model
transformer = Transformer(
    num_layers=NUM_LAYERS,
    d_model=D_MODEL,
    num_heads=NUM_HEADS,
    dff=DFF,
    input_vocab_size=VOCAB_SIZE, # Assuming same vocab for input and target for simplicity
    target_vocab_size=VOCAB_SIZE,
    max_length=MAX_LENGTH,
    dropout_rate=DROPOUT_RATE)

# --- Loss and Optimizer ---
def masked_loss(label, pred):
    mask = label != token_to_id['-'] # Assuming '-' is padding
    loss_object = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True, reduction='none')
    loss = loss_object(label, pred)
    mask = tf.cast(mask, dtype=loss.dtype)
    loss *= mask
    loss = tf.reduce_sum(loss)/tf.reduce_sum(mask)
    return loss

def masked_accuracy(label, pred):
    pred = tf.argmax(pred, axis=2)
    label = tf.cast(label, pred.dtype)
    match = label == pred
    mask = label != token_to_id['-']
    match = match & mask
    match = tf.cast(match, dtype=tf.float32)
    mask = tf.cast(mask, dtype=tf.float32)
    return tf.reduce_sum(match)/tf.reduce_sum(mask)

# Learning rate schedule (optional, but often helpful)
class CustomSchedule(tf.keras.optimizers.schedules.LearningRateSchedule):
    def __init__(self, d_model, warmup_steps=4000):
        super().__init__()
        self.d_model = tf.cast(d_model, tf.float32)
        self.warmup_steps = warmup_steps

    def __call__(self, step):
        step = tf.cast(step, dtype=tf.float32)
        arg1 = tf.math.rsqrt(step)
        arg2 = step * (self.warmup_steps**-1.5)
        return tf.math.rsqrt(self.d_model) * tf.math.minimum(arg1, arg2)

# learning_rate = CustomSchedule(D_MODEL)
# optimizer = tf.keras.optimizers.Adam(learning_rate, beta_1=0.9, beta_2=0.98, epsilon=1e-9)
optimizer = tf.keras.optimizers.Adam(learning_rate=0.001) # Simpler optimizer for demo

transformer.compile(optimizer=optimizer, loss=masked_loss, metrics=[masked_accuracy])

# --- Training ---
EPOCHS = 5 # Increase for real training (e.g., 20-100+)
# This will restart the kernel if you are using a notebook and run out of memory.
# Reduce BATCH_SIZE or MAX_LENGTH if that happens.
# For very long sequences like SARS-CoV-2 (29903bp), this model will require significant GPU memory.
# Consider strategies like:
# 1. Using a smaller MAX_LENGTH by focusing on variable regions if applicable.
# 2. Gradient accumulation if batch size needs to be very small.
# 3. More powerful GPUs.
# 4. Model parallelism (more complex).

print("\nStarting training...")
try:
    history = transformer.fit(
        train_dataset,
        epochs=EPOCHS,
        validation_data=val_dataset
    )
    print("Training completed.")
    # Save model weights
    # transformer.save_weights('./transformer_ancestor_predictor_weights') # HDF5 format
    # print("Model weights saved.")

except Exception as e:
    print(f"An error occurred during training: {e}")
    print("Ensure MAX_LENGTH is appropriate for your data and available memory.")
    print("For SARS-CoV-2 full genome (~30k bp), D_MODEL=128, NUM_LAYERS=2, BATCH_SIZE=16 can already be demanding.")

# Example of how to load weights:
# loaded_transformer = Transformer(...) # Instantiate with same parameters
# loaded_transformer.load_weights('./transformer_ancestor_predictor_weights')
# print("Model weights loaded.")
```

### Inference / Prediction

```python
class AncestorPredictor(tf.Module):
    def __init__(self, transformer_model, token_to_id_map, id_to_token_map, max_length):
        self.transformer = transformer_model
        self.token_to_id = token_to_id_map
        self.id_to_token = id_to_token_map
        self.max_length = max_length
        self.start_token_id = self.token_to_id['[START]']
        self.end_token_id = self.token_to_id['[END]']
        self.padding_token_id = self.token_to_id['-']


    def __call__(self, descendant_sequence_str, temperature=0.0): # temperature > 0 for stochasticity
        # Tokenize the input descendant sequence
        descendant_tokens = tokenize_sequence(descendant_sequence_str, self.token_to_id, self.max_length, add_start_end_tokens=False)
        encoder_input = tf.constant([descendant_tokens], dtype=tf.int64) # Add batch dimension

        # Start the decoder input with the [START] token
        output_array = tf.TensorArray(dtype=tf.int64, size=0, dynamic_size=True)
        output_array = output_array.write(0, self.start_token_id)

        for i in tf.range(self.max_length -1 ): # -1 because [START] is already added
            output = tf.transpose(output_array.stack()) # (1, current_length)

            # Pad the current output to MAX_LENGTH for the decoder's positional embedding
            # This is a simplified way; for efficiency, the positional embedding in decoder
            # could be sliced based on current length if model is modified.
            current_length = tf.shape(output)[1]
            padding_needed = self.max_length - current_length

            # Create padding tensor
            paddings = [[0,0], [0, padding_needed]] # Pad only on the sequence dimension at the end
            output_padded = tf.pad(output, paddings, "CONSTANT", constant_values=self.padding_token_id)

            predictions = self.transformer([encoder_input, output_padded], training=False) # (batch, target_len, vocab_size)

            # Select the last token from the seq_len dimension
            predictions = predictions[:, current_length-1, :]  # (batch_size, vocab_size)

            if temperature == 0.0:
                predicted_id = tf.argmax(predictions, axis=-1) # Greedy search
            else:
                # Sample (helps avoid repetitive loops sometimes)
                # Ensure predictions are not all -inf or NaN before applying temperature
                safe_predictions = tf.where(tf.math.is_finite(predictions), predictions, tf.zeros_like(predictions) - 1e9)
                predicted_id = tf.random.categorical(safe_predictions / temperature, num_samples=1)[0] # (1,)

            predicted_id = tf.cast(predicted_id[0], tf.int64) # Get scalar

            if tf.equal(predicted_id, self.end_token_id):
                break # Stop if [END] token is predicted

            output_array = output_array.write(i + 1, predicted_id)

        output = tf.transpose(output_array.stack())
        predicted_ancestor_str = detokenize_sequence(output[0].numpy(), self.id_to_token)

        return predicted_ancestor_str


# --- (Untested) Gemini: Example usage of the predictor ---
# Ensure the transformer model used here is the trained one.
# If you just ran training, 'transformer' variable holds the trained model.
# If loading from saved weights, instantiate and load weights first.

# predictor = AncestorPredictor(transformer, token_to_id, id_to_token, MAX_LENGTH)

# # Take a sample from validation set or provide a new sequence
# if X_val.shape[0] > 0:
#     sample_descendant_tokens = X_val[0] # This is tokenized and padded
#     # Detokenize for string input to predictor
#     sample_descendant_str = detokenize_sequence(sample_descendant_tokens, id_to_token)
#
#     print(f"\nInput Descendant (string): {sample_descendant_str}")
#
#     # Get the true ancestor for comparison
#     true_ancestor_tokens = y_val[0] # This is tokenized, padded, with START/END
#     true_ancestor_str = detokenize_sequence(true_ancestor_tokens, id_to_token)
#     print(f"True Ancestor (string):    {true_ancestor_str}")
#
#     predicted_ancestor = predictor(sample_descendant_str)
#     print(f"Predicted Ancestor (string): {predicted_ancestor}")
# else:
#     print("Validation set is empty, cannot demonstrate prediction with a sample.")
#     # Example with a dummy sequence if X_val is empty
#     # dummy_desc_seq = "ACGTN" # Replace with a more realistic sequence
#     # print(f"Input Descendant (string): {dummy_desc_seq}")
#     # predicted_ancestor = predictor(dummy_desc_seq)
#     # print(f"Predicted Ancestor (string): {predicted_ancestor}")

```
*Gemini: The inference loop above is a basic greedy decoder. More advanced techniques like beam search can improve results but are more complex to implement. The padding in the inference loop for `output_padded` is a simplification; a more optimal approach would involve adjusting positional encodings or model inputs dynamically.*

## Conclusion

This tutorial outlined a pipeline for processing SARS-CoV-2 genomic data using UShER and `matUtils` to extract ancestor-descendant sequence pairs, and then applying a Transformer model built with TensorFlow to predict ancestral sequences from descendant sequences. By following these steps, researchers can leverage advanced machine learning techniques to gain deeper insights into viral evolution and ancestral states.

The provided Transformer implementation is a starting point. Further refinements could include:
*   More extensive hyperparameter tuning.
*   Using more sophisticated learning rate schedules.
*   Implementing beam search for inference.
*   Training on much larger datasets derived from comprehensive UShER trees.
*   Adapting the model for different types of genomic input or prediction tasks (e.g., predicting specific mutations, recombination events).
*   Careful consideration of sequence length, GPU memory, and batch sizes, especially for full-length SARS-CoV-2 genomes.

## Further Reading and Resources

*   **UShER & matUtils:** [https://usher-wiki.readthedocs.io/](https://usher-wiki.readthedocs.io/)
*   **NCBI Datasets:** [https://www.ncbi.nlm.nih.gov/datasets/](https://www.ncbi.nlm.nih.gov/datasets/)
*   **TensorFlow Transformer Tutorial:** [https://www.tensorflow.org/text/tutorials/transformer](https://www.tensorflow.org/text/tutorials/transformer)
*   **Effective Transformers (Book/Blog by Lilian Weng):** [https://lilianweng.github.io/posts/2023-01-27-the-transformer-family-v2/](https://lilianweng.github.io/posts/2023-01-27-the-transformer-family-v2/) (For deeper understanding of Transformer variants)
*   **Attention Is All You Need (Original Transformer Paper):** Vaswani, A., et al. (2017). [https://arxiv.org/abs/1706.03762](https://arxiv.org/abs/1706.03762)
*   **DendroPy Documentation:** [https://dendropy.org/](https://dendropy.org/)

This comprehensive pipeline provides a foundation for applying deep learning to complex phylogenetic questions.
