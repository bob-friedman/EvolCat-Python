<a name="top"></a>
# Pipeline: SARS-CoV-2 Genomic Data Processing and Transformer Modeling

## Introduction
This tutorial provides a comprehensive, end-to-end guide for processing SARS-CoV-2 genomic data and implementing a Transformer-based model to predict ancestral sequences. We will cover data collection with standard bioinformatics tools, feature engineering in Python, and building a sequence-to-sequence model with TensorFlow/Keras.

This pipeline demonstrates how to leverage the power of Transformer networks to capture complex evolutionary patterns in genomic data.

### 🧬 Pipeline Overview
Our workflow is broken down into four main phases:

1.  **Data Acquisition & ASR:** Use tools like UShER and TreeTime to process a large phylogenetic tree and generate a dataset of (ancestor, descendant) sequence pairs.
2.  **Data Preprocessing:** Use Python to tokenize the DNA sequences and prepare them for the model.
3.  **Model Implementation:** Build a full Encoder-Decoder Transformer model using TensorFlow and Keras.
4.  **Training & Inference:** Train the model on the prepared data and use it to predict ancestral sequences.

### Prerequisites
<details>
<summary><b>Click to view required tools and knowledge</b></summary>

#### Tools:
*   **Conda:** For environment management.
*   **UShER:** For working with large phylogenetic trees (`conda install -c bioconda usher`).
*   **NCBI Datasets CLI:** For downloading reference genomes.
*   **BCFTools:** Utilities for VCF files (`conda install -c bioconda bcftools`).
*   **Python 3.8+** with the following libraries:
    *   `pandas`
    *   `numpy`
    *   `tensorflow`
    *   `biopython`
    *   `pysam`
    *   `dendropy`
    *   `bgzip`

#### Conceptual Knowledge:
*   Basic understanding of genomic data (FASTA, VCF).
*   Familiarity with phylogenetic concepts (ancestral reconstruction, mutations).
*   Basic knowledge of machine learning concepts (Transformers, sequence-to-sequence models).
*   Basic comfort with the command-line.
</details>

### Expected Outcomes
By the end of this tutorial, you will have:
*   A workflow for collecting and processing SARS-CoV-2 data.
*   A trained Transformer model in TensorFlow capable of predicting ancestral sequences.
*   A framework you can adapt for your own research questions in viral genomics.

---

## Phase 1: Data Acquisition & Ancestral Pair Generation
This phase focuses on using command-line bioinformatics tools to go from public databases to a clean, structured dataset of ancestor-descendant sequence pairs.

### Step 1.1: Download Core UShER Data
The UShER project provides a continuously updated global phylogenetic tree for SARS-CoV-2. We will download the core files.
```bash
# Navigate to: https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/
# Download the latest tree, VCF, and metadata files.
wget https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.all.masked.pb.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.metadata.tsv.gz

# Decompress for easier processing by some tools
gunzip public-latest.all.masked.pb.gz
gunzip public-latest.metadata.tsv.gz
```

### Step 1.2: Download the Reference Genome
The Wuhan-Hu-1 genome is the standard reference for SARS-CoV-2.
```bash
# Download reference genome using NCBI Datasets CLI
ncbi-datasets-cli download genome accession NC_045512.2 --include genome --filename SC2_ref.zip
unzip SC2_ref.zip
mv ncbi_dataset/data/NC_045512.2/NC_045512.2.fna ./SC2_ref.fasta
rm -rf ncbi_dataset SC2_ref.zip README.md # Clean up
```

### Step 1.3: Extract a Sub-Clade of Interest
Working with the full >8 million sample tree is computationally prohibitive. We will use `matUtils` (part of UShER) to extract a specific, smaller clade to work with. Here, we use `XBB.2.3` as an example.

```bash
# Extract all samples belonging to clade 'XBB.2.3' into a new VCF file
matUtils extract -i public-latest.all.masked.pb -c 'XBB.2.3' -v xbb23.vcf

# Compress the new VCF with bgzip and index it with bcftools for downstream tools
bgzip --force xbb23.vcf
bcftools index --tbi xbb23.vcf.gz
```
This gives us a manageable VCF file (`xbb23.vcf.gz`) containing only the mutations relevant to our clade of interest.

### Step 1.4: Perform Ancestral State Reconstruction (ASR)
Now we perform ASR on our clade's data to infer the sequences of ancestral nodes. The most direct way is to use `TreeTime`.

```bash
# This is a conceptual command. TreeTime requires a tree and an alignment.
# You would first convert your VCF (xbb23.vcf.gz) and reference (SC2_ref.fasta)
# into a multiple sequence alignment, and extract the corresponding subtree.

# A simplified conceptual workflow:
# 1. Create alignment from VCF: `bcftools consensus ...`
# 2. Extract subtree: `matUtils extract --prune ...`
# 3. Run TreeTime
treetime ancestral --aln clade_alignment.fasta \
                   --tree clade_subtree.newick \
                   --outdir treetime_asr_output
```
This step is complex. For this tutorial, we will assume this step has been run, and we have the key output files: `annotated_tree.nexus` and `ancestral_sequences.fasta`.

### Step 1.5: Generate Ancestor-Descendant Pairs
With the ASR results, we can now create our final training data: a table of (ancestor, descendant) sequence pairs. We use `DendroPy` to traverse the annotated tree and extract these pairs.

<details>
<summary><b>Click to view Python script for generating pairs</b></summary>

```python
import dendropy
import pandas as pd

# --- Configuration ---
TREE_FILE = 'treetime_asr_output/annotated_tree.nexus'
ANCESTRAL_SEQS_FILE = 'treetime_asr_output/ancestral_sequences.fasta'
OUTPUT_PAIRS_FILE = 'ancestor_descendant_pairs.tsv'

print(f"Loading tree from {TREE_FILE}...")
tree = dendropy.Tree.get(path=TREE_FILE, schema="nexus")

print(f"Loading sequences from {ANCESTRAL_SEQS_FILE}...")
seq_matrix = dendropy.DnaCharacterMatrix.get(path=ANCESTRAL_SEQS_FILE, schema="fasta")
seq_dict = {taxon.label: str(seq_matrix[taxon]) for taxon in seq_matrix.taxon_namespace}

print("Generating ancestor-descendant pairs...")
pair_data = []
for node in tree.postorder_node_iter():
    # An ancestor node is any node that is not a leaf
    if not node.is_leaf():
        ancestor_id = node.label
        if ancestor_id not in seq_dict:
            continue
        
        ancestor_seq = seq_dict[ancestor_id]

        # Iterate through its direct children (descendants)
        for child_node in node.child_nodes():
            descendant_id = child_node.label
            if descendant_id not in seq_dict:
                continue
            
            descendant_seq = seq_dict[descendant_id]
            
            # Add the pair to our list
            pair_data.append({
                "ancestor_id": ancestor_id,
                "descendant_id": descendant_id,
                "ancestor_seq": ancestor_seq,
                "descendant_seq": descendant_seq
            })

# Create and save a DataFrame
df_pairs = pd.DataFrame(pair_data)
# Filter out pairs where ancestor and descendant sequences are identical
df_filtered = df_pairs[df_pairs['ancestor_seq'] != df_pairs['descendant_seq']].copy()

df_filtered.to_csv(OUTPUT_PAIRS_FILE, sep='\t', index=False)
print(f"Processing complete. Wrote {len(df_filtered)} non-identical pairs to {OUTPUT_PAIRS_FILE}.")
```
</details>

---

## Phase 2: Data Preprocessing for the Model
Now we take our `ancestor_descendant_pairs.tsv` file and prepare it for the Transformer. This involves tokenization and padding.

<details>
<summary><b>Click to view Python script for data preprocessing</b></summary>

```python
import pandas as pd
from tensorflow import keras

# --- Configuration ---
PAIRS_FILE = 'ancestor_descendant_pairs.tsv'
MAX_SEQ_LENGTH = 30000 # SARS-CoV-2 is ~29.9k bp. Pad to a consistent length.

# --- Load Data ---
df = pd.read_csv(PAIRS_FILE, sep='\t')
# For demonstration, use a smaller subset. Remove this for a full run.
df = df.sample(n=min(10000, len(df)), random_state=42)

# --- Create Vocabulary and Tokenizer ---
# Find all unique characters in both ancestor and descendant sequences
vocab_chars = set("".join(df['ancestor_seq']) + "".join(df['descendant_seq']))
vocab = sorted(list(vocab_chars))

# Add special tokens
special_tokens = ["[PAD]", "[START]", "[END]"]
vocab = special_tokens + vocab
VOCAB_SIZE = len(vocab)
char_to_token = {char: i for i, char in enumerate(vocab)}

print(f"Vocabulary size: {VOCAB_SIZE}")
print(f"Vocabulary: {vocab}")

# --- Tokenization and Padding Function ---
def tokenize_and_pad(sequences, is_target=False):
    tokenized_list = []
    for seq in sequences:
        tokens = [char_to_token.get(char, char_to_token["[PAD]"]) for char in seq]
        if is_target:
            tokens = [char_to_token["[START]"]] + tokens + [char_to_token["[END]"]]
        tokenized_list.append(tokens)
    
    return keras.preprocessing.sequence.pad_sequences(
        tokenized_list, maxlen=MAX_SEQ_LENGTH, padding="post", value=char_to_token["[PAD]"]
    )

# --- Apply Preprocessing ---
# Encoder input: descendant sequences
encoder_input_data = tokenize_and_pad(df["descendant_seq"].values, is_target=False)
# Decoder input/output: ancestor sequences
decoder_target_data = tokenize_and_pad(df["ancestor_seq"].values, is_target=True)

# Prepare inputs for "teacher forcing"
decoder_input_data = decoder_target_data[:, :-1]
decoder_output_data = decoder_target_data[:, 1:]

print(f"\nShape of encoder input vectors: {encoder_input_data.shape}")
print(f"Shape of decoder input vectors: {decoder_input_data.shape}")
print(f"Shape of decoder output vectors: {decoder_output_data.shape}")
```
</details>

---

## Phase 3: Building the Transformer Model
We will implement a full Encoder-Decoder Transformer model using Keras custom layers.

<details>
<summary><b>Click to view the full TensorFlow/Keras implementation</b></summary>

The full implementation involves creating several custom Keras layers for the key components of the Transformer architecture.

#### Transformer Building Blocks (Custom Layers)
```python
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import numpy as np

# --- Model Hyperparameters ---
EMBED_DIM = 128      # Embedding dimension for each token
NUM_HEADS = 8        # Number of attention heads
FF_DIM = 512         # Hidden layer size in feed forward network
NUM_ENCODER_LAYERS = 3
NUM_DECODER_LAYERS = 3
DROPOUT_RATE = 0.1

# --- Positional Embedding Layer ---
class PositionalEmbedding(layers.Layer):
    def __init__(self, vocab_size, embed_dim, maxlen, **kwargs):
        super().__init__(**kwargs)
        self.token_emb = layers.Embedding(input_dim=vocab_size, output_dim=embed_dim)
        self.pos_emb = self.get_positional_encoding(maxlen, embed_dim)

    def get_positional_encoding(self, maxlen, embed_dim):
        positions = np.arange(maxlen)[:, np.newaxis]
        div_term = np.exp(np.arange(0, embed_dim, 2) * -(np.log(10000.0) / embed_dim))
        pos_encoding = np.zeros((maxlen, embed_dim))
        pos_encoding[:, 0::2] = np.sin(positions * div_term)
        pos_encoding[:, 1::2] = np.cos(positions * div_term)
        return tf.cast(pos_encoding[np.newaxis, ...], dtype=tf.float32)

    def call(self, x):
        maxlen = tf.shape(x)[-1]
        x = self.token_emb(x)
        return x + self.pos_emb[:, :maxlen, :]

# --- Transformer Encoder Layer ---
class TransformerEncoder(layers.Layer):
    def __init__(self, embed_dim, ff_dim, num_heads, **kwargs):
        super().__init__(**kwargs)
        self.attn = layers.MultiHeadAttention(num_heads=num_heads, key_dim=embed_dim)
        self.ffn = keras.Sequential([layers.Dense(ff_dim, activation="relu"), layers.Dense(embed_dim)])
        self.layernorm1 = layers.LayerNormalization(epsilon=1e-6)
        self.layernorm2 = layers.LayerNormalization(epsilon=1e-6)
        self.dropout1 = layers.Dropout(DROPOUT_RATE)
        self.dropout2 = layers.Dropout(DROPOUT_RATE)

    def call(self, inputs, training=False):
        attn_output = self.attn(inputs, inputs)
        attn_output = self.dropout1(attn_output, training=training)
        out1 = self.layernorm1(inputs + attn_output)
        ffn_output = self.ffn(out1)
        ffn_output = self.dropout2(ffn_output, training=training)
        return self.layernorm2(out1 + ffn_output)

# --- Transformer Decoder Layer ---
class TransformerDecoder(layers.Layer):
    def __init__(self, embed_dim, ff_dim, num_heads, **kwargs):
        super().__init__(**kwargs)
        self.self_attn = layers.MultiHeadAttention(num_heads=num_heads, key_dim=embed_dim)
        self.cross_attn = layers.MultiHeadAttention(num_heads=num_heads, key_dim=embed_dim)
        self.ffn = keras.Sequential([layers.Dense(ff_dim, activation="relu"), layers.Dense(embed_dim)])
        self.layernorm1 = layers.LayerNormalization(epsilon=1e-6)
        self.layernorm2 = layers.LayerNormalization(epsilon=1e-6)
        self.layernorm3 = layers.LayerNormalization(epsilon=1e-6)
        self.dropout1 = layers.Dropout(DROPOUT_RATE)
        self.dropout2 = layers.Dropout(DROPOUT_RATE)
        self.dropout3 = layers.Dropout(DROPOUT_RATE)

    def call(self, inputs, encoder_outputs, training=False):
        causal_mask = self.get_causal_attention_mask(inputs)
        self_attn_output = self.self_attn(query=inputs, value=inputs, key=inputs, attention_mask=causal_mask)
        self_attn_output = self.dropout1(self_attn_output, training=training)
        out1 = self.layernorm1(inputs + self_attn_output)

        cross_attn_output = self.cross_attn(query=out1, value=encoder_outputs, key=encoder_outputs)
        cross_attn_output = self.dropout2(cross_attn_output, training=training)
        out2 = self.layernorm2(out1 + cross_attn_output)

        ffn_output = self.ffn(out2)
        ffn_output = self.dropout3(ffn_output, training=training)
        return self.layernorm3(out2 + ffn_output)
        
    def get_causal_attention_mask(self, inputs):
        input_shape = tf.shape(inputs)
        batch_size, sequence_length = input_shape[0], input_shape[1]
        i = tf.range(sequence_length)[:, tf.newaxis]
        j = tf.range(sequence_length)
        return tf.cast(i >= j, dtype="int32")
```

#### Build the Full Transformer
```python
# --- Define Inputs ---
encoder_inputs = keras.Input(shape=(None,), dtype="int32", name="descendant_sequence")
decoder_inputs = keras.Input(shape=(None,), dtype="int32", name="ancestor_sequence_input")

# --- Encoder Path ---
encoder_embedding = PositionalEmbedding(VOCAB_SIZE, EMBED_DIM, MAX_SEQ_LENGTH)(encoder_inputs)
x = encoder_embedding
for _ in range(NUM_ENCODER_LAYERS):
    x = TransformerEncoder(EMBED_DIM, FF_DIM, NUM_HEADS)(x)
encoder_outputs = x

# --- Decoder Path ---
decoder_embedding = PositionalEmbedding(VOCAB_SIZE, EMBED_DIM, MAX_SEQ_LENGTH)(decoder_inputs)
x = decoder_embedding
for _ in range(NUM_DECODER_LAYERS):
    x = TransformerDecoder(EMBED_DIM, FF_DIM, NUM_HEADS)(x, encoder_outputs)

# --- Final Output Layer ---
output_logits = layers.Dense(VOCAB_SIZE, name="logits")(x)

# --- Create the Model ---
transformer = keras.Model([encoder_inputs, decoder_inputs], output_logits, name="viral_transformer")
transformer.summary()
```
</details>

---

## Phase 4: Training & Inference
With the model built and data prepared, we can now train it and use it for predictions.

### Training the Model
```python
# Create a tf.data.Dataset for efficient training
dataset = tf.data.Dataset.from_tensor_slices(
    ((encoder_input_data, decoder_input_data), decoder_output_data)
)
dataset = dataset.batch(BATCH_SIZE).shuffle(buffer_size=1024).prefetch(tf.data.AUTOTUNE)

# Compile the model
transformer.compile(
    optimizer=keras.optimizers.Adam(learning_rate=1e-4),
    loss=keras.losses.SparseCategoricalCrossentropy(from_logits=True),
    metrics=[keras.metrics.SparseCategoricalAccuracy()],
)

# Train the model (a GPU is highly recommended)
EPOCHS = 10
history = transformer.fit(dataset, epochs=EPOCHS)
```

### Inference (Predicting Ancestors)
Inference requires an auto-regressive decoding loop where we predict one token at a time.

<details>
<summary><b>Click to view Python code for inference</b></summary>

```python
# A map to convert token IDs back to characters
token_to_char = {i: char for char, i in char_to_token.items()}

def decode_sequence(input_sentence):
    # Tokenize the input descendant sequence
    tokenized_input = tokenize_and_pad([input_sentence])[0]
    
    # Initialize the decoder's input with the [START] token
    decoded_sentence = [char_to_token["[START]"]]
    
    for i in range(MAX_SEQ_LENGTH):
        # Prepare inputs for the model
        encoder_input = tf.constant([tokenized_input], dtype=tf.int32)
        decoder_input = tf.constant([decoded_sentence], dtype=tf.int32)
        
        # Get the model's prediction (logits for the next token)
        predictions = transformer([encoder_input, decoder_input], training=False)
        
        # Sample the token with the highest probability (greedy decoding)
        next_token_id = tf.argmax(predictions[0, i, :]).numpy()
        decoded_sentence.append(next_token_id)
        
        # Stop if the model predicts the [END] token
        if next_token_id == char_to_token["[END]"]:
            break
            
    # Convert the token sequence back to characters, removing special tokens
    return "".join(token_to_char.get(token, "?") for token in decoded_sentence[1:-1])

# --- Test with an example ---
test_idx = 0
input_descendant = df["descendant_seq"].iloc[test_idx]
true_ancestor = df["ancestor_seq"].iloc[test_idx]

predicted_ancestor = decode_sequence(input_descendant)

print("-" * 50)
print(f"INPUT DESCENDANT: {input_descendant[:80]}...")
print(f"TRUE ANCESTOR:    {true_ancestor[:80]}...")
print(f"PREDICTED ANCESTOR: {predicted_ancestor[:80]}...")
print("-" * 50)
```
</details>

---
## Conclusion
This tutorial outlined a comprehensive pipeline for applying a Transformer model to a phylogenetic task. We covered data acquisition with bioinformatics tools, sequence processing, model building, and inference. This framework provides a powerful starting point for leveraging deep learning to gain deeper insights into viral evolution.

## Further Reading
*   **UShER & matUtils:** [https://usher-wiki.readthedocs.io/](https://usher-wiki.readthedocs.io/)
*   **The Viral Chase (Paper):** [https://doi.org/10.20944/preprints202506.0456.v1](https://doi.org/10.20944/preprints202506.0456.v1)
*   **TensorFlow Transformer Tutorial:** [https://www.tensorflow.org/text/tutorials/transformer](https://www.tensorflow.org/text/tutorials/transformer)
*   **Attention Is All You Need:** Vaswani, A., et al. (2017). [https://arxiv.org/abs/1706.03762](https://arxiv.org/abs/1706.03762)

## Credits
Jules AI and Gemini Pro helped develop this documentation page and the code.
