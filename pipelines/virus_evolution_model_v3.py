# ==============================================================================
# Part 1: Environment Setup
# ==============================================================================
# This section should be in its own cell. After it runs, the kernel will restart.
# You must then run the subsequent cells.

# Installation of Conda for Google Colab
!pip install -q condacolab
import condacolab
condacolab.install()

# After this cell runs and the kernel restarts, proceed to the next cell.

# ==============================================================================
# Part 2: Install Dependencies
# ==============================================================================
# This cell runs after the kernel restart from the Conda installation.

import os
from google.colab import drive

# Mount Google Drive
drive.mount('/content/drive')
os.chdir('/content/drive/MyDrive')
print(f"Current working directory: {os.getcwd()}")

# Install all necessary bioinformatics tools into the Conda environment
# This single command is more efficient than multiple separate installs.
# Verify installations with !mafft --version; !iqtree --version; !mmseqs -h; !matUtils --version
!conda install -y -c bioconda -c conda-forge usher mafft iqtree mmseqs2 entrez-direct pysam bcftools

# ==============================================================================
# Part 3: Data Collection and S-Gene Retrieval
# ==============================================================================
import os
import sys
import pandas as pd
import gc

# --- Helper Function ---
def ensure_file_exists(file_path, command_to_run):
    """Checks if a file exists. If not, runs a command to generate it."""
    print(f"Checking for: {os.path.basename(file_path)}...")
    if os.path.isfile(file_path):
        print(f"--> Found. Skipping generation.\n")
    else:
        print(f"--> Not found. Running command:\n    $ {command_to_run}")
        os.system(command_to_run)
        if os.path.isfile(file_path):
            print(f"--> Verified: '{os.path.basename(file_path)}' created successfully.\n")
        else:
            print(f"--> FATAL: Command finished, but file '{os.path.basename(file_path)}' was not found.\n")
            sys.exit(1) # Exit if a critical file is not created

# --- File Definitions ---
BASE_PATH = "/content/drive/MyDrive"
PB_FILE = os.path.join(BASE_PATH, "public-latest.all.masked.pb")
METADATA_FILE_FULL = os.path.join(BASE_PATH, "public-latest.metadata.tsv")
REF_FASTA = os.path.join(BASE_PATH, "NC_045512v2.fa")

CLADE_NAME = "XBB.2.3"
CLADE_VCF_FILE = os.path.join(BASE_PATH, "xbb23.vcf")
CLADE_METADATA_FILE = os.path.join(BASE_PATH, "public-xbb23-latest.metadata.tsv")
ACCESSION_FILE = os.path.join(BASE_PATH, "xbb23_accession_list.txt")
GENE_S_CDS_RAW = os.path.join(BASE_PATH, "gene_sequences_cds_raw.fas")
GENE_S_FASTA = os.path.join(BASE_PATH, "spike_sequences_raw.fas")

# --- Download Core Data ---
ensure_file_exists(PB_FILE, "wget http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.all.masked.pb.gz && gunzip -f public-latest.all.masked.pb.gz")
ensure_file_exists(METADATA_FILE_FULL, "wget http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.metadata.tsv.gz && gunzip -f public-latest.metadata.tsv.gz")
ensure_file_exists(REF_FASTA, "wget https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/chromosomes/NC_045512v2.fa.gz && gunzip -f NC_045512v2.fa.gz")

# --- Extract Clade-Specific Data ---
print(f"Extracting VCF for clade '{CLADE_NAME}'...")
# matUtils is fragile in Colab, so it may not report an error as it resides within this function
# modify names manually below if changing the "File Definitions" above this code block
ensure_file_exists(CLADE_VCF_FILE, f"matUtils extract -i public-latest.all.masked.pb -c 'XBB.2.3' -v xbb23.vcf")

# --- Generate Accession List ---
print("Generating accession list for VCF samples...")
if not os.path.isfile(ACCESSION_FILE):
    # 1. Get sample names directly from the VCF
    !bcftools query -l {CLADE_VCF_FILE} > vcf_samples.txt
    with open("vcf_samples.txt", "r") as f:
        vcf_samples_set = set(line.strip() for line in f)

    # 2. Load only necessary metadata columns for memory efficiency
    print(f"Loading relevant columns from {os.path.basename(METADATA_FILE_FULL)}...")
    metadata_df = pd.read_csv(
        METADATA_FILE_FULL,
        sep='\t',
        usecols=['strain', 'genbank_accession'],
        dtype=str
    )

    # 3. Filter metadata for samples present in our VCF
    print("Filtering metadata for VCF samples...")
    clade_metadata_df = metadata_df[metadata_df['strain'].isin(vcf_samples_set)].copy()
    clade_metadata_df.dropna(subset=['genbank_accession'], inplace=True)
    clade_metadata_df = clade_metadata_df[clade_metadata_df['genbank_accession'] != '']

    accession_list = clade_metadata_df['genbank_accession'].unique()

    # 4. Save the list to a file
    with open(ACCESSION_FILE, 'w') as f:
        for acc in accession_list:
            f.write(f"{acc}\n")

    print(f"--> Successfully created '{os.path.basename(ACCESSION_FILE)}' with {len(accession_list)} unique accession numbers.\n")

    # Clean up memory
    del metadata_df, clade_metadata_df, vcf_samples_set
    gc.collect()
else:
    print(f"--> Found '{os.path.basename(ACCESSION_FILE)}'. Skipping generation.\n")

# Commented out IPython magic to ensure Python compatibility.
# %%shell
# # This is a shell only script (no Python code allowed in this Colab cell) for NCBI sequence retrieval.
# # 'grep' not ideal for parsing multiline fasta files. Implemented instead in 'awk'.
# 
# # === VARIABLES ===
# ACCESSION_FILE="xbb23_accession_list.txt" # Your accession file
# GENE_PATTERN="[gene=S]"
# GENE_S_FASTA="s_gene_sequences.fasta"
# LOG_FILE="bare_metal_pipeline.log"
# SCRIPT_FILE="bare_metal_run.sh"
# 
# # === SCRIPT CREATION ===
# # This script abandons efetch for retrieval and uses curl directly.
# cat << 'EOF' > ${SCRIPT_FILE}
# #!/bin/bash
# set -e
# 
# # --- Get Variables from Outside ---
# ACCESSION_FILE="$1"
# GENE_PATTERN="$2"
# GENE_S_FASTA="$3"
# BATCH_SIZE=100
# 
# echo "--- Preparing for 'bare metal' controlled batch processing using curl ---"
# 
# TEMP_DIR=$(mktemp -d)
# trap 'rm -rf -- "$TEMP_DIR"' EXIT
# rm -f "${GENE_S_FASTA}"
# 
# echo "Splitting ${ACCESSION_FILE} into chunks of ${BATCH_SIZE}..."
# split -l "${BATCH_SIZE}" "${ACCESSION_FILE}" "${TEMP_DIR}/batch_"
# echo "Splitting complete."
# echo ""
# 
# # --- Main Logic: Loop Over Each Batch File ---
# for batch_file in ${TEMP_DIR}/batch_*
# do
#     echo "--- Processing batch: $(basename ${batch_file}) ---"
# 
#     # 1. Read the 100 accessions from the batch file and join them with a comma.
#     id_list=$(paste -sd ',' "${batch_file}")
# 
#     # 2. ABANDON EFETCH. Use curl to build the raw POST request ourselves. There is no hidden logic.
#     curl -s -X POST "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi" \
#          -d "db=nuccore" \
#          -d "rettype=fasta_cds_na" \
#          -d "retmode=text" \
#          -d "id=${id_list}" | \
#     awk -v pattern="${GENE_PATTERN}" '
#         BEGIN { RS = ">"; OFS = "" }
#         NF > 0 {
#             header_end = index($0, "\n")
#             header = substr($0, 1, header_end - 1)
#             sequence = substr($0, header_end + 1)
#             gsub(/\n/, "", sequence)
# 
#             if (index(header, pattern)) {
#                 print ">" header, "\n" sequence
#             }
#         }
#     ' >> "${GENE_S_FASTA}"
# 
#     echo "Batch $(basename ${batch_file}) finished. Sleeping for 1 second."
#     sleep 1
# done
# 
# echo ""
# echo "--- All batches processed. ---"
# 
# 
# # --- Final Report ---
# if [ -s "${GENE_S_FASTA}" ]; then
#     echo "SUCCESS: Pipeline finished."
#     echo "Found $(grep -c "^>" ${GENE_S_FASTA}) sequences in ${GENE_S_FASTA}."
# else
#     echo "ERROR: Pipeline finished, but the output file is empty."
#     exit 1
# fi
# EOF
# 
# # === SCRIPT EXECUTION ===
# bash ${SCRIPT_FILE} "${ACCESSION_FILE}" "${GENE_PATTERN}" "${GENE_S_FASTA}" 2>&1 | tee "${LOG_FILE}"

# ==============================================================================
# Part 4: MSA, Phylogeny, and Ancestral State Reconstruction (ASR)
# ==============================================================================
# Tested with small sample size since IQ-TREE is cpu limited with large datasets

# Use the system package manager to install seqtk
# The -y flag automatically answers "yes" to prompts
!sudo apt-get update && sudo apt-get install -y seqtk

# --- File Definitions for this section ---
GENE_S_FASTA = "s_gene_sequences.fasta"
REPRESENTATIVE_SEQS = "spike_representatives.fas"
REPRESENTATIVE_SEQS2 = "spike_reps_300.fas"
ALIGNED_SEQS = "spike_aligned.fas"
IQTREE_PREFIX = "iqtree_asr_run" # prefix to all iqtree output files

# (Optional) Simplify the FASTA Headers for debugging any compatibility problems.
# May be needed for downstream analysis, however.
# Also, would have to edit REPRESENTATIVE_SEQS to spike_representatives_simplified.fas
# !awk '/^>/ {print ">seq" ++i; next} {print}' spike_representatives.fas > spike_representatives_simplified.fas

# Step 1: Cluster sequences to retieve unique representatives.
# This speeds up alignment by removing 100% identical sequences.
# The `_rep_seq.fasta` file below used for alignment and tree-building.
print("\nClustering sequences to find unique representatives...")
!mmseqs easy-cluster {GENE_S_FASTA} cluster_results /content/tmp_mmseqs --min-seq-id 1.0 -c 0.8 --cov-mode 1
# The key output is `cluster_results_rep_seq.fasta`. Renamed for clarity.
!mv cluster_results_rep_seq.fasta {REPRESENTATIVE_SEQS}

# (Optional) Create a random sample of 300 sequences using seqtk
!seqtk sample -s100 {REPRESENTATIVE_SEQS} 300 > {REPRESENTATIVE_SEQS2}

# Step 2: Create a Multiple Sequence Alignment (MSA) of the representative sequences.
# MAFFT will automatically choose a suitable algorithm (like FFT-NS-2).
print(f"\nAligning representative sequences with MAFFT...")
# Run the command directly in the cell, using the variables, otherwise standard output fails
!mafft --auto {REPRESENTATIVE_SEQS2} > {ALIGNED_SEQS}

# Step 3: IQ-TREE for phylogenetic inference and ancestral state reconstruction.
# -s: input alignment
# -m MFP: ModelFinder Plus automatically selects the best substitution model.
# -B 5000: Ultrafast bootstrap with 5000 replicates to assess branch support.
# -asr: Perform ancestral state reconstruction.
# -pre: Set a prefix for all output files (e.g., iqtree_asr_run.treefile, .log, .state).
# -T AUTO: Use an optimal number of threads.
print(f"\nRunning IQ-TREE for ASR...")
# Refer to other pipeline version for limitations on this step (see "rapid radiation" as a confounding factor).
# If pipeline fails or to debug this step, then revert to the two commands below instead of within function below.
# !iqtree -s spike_aligned_300.fas -m MFP -B 5000 -T AUTO -safe -redo
# !iqtree -s spike_aligned_300.fas -t spike_aligned_300.fas.treefile -asr -undo
ensure_file_exists(
    f"{IQTREE_PREFIX}.state", # check for .state file (ancestral/descendant predictions)
    # parameter 'redo' to overwrite any checkpoints from previous runs (checkpoints allow for a "resume" feature)
    # however, for asr as a separate step, then may use -undo to rely on a checkpoint (see its logs)
    f"iqtree -s {ALIGNED_SEQS} -m MFP -B 5000 -asr -pre {IQTREE_PREFIX} -T AUTO -safe -undo"
)

print("\nIQ-TREE run complete. Key output files:")
!ls -lh {IQTREE_PREFIX}*

# ==============================================================================
# Part 5: The Bridge - Create Ancestor-Descendant Pairs from IQ-TREE Output
# ==============================================================================

# The code reads the tree structure from `iqtree_asr_run.treefile` and the ancestral sequences from
# `iqtree_asr_run.state`, then pairs them up into the `ancestor_descendant_pairs.tsv` file.

import os
from google.colab import drive

# Mount Google Drive
drive.mount('/content/drive')
os.chdir('/content/drive/MyDrive')
print(f"Current working directory: {os.getcwd()}")

!pip install -q ete3 biopython

from ete3 import Tree
from Bio import SeqIO
import pandas as pd

# --- File Definitions ---
IQTREE_PREFIX = "iqtree_asr_run" # prefix to all iqtree output files
TREE_FILE = f"{IQTREE_PREFIX}.treefile"
STATE_FILE = f"{IQTREE_PREFIX}.state" # This is a tab-delimited text file of all node sequences
OUTPUT_PAIRS_FILE = "ancestor_descendant_pairs.tsv"

print(f"\nParsing tree from '{TREE_FILE}' and sequences from '{STATE_FILE}'...")

# 1. Load all sequences (tips and ancestral nodes) into a dictionary.
# The .state file contains sequences for both terminal leaves and internal nodes.
#### verify the format of the .state file and adapt the sequence names for your study
#### this file has been updated for parsing although not confirmed as robust to error
#### sequences = {record.id: str(record.seq) for record in SeqIO.parse(STATE_FILE, "fasta")}
print(f"Loaded {len(sequences)} total sequences (leaves + internal nodes).")

# 2. Load the phylogenetic tree.
# The format=1 ensures branch lengths with scientific notation are read correctly.
t = Tree(TREE_FILE, format=1)
print(f"Loaded tree with {len(t)} leaves and {t.get_tree_root().name} as the root.")

# 3. Traverse the tree and create pairs.
pairs = []
for node in t.traverse("preorder"): # "preorder" starts from the root and goes down
    # The root has no ancestor, so we skip it.
    if node.is_root():
        continue

    ancestor_node = node.up
    ancestor_name = ancestor_node.name
    descendant_name = node.name

    # Ensure both ancestor and descendant have sequences in our dictionary
    if ancestor_name in sequences and descendant_name in sequences:
        pairs.append({
            "ancestor_name": ancestor_name,
            "ancestor_seq": sequences[ancestor_name],
            "descendant_name": descendant_name,
            "descendant_seq": sequences[descendant_name]
        })

# 4. Convert to a DataFrame and save as a TSV file.
df_pairs = pd.DataFrame(pairs)
df_pairs.to_csv(OUTPUT_PAIRS_FILE, sep='\t', index=False)

print("-" * 50)
print(f"Successfully created '{OUTPUT_PAIRS_FILE}' with {len(df_pairs)} pairs.")
print("Columns:", df_pairs.columns.tolist())
print("\nFirst 5 rows of the generated file:")
print(df_pairs.head())
print("-" * 50)

# ==============================================================================
# Part 6: Filter Pairs for Training (with "Identical Pairs" Logic)
# ==============================================================================
import pandas as pd

PAIRS_FILE = 'ancestor_descendant_pairs.tsv'

# Load the full set of generated pairs
df_pairs = pd.read_csv(PAIRS_FILE, sep='\t')

# Identify the mutated and identical pairs
is_mutated = df_pairs['ancestor_seq'] != df_pairs['descendant_seq']
df_mutated = df_pairs[is_mutated]
df_identical = df_pairs[~is_mutated]

print(f"Found {len(df_mutated)} pairs with at least one mutation.")
print(f"Found {len(df_identical)} identical pairs (stasis).")

# --- The Strategy: Balance the dataset ---
# Keep all mutated pairs and a representative subsample of identical pairs.
# This teaches the model that stasis is common, preventing it from over-predicting mutations.
if not df_mutated.empty and len(df_identical) > len(df_mutated):
    # Sample a number of identical pairs equal to the number of mutated pairs
    df_identical_sample = df_identical.sample(n=len(df_mutated), random_state=42)
else:
    # If there are few or no mutated pairs, or fewer identical than mutated,
    # just use all the identical pairs.
    df_identical_sample = df_identical

# Combine them into the final training dataset and shuffle
df_final_training = pd.concat([df_mutated, df_identical_sample]).sample(frac=1).reset_index(drop=True)

# Save the final training set for the next notebook/script if desired
df_final_training.to_csv('final_training_pairs.tsv', sep='\t', index=False)

print(f"\nFinal training set size: {len(df_final_training)}")
num_mutated_final = len(df_final_training[df_final_training['ancestor_seq'] != df_final_training['descendant_seq']])
num_identical_final = len(df_final_training[df_final_training['ancestor_seq'] == df_final_training['descendant_seq']])
print(f"Composition: {num_mutated_final} mutated pairs, {num_identical_final} identical pairs.")

# This `df_final_training` DataFrame is now structured to be the input
# for Transformer model's tokenization and training steps.

# ==============================================================================
# Part 7: Transformer (Seq2Seq Model with encoder/decoder)
# ==============================================================================

import os
from google.colab import drive

# Mount Google Drive
drive.mount('/content/drive')
os.chdir('/content/drive/MyDrive')
print(f"Current working directory: {os.getcwd()}")

# Setup and Data Loading
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import pandas as pd
import numpy as np
import os

# --- Configuration ---
PAIRS_FILE = 'final_training_pairs.tsv'
BATCH_SIZE = 8  # smaller batches use less RAM at once
EPOCHS = 5
MAX_SEQ_LENGTH = 1024 # this is within Colab limits, but then the (gene) sequence pair must be 512/512 in length

# Model Hyperparameters (probably insufficient for any robust model of evolution, but adapted to Colab's limits)
EMBED_DIM = 64
NUM_HEADS = 4
FF_DIM = 256
NUM_ENCODER_LAYERS = 2
NUM_DECODER_LAYERS = 2
DROPOUT_RATE = 0.1

# --- Load and Preprocess Data ---
if not os.path.exists(PAIRS_FILE):
    print(f"ERROR: The data file '{PAIRS_FILE}' was not found.")
    print("Please ensure the previous data preparation steps have been completed successfully.")
else:
    print(f"Loading data from '{PAIRS_FILE}'...")
    df = pd.read_csv(PAIRS_FILE, sep='\t')

    # For demonstration, let's use a smaller subset. Remove this for a full run.
    df = df.sample(n=min(10000, len(df)), random_state=42)
    print(f"Using a subset of {len(df)} pairs for this demonstration.")

    # Create the vocabulary
    # Find all unique characters in both ancestor and descendant sequences
    vocab = set()
    for seq in pd.concat([df['ancestor_seq'], df['descendant_seq']]):
        vocab.update(list(seq))
    vocab = sorted(list(vocab))

    # Add special tokens
    special_tokens = ["[PAD]", "[START]", "[END]"]
    vocab = special_tokens + vocab
    VOCAB_SIZE = len(vocab)

    print(f"Vocabulary size: {VOCAB_SIZE}")
    print(f"Vocabulary: {vocab}")

    # Create token mapping dictionaries
    char_to_token = {char: i for i, char in enumerate(vocab)}
    token_to_char = {i: char for i, char in enumerate(vocab)}

    # --- Tokenization and Padding Function ---
    def tokenize_and_pad(sequences):
        tokenized = []
        for seq in sequences:
            # Add start and end tokens
            current_tokens = [char_to_token["[START]"]]
            current_tokens.extend([char_to_token.get(char, 0) for char in seq])
            current_tokens.append(char_to_token["[END]"])
            tokenized.append(current_tokens)

        # Pad sequences to the max length
        return keras.preprocessing.sequence.pad_sequences(
            tokenized, maxlen=MAX_SEQ_LENGTH, padding="post", value=char_to_token["[PAD]"]
        )

    # Apply the function
    ancestor_vectors = tokenize_and_pad(df["ancestor_seq"].values)
    descendant_vectors = tokenize_and_pad(df["descendant_seq"].values)

    print(f"\nShape of ancestor vectors: {ancestor_vectors.shape}")
    print(f"Shape of descendant vectors: {descendant_vectors.shape}")

    # --- Create tf.data.Dataset ---
    # Prepare inputs and outputs for the model
    # Decoder input is the sequence shifted right (starts with [START])
    # Decoder output is the sequence as is (ends with [END])
    decoder_inputs = descendant_vectors[:, :-1]
    decoder_outputs = descendant_vectors[:, 1:]

    dataset = tf.data.Dataset.from_tensor_slices(
        ((ancestor_vectors, decoder_inputs), decoder_outputs)
    )
    dataset = dataset.batch(BATCH_SIZE).shuffle(buffer_size=1024).prefetch(tf.data.AUTOTUNE)

# Transformer Building Blocks (Custom Layers)

# --- Positional Embedding Layer ---
# Combines token embedding with a fixed positional encoding.
class PositionalEmbedding(layers.Layer):
    def __init__(self, vocab_size, embed_dim, maxlen, **kwargs):
        super().__init__(**kwargs)
        self.maxlen = maxlen
        self.embed_dim = embed_dim
        self.vocab_size = vocab_size
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
        # Add the positional encoding, slicing to match the input length
        return x + self.pos_emb[:, :maxlen, :]

# --- Transformer Encoder Layer ---
class TransformerEncoder(layers.Layer):
    def __init__(self, embed_dim, ff_dim, num_heads, **kwargs):
        super().__init__(**kwargs)
        self.attn = layers.MultiHeadAttention(num_heads=num_heads, key_dim=embed_dim)
        self.ffn = keras.Sequential(
            [layers.Dense(ff_dim, activation="relu"), layers.Dense(embed_dim),]
        )
        self.layernorm1 = layers.LayerNormalization(epsilon=1e-6)
        self.layernorm2 = layers.LayerNormalization(epsilon=1e-6)
        self.dropout1 = layers.Dropout(DROPOUT_RATE)
        self.dropout2 = layers.Dropout(DROPOUT_RATE)

    def call(self, inputs, training=False):
        # Multi-head self-attention
        attn_output = self.attn(inputs, inputs)
        attn_output = self.dropout1(attn_output, training=training)
        out1 = self.layernorm1(inputs + attn_output)

        # Feed-forward network
        ffn_output = self.ffn(out1)
        ffn_output = self.dropout2(ffn_output, training=training)
        return self.layernorm2(out1 + ffn_output)

# --- Transformer Decoder Layer ---
class TransformerDecoder(layers.Layer):
    def __init__(self, embed_dim, ff_dim, num_heads, **kwargs):
        super().__init__(**kwargs)
        self.self_attn = layers.MultiHeadAttention(num_heads=num_heads, key_dim=embed_dim)
        self.cross_attn = layers.MultiHeadAttention(num_heads=num_heads, key_dim=embed_dim)
        self.ffn = keras.Sequential(
            [layers.Dense(ff_dim, activation="relu"), layers.Dense(embed_dim),]
        )
        self.layernorm1 = layers.LayerNormalization(epsilon=1e-6)
        self.layernorm2 = layers.LayerNormalization(epsilon=1e-6)
        self.layernorm3 = layers.LayerNormalization(epsilon=1e-6)
        self.dropout1 = layers.Dropout(DROPOUT_RATE)
        self.dropout2 = layers.Dropout(DROPOUT_RATE)
        self.dropout3 = layers.Dropout(DROPOUT_RATE)

    def call(self, inputs, encoder_outputs, training=False):
        # Causal (masked) self-attention
        causal_mask = self.get_causal_attention_mask(inputs)
        self_attn_output = self.self_attn(
            query=inputs, value=inputs, key=inputs, attention_mask=causal_mask
        )
        self_attn_output = self.dropout1(self_attn_output, training=training)
        out1 = self.layernorm1(inputs + self_attn_output)

        # Cross-attention (attends to encoder output)
        cross_attn_output = self.cross_attn(
            query=out1, value=encoder_outputs, key=encoder_outputs
        )
        cross_attn_output = self.dropout2(cross_attn_output, training=training)
        out2 = self.layernorm2(out1 + cross_attn_output)

        # Feed-forward network
        ffn_output = self.ffn(out2)
        ffn_output = self.dropout3(ffn_output, training=training)
        return self.layernorm3(out2 + ffn_output)

    def get_causal_attention_mask(self, inputs):
        input_shape = tf.shape(inputs)
        batch_size, sequence_length = input_shape[0], input_shape[1]
        i = tf.range(sequence_length)[:, tf.newaxis]
        j = tf.range(sequence_length)
        mask = tf.cast(i >= j, dtype="int32")
        return mask

# Build and Train the Transformer Model

# --- Build the full model ---
encoder_inputs = keras.Input(shape=(None,), dtype="int32", name="ancestor")
decoder_inputs = keras.Input(shape=(None,), dtype="int32", name="descendant")

# Encoder path
encoder_embedding = PositionalEmbedding(VOCAB_SIZE, EMBED_DIM, MAX_SEQ_LENGTH)(encoder_inputs)
x = encoder_embedding
for _ in range(NUM_ENCODER_LAYERS):
    x = TransformerEncoder(EMBED_DIM, FF_DIM, NUM_HEADS)(x)
encoder_outputs = x

# Decoder path
decoder_embedding = PositionalEmbedding(VOCAB_SIZE, EMBED_DIM, MAX_SEQ_LENGTH)(decoder_inputs)
x = decoder_embedding
for _ in range(NUM_DECODER_LAYERS):
    x = TransformerDecoder(EMBED_DIM, FF_DIM, NUM_HEADS)(x, encoder_outputs)

# Final output layer
# This projects the decoder's output back to the vocabulary space to retrieve probabilities for each token
output_logits = layers.Dense(VOCAB_SIZE, name="logits")(x)

# Create the model
transformer = keras.Model([encoder_inputs, decoder_inputs], output_logits, name="viral_transformer")
transformer.summary()

# --- Compile and Train ---
transformer.compile(
    optimizer=keras.optimizers.Adam(learning_rate=1e-4),
    loss=keras.losses.SparseCategoricalCrossentropy(from_logits=True),
    metrics=[keras.metrics.SparseCategoricalAccuracy()],
)

# Optional: Add a callback to save the model
model_checkpoint_callback = keras.callbacks.ModelCheckpoint(
    filepath="viral_transformer_checkpoint.weights.h5",
    save_weights_only=True,
    monitor='loss',
    mode='max',
    save_best_only=True)

# Note: Training slow on a CPU. GPU is highly recommended.
history = transformer.fit(
    dataset,
    epochs=EPOCHS,
    callbacks=[model_checkpoint_callback]
    # To use a validation split, create a separate validation dataset
    # validation_data=val_dataset
)

# Inference and Prediction
# Load best weights if you trained with the callback
# transformer.load_weights("viral_transformer_checkpoint.weights.h5")

def decode_sequence(input_sentence):
    # Tokenize the input ancestor sequence
    tokenized_input = tokenize_and_pad([input_sentence])[0]

    # Initialize the decoder's input with the [START] token
    decoded_sentence = [char_to_token["[START]"]]

    for i in range(MAX_SEQ_LENGTH):
        # Prepare inputs for the model
        ancestor_input = tf.constant([tokenized_input], dtype=tf.int32)
        decoder_input = tf.constant([decoded_sentence], dtype=tf.int32)

        # Get the model's prediction (logits for the next token)
        predictions = transformer([ancestor_input, decoder_input], training=False)

        # Sample the token with the highest probability
        next_token_id = tf.argmax(predictions[0, i, :]).numpy()

        # Append the predicted token to the sequence
        decoded_sentence.append(next_token_id)

        # If the model predicts the [END] token, stop
        if next_token_id == char_to_token["[END]"]:
            break

    # Convert the token sequence back to characters
    return "".join(token_to_char.get(token, "?") for token in decoded_sentence[1:-1])

def compare_sequences(ancestor, actual, predicted):
    print("-" * 50)
    print(f"Ancestor:  {ancestor[:80]}...")
    print(f"Actual:    {actual[:80]}...")
    print(f"Predicted: {predicted[:80]}...")
    print("-" * 50)
    mutations = []
    for i, (a, p) in enumerate(zip(ancestor, predicted)):
        if a != p:
            mutations.append(f"  - Position {i+1}: {a} -> {p}")

    if mutations:
        print("Predicted Mutations:")
        print("\n".join(mutations))
    else:
        print("Predicted sequence is identical to the ancestor.")
    print("-" * 50)


# --- Test with a couple of examples from our dataset ---
for i in range(2):
    idx = np.random.randint(0, len(df))
    ancestor_seq = df["ancestor_seq"].iloc[idx]
    actual_descendant_seq = df["descendant_seq"].iloc[idx]

    predicted_descendant_seq = decode_sequence(ancestor_seq)

    compare_sequences(ancestor_seq, actual_descendant_seq, predicted_descendant_seq)
