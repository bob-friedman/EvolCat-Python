# --- New Cell ---

# 1. Environment Setup
# After this cell runs, the kernel will restart automatically.
# You must then run the subsequent cells manually.
!pip install -q condacolab
import condacolab
condacolab.install()

# --- New Cell ---

# 2. Install Dependencies
# This cell runs after the kernel has restarted.
import os
from google.colab import drive

# Mount Google Drive to save our work
drive.mount('/content/drive')
os.chdir('/content/drive/MyDrive')
print(f"Current working directory: {os.getcwd()}")

# Install required bioinformatics tools
!conda install -y -c bioconda -c conda-forge mafft iqtree mmseqs2
!pip install -q biopython ete3

# --- New Cell ---

# 3. Targeted Sequence Retrieval
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

# --- Configuration ---
Entrez.email = "your.email@example.com" # Required by NCBI
samples_per_group = 5
search_groups = {
    "Ancestral_Wuhan": '("2019/12/01"[PDAT] : "2020/04/30"[PDAT])',
    "Alpha_Variant": '"Alpha variant"[All Fields]',
    "Delta_Variant": '"Delta variant"[All Fields]',
    "Early_Omicron_BA.2": '"Omicron"[All Fields] AND ("2022/03/01"[PDAT] : "2022/06/30"[PDAT])',
    "Late_Omicron_JN.1": '"Omicron JN.1"[All Fields]'
}
base_query = '("SARS-CoV-2"[Organism]) AND "complete genome"[Title] AND 29000:30000[SLEN]'
all_accessions = []

# --- Retrieve Accession Numbers ---
print("Fetching accession numbers...")
for group, term in search_groups.items():
    full_query = f"({base_query}) AND ({term})"
    handle = Entrez.esearch(db="nucleotide", term=full_query, retmax=samples_per_group, sort="relevance")
    record = Entrez.read(handle)
    handle.close()
    if id_list := record["IdList"]:
        summary_handle = Entrez.esummary(db="nucleotide", id=",".join(id_list))
        summary_record = Entrez.read(summary_handle)
        summary_handle.close()
        all_accessions.extend([item["AccessionVersion"] for item in summary_record])

# --- Retrieve GenBank Records & Extract Spike Proteins ---
print("Fetching GenBank records and extracting Spike proteins...")
handle = Entrez.efetch(db="nucleotide", id=all_accessions, rettype="gb", retmode="text")
spike_records = []
for record in SeqIO.parse(handle, "genbank"):
    for feature in record.features:
        if feature.type == "CDS" and 'S' in feature.qualifiers.get("gene", []):
            spike_seq = feature.qualifiers['translation'][0]
            spike_records.append(SeqRecord(spike_seq, id=record.id, description=""))
            break
handle.close()

# Save the extracted sequences to a file
SeqIO.write(spike_records, "spike_proteins.fasta", "fasta")
print(f"--> Extracted {len(spike_records)} Spike sequences to spike_proteins.fasta")

# --- New Cell ---

import os
from google.colab import drive

# Mount Google Drive to save our work
drive.mount('/content/drive')
os.chdir('/content/drive/MyDrive')
print(f"Current working directory: {os.getcwd()}")

# Install necessary bioinformatics tools
!conda install -y -c bioconda -c conda-forge iqtree
# !pip install -q biopython ete3

# 4. MSA, Phylogeny, and Ancestral State Reconstruction (ASR)
# --- Step 1: Cluster sequences to get unique representatives ---
# This removes any 100% identical sequences before alignment.
print("\nClustering sequences to find unique representatives...")
!mmseqs easy-cluster spike_proteins.fasta cluster_results /content/tmp_mmseqs --min-seq-id 1.0
!mv cluster_results_rep_seq.fasta sequence_representatives.fas

# --- Step 2: Create a Multiple Sequence Alignment (MSA) ---
print(f"\nAligning representative sequences with MAFFT...")
!mafft --auto sequence_representatives.fas > sequence_aligned.fas

# --- Step 3: Run IQ-TREE for phylogeny and bootstrap values ---
# The '-B 1000' flag performs 1000 bootstrap replicates to assess branch support.
print("Running IQ-TREE with 1000 bootstrap replicates...")
!iqtree -s sequence_aligned.fas -m AUTO -B 1000 -T AUTO -asr --redo -safe

print("\nBootstrap analysis complete.")
print("The results with support values are typically in 'sequence_aligned.fas.contree'.")

print("\nIQ-TREE run complete.")
!ls -lh sequence_aligned.fas*

# 5. Create Ancestor-Descendant Pairs (Method using Pandas)
from ete3 import Tree
from Bio import SeqIO
import pandas as pd

# --- File Definitions ---
# The .treefile contains the correct tree topology AND the 'NodeX' labels.
TREE_FILE = "sequence_aligned.fas.treefile"
# The .state file contains the ancestral sequence data in a table format.
STATE_FILE = "sequence_aligned.fas.state"
# The original FASTA file contains the tip (leaf) sequences.
TIPS_FILE = "sequence_aligned.fas"
OUTPUT_PAIRS_FILE = "ancestor_descendant_pairs.tsv"

print(f"\nParsing tree from '{TREE_FILE}' and tabular data from '{STATE_FILE}'...")

# 1. Parse the .state file using pandas, which is designed for tables.
# This is for reading the ancestral sequence data.
try:
    df_state = pd.read_csv(STATE_FILE, sep='\t', comment='#')
except Exception as e:
    raise IOError(f"Pandas failed to read the state file. Error: {e}")

# Reconstruct the full ancestral sequences from the table.
# We group by the 'Node' column and join the 'State' characters.
ancestral_seqs = df_state.groupby('Node')['State'].apply(''.join).to_dict()

# 2. Retrieve the sequences for the tips (leaves) of the tree.
tip_seqs = {record.id: str(record.seq) for record in SeqIO.parse(TIPS_FILE, "fasta")}

# 3. Combine both dictionaries into a single master dictionary of all sequences.
all_sequences = {**ancestral_seqs, **tip_seqs}

if not ancestral_seqs:
    raise RuntimeError("Successfully parsed .state file, but no ancestral sequences were reconstructed.")

# 4. Load the tree directly from the .treefile. The node names will now match.
t = Tree(TREE_FILE, format=1)

# 5. Create pairs by matching names, which are now guaranteed to be consistent.
pairs = []
for node in t.traverse("preorder"):
    if not node.is_root():
        ancestor_name = node.up.name
        descendant_name = node.name

        # This check will now succeed because all names and sequences are present.
        if ancestor_name in all_sequences and descendant_name in all_sequences:
            pairs.append({
                "ancestor_name": ancestor_name,
                "ancestor_seq": all_sequences[ancestor_name],
                "descendant_name": descendant_name,
                "descendant_seq": all_sequences[descendant_name]
            })

# 6. Save the results.
df_pairs = pd.DataFrame(pairs)
df_pairs.to_csv(OUTPUT_PAIRS_FILE, sep='\t', index=False)

if len(df_pairs) > 0:
    print(f"--> Successfully created '{OUTPUT_PAIRS_FILE}' with {len(df_pairs)} pairs.")
else:
    print(f"--> ERROR: Failed to create pairs. A critical mismatch still exists.")

# --- New Cell ---

# 6. Filter Pairs for Transformer Training
import pandas as pd

df_pairs = pd.read_csv('ancestor_descendant_pairs.tsv', sep='\t')
df_mutated = df_pairs[df_pairs['ancestor_seq'] != df_pairs['descendant_seq']]
df_identical = df_pairs[df_pairs['ancestor_seq'] == df_pairs['descendant_seq']]

print(f"Found {len(df_mutated)} mutated pairs and {len(df_identical)} identical pairs.")
if not df_mutated.empty and len(df_identical) > len(df_mutated):
    df_identical_sample = df_identical.sample(n=len(df_mutated), random_state=42)
else:
    df_identical_sample = df_identical

df_final_training = pd.concat([df_mutated, df_identical_sample]).sample(frac=1).reset_index(drop=True)
df_final_training.to_csv('final_training_pairs.tsv', sep='\t', index=False)
print(f"--> Final balanced training set created with {len(df_final_training)} pairs.")

### Interpretability & Visualization
# This is a modified code block specifically to retrieve Attention Scores
# Depends on the final pairs of sequence data only

# --- New Cell ---

# 1. Setup and Imports
from google.colab import drive
import os

# Mount Google Drive and set the working directory
drive.mount('/content/drive')
os.chdir('/content/drive/MyDrive')
print(f"Current working directory: {os.getcwd()}")

!pip install -q biopython

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.PDB import PDBParser, PDBIO
import urllib.request

print("\n--- Starting Interpretability Analysis ---")

# --- 2. Configuration & Rebuild Vocabulary ---
PAIRS_FILE = 'final_training_pairs.tsv'
MODEL_CHECKPOINT_FILE = 'viral_transformer.weights.h5'
MAX_SEQ_LENGTH = 1300
EMBED_DIM = 64
NUM_HEADS = 4
FF_DIM = 256
DROPOUT_RATE = 0.1
NUM_ENCODER_LAYERS = 2
NUM_DECODER_LAYERS = 2

# Check if necessary files exist
if not os.path.exists(PAIRS_FILE) or not os.path.exists(MODEL_CHECKPOINT_FILE):
    raise FileNotFoundError("Ensure 'final_training_pairs.tsv' and 'viral_transformer.weights.h5' are in the current directory.")

df_for_vocab = pd.read_csv(PAIRS_FILE, sep='\t')
vocab = set()
for seq in pd.concat([df_for_vocab['ancestor_seq'].dropna(), df_for_vocab['descendant_seq'].dropna()]):
    vocab.update(list(str(seq)))
vocab = sorted(list(vocab))
special_tokens = ["[PAD]", "[START]", "[END]"]
vocab = special_tokens + vocab
VOCAB_SIZE = len(vocab)
char_to_token = {char: i for i, char in enumerate(vocab)}

# --- 3. Rebuild the FULL Model and Load Weights ---
# We define the model components exactly as during training.
class PositionalEmbedding(layers.Layer):
    def __init__(self, vocab_size, embed_dim, maxlen, **kwargs):
        super().__init__(**kwargs)
        self.token_emb = layers.Embedding(input_dim=vocab_size, output_dim=embed_dim)
        pos_encoding = np.zeros((maxlen, embed_dim))
        positions = np.arange(maxlen)[:, np.newaxis]
        div_term = np.exp(np.arange(0, embed_dim, 2) * -(np.log(10000.0) / embed_dim))
        pos_encoding[:, 0::2] = np.sin(positions * div_term)
        pos_encoding[:, 1::2] = np.cos(positions * div_term)
        self.pos_emb = tf.cast(pos_encoding[np.newaxis, ...], dtype=tf.float32)
    def call(self, x):
        return self.token_emb(x) + self.pos_emb[:, :tf.shape(x)[-1], :]

class TransformerEncoder(layers.Layer):
    def __init__(self, embed_dim, ff_dim, num_heads, **kwargs):
        super().__init__(**kwargs)
        self.attn = layers.MultiHeadAttention(num_heads=num_heads, key_dim=embed_dim)
        self.ffn = keras.Sequential([layers.Dense(ff_dim, activation="relu"), layers.Dense(embed_dim)])
        self.layernorm1 = layers.LayerNormalization(epsilon=1e-6)
        self.layernorm2 = layers.LayerNormalization(epsilon=1e-6)
        self.dropout1 = layers.Dropout(DROPOUT_RATE)
        self.dropout2 = layers.Dropout(DROPOUT_RATE)
    def call(self, inputs, training=False, return_attention=False):
        # Always run with return_attention_scores=True to get both outputs
        attn_output, attn_scores = self.attn(inputs, inputs, return_attention_scores=True)
        attn_output = self.dropout1(attn_output, training=training)
        out1 = self.layernorm1(inputs + attn_output)
        ffn_output = self.ffn(out1)
        ffn_output = self.dropout2(ffn_output, training=training)
        final_output = self.layernorm2(out1 + ffn_output)

        # Return scores only if requested
        if return_attention:
            return final_output, attn_scores
        return final_output

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
    def get_causal_attention_mask(self, inputs):
        input_shape = tf.shape(inputs)
        i = tf.range(input_shape[1])[:, tf.newaxis]
        j = tf.range(input_shape[1])
        mask = tf.cast(i >= j, dtype="int32")
        return mask
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

# Build the FULL original model structure
keras.backend.clear_session()
encoder_inputs = keras.Input(shape=(None,), dtype="int32", name="ancestor")
decoder_inputs = keras.Input(shape=(None,), dtype="int32", name="descendant")
encoder_embedding = PositionalEmbedding(VOCAB_SIZE, EMBED_DIM, MAX_SEQ_LENGTH, name="encoder_embedding")(encoder_inputs)
x = encoder_embedding
for i in range(NUM_ENCODER_LAYERS):
    x = TransformerEncoder(EMBED_DIM, FF_DIM, NUM_HEADS, name=f"encoder_{i}")(x)
encoder_outputs = x
decoder_embedding = PositionalEmbedding(VOCAB_SIZE, EMBED_DIM, MAX_SEQ_LENGTH, name="decoder_embedding")(decoder_inputs)
x = decoder_embedding
for i in range(NUM_DECODER_LAYERS):
    x = TransformerDecoder(EMBED_DIM, FF_DIM, NUM_HEADS, name=f"decoder_{i}")(x, encoder_outputs)
output_logits = layers.Dense(VOCAB_SIZE, name="logits")(x)
full_model = keras.Model([encoder_inputs, decoder_inputs], output_logits, name="viral_transformer")

# Load the weights into the full model
full_model.load_weights(MODEL_CHECKPOINT_FILE)
print(f"Successfully loaded weights into the full model from '{MODEL_CHECKPOINT_FILE}'.")

# --- Utility Function for Tokenization ---
def tokenize_and_pad(sequences, maxlen, char_map):
    tokenized = []
    for seq in sequences:
        current_tokens = [char_map["[START]"]]
        current_tokens.extend([char_map.get(char, 0) for char in str(seq)])
        current_tokens.append(char_map["[END]"])
        tokenized.append(current_tokens)
    return keras.preprocessing.sequence.pad_sequences(tokenized, maxlen=maxlen, padding="post", dtype='int32')

# --- 4. Extract and Save Attention Scores from the Main Model ---
print("Extracting attention scores from the final encoder layer...")

# Get the layers we need directly from the successfully loaded model
embedding_layer = full_model.get_layer("encoder_embedding")
encoder_0 = full_model.get_layer("encoder_0")
encoder_1 = full_model.get_layer("encoder_1") # This is the last encoder

# Retrieve a sample sequence (existing code does this)
df = pd.read_csv(PAIRS_FILE, sep='\t').dropna()
mutated_df = df[df['ancestor_seq'] != df['descendant_seq']]
sample_ancestor = mutated_df["ancestor_seq"].iloc[0] if not mutated_df.empty else df["ancestor_seq"].iloc[0]
tokenized_ancestor = tokenize_and_pad([sample_ancestor], MAX_SEQ_LENGTH, char_to_token)

# Manually pass the data through the encoder layers
embedded_input = embedding_layer(tokenized_ancestor)
x = encoder_0(embedded_input)
# Call the final encoder with our new 'return_attention=True' flag
_, last_layer_attention = encoder_1(x, return_attention=True)

# Process the scores to get one value per residue
# last_layer_attention shape: (batch, heads, seq_len, seq_len)
# We average across the heads
avg_attention_matrix = np.mean(last_layer_attention[0], axis=0)

# We want the attention each token *receives*, so average over the columns (axis=0)
per_residue_attention_full = np.mean(avg_attention_matrix, axis=0)

# --- 5. Save the final attention scores for the FULL sequence length to CSV ---
attention_df = pd.DataFrame({
    'Position': np.arange(1, len(per_residue_attention_full) + 1),
    'Attention': per_residue_attention_full
})
attention_df.to_csv('attention_scores.csv', index=False)

print("\nSuccessfully created 'attention_scores.csv'.")
print(attention_df.head())

# (Optional) Method to Find Vocabulary Size of a Model
# Prints shapes which defines a model

from google.colab import drive
import os

# Mount Google Drive and set the working directory
drive.mount('/content/drive')
os.chdir('/content/drive/MyDrive')
print(f"Current working directory: {os.getcwd()}")

import h5py

def print_h5_structure(h5_file_path):
    """
    Prints the hierarchical structure of an H5 file, including dataset shapes.
    """
    try:
        with h5py.File(h5_file_path, 'r') as f:
            print(f"--- Structure of {h5_file_path} ---")
            def print_attrs(name, obj):
                # This function is called for every object in the file
                indent = '  ' * name.count('/')
                if isinstance(obj, h5py.Dataset):
                    print(f"{indent}Dataset: {name}, Shape: {obj.shape}")
                else: # It's a Group
                    print(f"{indent}Group: {name}")

            f.visititems(print_attrs)
        print("--- End of Structure ---")

    except FileNotFoundError:
        print(f"ERROR: Model file not found at '{h5_file_path}'")
    except Exception as e:
        print(f"An error occurred: {e}")

# --- Inspect the model file ---
h5_path = 'viral_transformer.weights.h5'
print_h5_structure(h5_path)