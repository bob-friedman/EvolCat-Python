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
*   **Bgzip:** Required method of compression for downstream tasks. (`pip install bgzip`)
*   **BCFTools:** Utilities for variant calling and manipulating VCF files. (`conda install bcftools`)

### Conceptual Knowledge:
*   Basic understanding of genomic data (FASTA, VCF formats).
*   Familiarity with phylogenetic concepts (ancestral sequence reconstruction, phylogenetic trees, mutations).
*   Basic knowledge of machine learning concepts (Transformers, sequence-to-sequence models, embeddings, attention mechanisms).
*   Understanding of basic Linux/command-line operations.

## Expected Outcomes

This tutorial goal is for:
*   Understanding the workflow for collecting and processing SARS-CoV-2 genomic data, specifically focusing on creating ancestor-descendant sequence pairs.
*   Learning how to perform feature engineering using Pandas to prepare data for a machine learning model.
*   Gaining insights into the conceptual framework of using Transformer models for predicting ancestral genomic sequences.
*   Have a trained Transformer model in TensorFlow capable of predicting ancestral sequences from descendant sequences.
*   Be able to adapt this pipeline for your own research questions in viral genomics or similar sequence-based tasks.
*   Understanding how to prepare and preprocess data for a nucleotide-based Transformer model.

## Scaling to Full Dataset of SARS-CoV-2 Genomes

Loop & Sample: A Monte Carlo Approach with parallelism and efficiency:
1. Execute a loop with a large number of iterations.
2. Randomly Select a Node: In each iteration, pick a random internal node from the clade's tree.
3. Filter by Size: Check the number of leaves descending from this randomly selected node, if it is within a predefined manageable range (i.e., 500 - 5,000 leaves) then proceed, otherwise, discard and select another node.
4. Extract sub-subtree for this manageable internal node.
5. Save inferred ancestral sequences for the nodes in this analyzed subtree. Verify node identifiers are unique and can be mapped back to the master clade tree.
6. Assess if this iterative analysis is still yielding new or significantly different ancestral state information for the nodes within the target clade. Stop when results appear to converge or after a sufficient number of iterations.

## Data Collection

This section outlines the steps to collect and preprocess SARS-CoV-2 data to create a dataset of (ancestor, descendant) sequence pairs.

### 1. Database Construction with UShER

The initial step involves constructing or updating a comprehensive SARS-CoV-2 database using UShER (>8 million samples). This database typically includes a global phylogenetic tree (`public-latest.all.masked.pb.gz`) and corresponding mutation-annotated VCF files.

```bash
# Download the latest UShER tree and VCF (example paths)
# Retrieve files at: https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/
wget https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.all.masked.pb.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.all.masked.vcf.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.metadata.tsv.gz

# It's good practice to decompress the VCF for easier processing by some tools,
gunzip public-latest.all.masked.vcf.gz
gunzip public-latest.all.masked.pb.gz
gunzip public-latest.metadata.tsv.gz
```
*Note: Ensure UShER is properly installed and configured. The paths to the database files should be specified.*

### 2. Downloading Reference Genome

We need a reference genome for alignment and variant calling. The Wuhan-Hu-1 genome (NC_045512.2) is commonly used for SARS-CoV-2.

```bash
# Download reference genome and annotation using NCBI Datasets CLI
ncbi-datasets-cli download genome accession NC_045512.2 --include genome --filename SC2_ref.zip
unzip SC2_ref.zip
mv ncbi_dataset/data/NC_045512.2/NC_045512.2.fna ./SC2_ref.fasta
rm -rf ncbi_dataset SC2_ref.zip README.md # Clean up
```

### 3. Extracting Ancestor-Descendant Sequence Pairs using UShER and `matUtils`

## Retrieving Subset of the Full Dataset

Examples:
*   Retrieving all clade names in virus data:
```bash
matUtils summary --input-mat public-latest.all.masked.pb --clades clades.tsv #
```

*   Retrieving all samples in virus data:
```bash
matUtils summary --input-mat public-latest.all.masked.pb --samples samples.txt
```

*   Retrieving Clade 23E (XBB.2.3) in virus data (~16,362 samples of >8M in full data):
```bash
matUtils extract -i public-latest.all.masked.pb -c 'XBB.2.3' -v xbb23.vcf
```

*   Verify sample count in metadata:
```bash
cat public-latest.metadata.tsv | grep -nr "XBB.2.3" | wc -l
```

*   Retrieve subset of data for this Clade:
```bash
cat public-latest.metadata.tsv | grep -nr "XBB.2.3" > public-xbb23-latest.metadata.tsv
echo "Creating accession list..."
cat public-xbb23-latest.metadata.tsv | cut -f2 | grep -v "genbank_accession" | uniq > accession_list.txt
echo "accession_list.txt created with $(wc -l < accession_list.txt) entries."
```

*   Tool for compressing VCF file for downstream tasks:
```bash
bgzip --force xbb23.vcf
```

BCFTools is available in Conda:
*   And an index file that accompanies the compressed VCF file:
```bash
bcftools index --tbi xbb23.vcf.gz
```

Retrieval method for data in VCF file (Conceptual):
```python
# Retrieve data in vcf file
conda install -q -y pysam 
import pysam
from pysam import VariantFile

vcf_filepath = 'xbb23.vcf.gz'
vcf = VariantFile(vcf_filepath)

print(f"Successfully opened {vcf_filepath}")
print("-" * 30)
print(f"Reference sequence in VCF: {vcf.header.contigs}")
print(f"Samples in VCF: {list(vcf.header.samples)}")
print("-" * 30)

variant_count = 0
log_interval = 500
for record in vcf.fetch():
    variant_count += 1
    # Examine samples for this variant
    # The 'GT' field in the FORMAT column tells us the genotype (0=ref, 1=alt)
    samples_with_variant = []
    for sample in record.samples.values():
        # sample['GT'] returns a tuple, e.g., (1, 1) for a homozygous alt
        # We check if the first alternate allele (1) is present in the genotype
        if 1 in sample['GT']:
            samples_with_variant.append(sample.name)

    # Only log detailed info periodically
    if variant_count % log_interval == 0:
        chrom = record.chrom
        pos = record.pos
        ref_allele = record.ref
        alt_alleles = record.alts

        print(f"Processed {variant_count} variants...")
        print(f"\nVariant at {chrom}:{pos} | REF: {ref_allele} | ALT: {alt_alleles[0] if alt_alleles else ''}")

        if samples_with_variant:
            print(f"  > Found in {len(samples_with_variant)} samples: {samples_with_variant[:3]}...") # Print first 3
        else:
            print("  > Not found in any sample in this VCF (might be an ancestral variant).")

print(f"\nProcessing complete. Analyzed {variant_count} total variants.")
vcf.close()
```

## Retrieval of NCBI Genome Data (Conceptual)

*   Create accession_list.txt of accession numbers in a list, for example:
GCF_000864765.1
GCF_000840245.1
Download genome data (zip file) in the accession list
```bash
ncbi-datasets-cli download genome accession --inputfile accession_list.txt --filename downloaded_genomes.zip
```

Or download by taxon:
```bash
ncbi-datasets-cli download virus genome taxon 'Severe acute respiratory syndrome coronavirus 2' --reference --include genome --filename sars-cov-2-genomes.zip
unzip -o sars-cov-2-genomes.zip # -o for overwrite without prompt
```

Find all FASTA files (.fna) in the download directory, concatenates into a single file for downstream tools:
```bash
echo "Consolidating downloaded sequences..."
find ncbi_dataset/data -name "*.fna" -exec cat {} + > retrieved_sequences.fasta
echo "All sequences consolidated into retrieved_sequences.fasta"
```

## Cluster Sequence Data by Similarity (Conceptual)

Construct multiple sequence alignment of sequences:
*   Cluster sequences with 100% identity (fast search of large datasets)
*   Alternatively, for very large dataset cluster at lower than 1.0, such as ~0.99 (min-seq-id)
*   Uses Fast Fourier Transform approximation for rapid homologous segment identification
```bash
mmseqs easy-cluster --min-seq-id 1.0
```

## Multiple Sequence Alignment (MSA) of above Data Output (Conceptual)
MAFFT selects a fast method of sequence alignment:
```bash
mafft --auto retrieved_sequences.fasta > aligned_sequences.fasta
```

## Construct Phylogenetic Tree of MSA (Conceptual)
Construct phylogenetic tree by IQ-TREE:
*   Scales tree by divergence (substitutions per site; GTR+F+G4 model is common for viruses
*   Use distinct prefix for output files: 'iqtree_run'; tree file is iqtree_run.treefile
*   Leaf names in tree must exactly match sequence names in aligned_sequences.fasta
```bash
iqtree -s aligned_sequences.fasta -m GTR+F+G4 -nt AUTO --threads-max AUTO -mem 16G -pre iqtree_run
```

## Ancestral Sequence Reconstruction by TreeTime (Conceptual):
*   Dependent on Newick tree format and Fasta multiple sequence alignment (set for automatic resolution of polytomies)
*   Output includes: annotated_tree.nexus with inferred ancestral sequences and mutations annotated on branches (FigTree) and ancestral_sequences.fasta of inferred sequences at internal nodes
```bash
treetime ancestral --aln aligned_sequences.fasta --tree iqtree_run.treefile --outdir treetime_asr_output
```

## Construct Data of Ancestral/Descendant Pairs of Sequences (Conceptual):
First, install the necessary library
```bash
!pip install -q dendropy
```

```python
import dendropy

# --- Configuration ---
# These are the output files from the previous TreeTime step
TREE_FILE = 'treetime_asr_output/annotated_tree.nexus'
ANCESTRAL_SEQS_FILE = 'treetime_asr_output/ancestral_sequences.fasta'
OUTPUT_PAIRS_FILE = 'ancestor_descendant_pairs.tsv'

print(f"Loading tree from {TREE_FILE}...")
# Load the tree, specifying the schema
tree = dendropy.Tree.get(path=TREE_FILE, schema="nexus")

print(f"Loading sequences from {ANCESTRAL_SEQS_FILE}...")
# Load the ancestral sequences into a dictionary for easy lookup
# DendroPy's FastaReader is a good choice
seqs = dendropy.DnaCharacterMatrix.get(path=ANCESTRAL_SEQS_FILE, schema="fasta")
seq_dict = {taxon.label: str(seqs[taxon]) for taxon in seqs.taxon_namespace}

print("Generating ancestor-descendant pairs...")
pairs_count = 0
with open(OUTPUT_PAIRS_FILE, 'w') as f_out:
    # Write a header
    f_out.write("ancestor_id\tdescendant_id\tancestor_seq\tdescendant_seq\n")

    # Iterate over every node in the tree
    for node in tree.postorder_node_iter():
        # A node is an ancestor if it is not a leaf (i.e., it has children)
        if not node.is_leaf():
            ancestor_id = node.label
            
            # Skip if the ancestor node's sequence is not available for some reason
            if ancestor_id not in seq_dict:
                continue
            
            ancestor_seq = seq_dict[ancestor_id]

            # Now, iterate through its direct children (descendants)
            for child_node in node.child_nodes():
                descendant_id = child_node.label
                
                # The descendant could be another internal node or a leaf (tip)
                if descendant_id not in seq_dict:
                    continue

                descendant_seq = seq_dict[descendant_id]
                
                # Write the pair to our output file
                f_out.write(f"{ancestor_id}\t{descendant_id}\t{ancestor_seq}\t{descendant_seq}\n")
                pairs_count += 1

print(f"Processing complete. Wrote {pairs_count} pairs to {OUTPUT_PAIRS_FILE}.")
```

## Feature Engineering with Pandas

Now, we'll process the files generated by `matUtils` into a structured Pandas DataFrame. This DataFrame will be the input for our Transformer model.

### Script to Create DataFrame for Ancestor-Descendant Pairs (Conceptual)

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

## Implement the Transformer based method #1 (Conceptual)

This section provides a Keras/TensorFlow implementation of the Transformer model.

```python
# Setup and Data Loading
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import pandas as pd
import numpy as np
import os

# --- Configuration ---
PAIRS_FILE = 'ancestor_descendant_pairs.tsv'
BATCH_SIZE = 64
EPOCHS = 10 # Start with a few epochs, increase for better results
MAX_SEQ_LENGTH = 30000 # SARS-CoV-2 is ~29.9k bp. Pad to this length.

# Model Hyperparameters
EMBED_DIM = 128      # Embedding dimension for each token
NUM_HEADS = 8        # Number of attention heads
FF_DIM = 512         # Hidden layer size in feed forward network inside transformer
NUM_ENCODER_LAYERS = 3
NUM_DECODER_LAYERS = 3
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
    # Prepare inputs and outputs for the "teacher forcing" model
    # Decoder input is the sequence shifted right (starts with [START])
    # Decoder output is the sequence as is (ends with [END])
    decoder_inputs = descendant_vectors[:, :-1]
    decoder_outputs = descendant_vectors[:, 1:]

    dataset = tf.data.Dataset.from_tensor_slices(
        ((ancestor_vectors, decoder_inputs), decoder_outputs)
    )
    dataset = dataset.batch(BATCH_SIZE).shuffle(buffer_size=1024).prefetch(tf.data.AUTOTUNE)

```

```python
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
```

```python
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
# This projects the decoder's output back to the vocabulary space to get probabilities for each token
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
    filepath="viral_transformer_checkpoint.h5",
    save_weights_only=True,
    monitor='val_sparse_categorical_accuracy',
    mode='max',
    save_best_only=True)

# Note: Training will be slow on a CPU. A GPU is highly recommended.
history = transformer.fit(
    dataset,
    epochs=EPOCHS,
    callbacks=[model_checkpoint_callback]
    # To use a validation split, you'd need to create a separate validation dataset
    # validation_data=val_dataset 
)
```

```python
# Inference and Prediction
# Load best weights if you trained with the callback
# transformer.load_weights("viral_transformer_checkpoint.h5")

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

# --- Test with a few examples from our dataset ---
for i in range(5):
    idx = np.random.randint(0, len(df))
    ancestor_seq = df["ancestor_seq"].iloc[idx]
    actual_descendant_seq = df["descendant_seq"].iloc[idx]
    predicted_descendant_seq = decode_sequence(ancestor_seq)
    compare_sequences(ancestor_seq, actual_descendant_seq, predicted_descendant_seq)
```

## Implement the Transformer based method #2 (Conceptual)

This section provides a Python/TensorFlow implementation of the Transformer model.
*This implementation is based on the TensorFlow tutorial "Transformer model for language understanding" ([https://www.tensorflow.org/text/tutorials/transformer](https://www.tensorflow.org/text/tutorials/transformer)), adapted for nucleotide sequences.*

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

## Credits

Jules AI and Gemini 2.5 Pro helped develop this documentation page and the code.
