# Retrieve code corresponding to Hie et al (2021)
!git clone https://github.com/brianhie/viral-mutation

# '%cd' command to change directory
# %cd viral-mutation

### **Step 1: Install Required Libraries**

!pip install -q transformers[torch] biopython

### **Step 2: Load Pre-trained Protein Model**
# Here, we will load a protein language model (ESM-2) from Hugging Face. It is
# specifically trained on protein sequences.

import torch
from transformers import AutoTokenizer, AutoModelForMaskedLM # <-- The key change

# Use a smaller, efficient version of the ESM-2 model
model_name = "facebook/esm2_t6_8M_UR50D"

# Load tokenizer and model with language modeling head
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = AutoModelForMaskedLM.from_pretrained(model_name) # <-- Use AutoModelForMaskedLM

# Move the model to the GPU if available for faster processing
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model.to(device)

print(f"Model '{model_name}' with Language Modeling head loaded successfully on {device}.")

### **Step 3: Load Baseline Sequence**
# Modify this cell to only load single baseline sequence.

import Bio.SeqIO

# Define paths to FASTA files
baseline_file = "/content/viral-mutation/examples/example_wt.fa"
target_file = "/content/viral-mutation/examples/example_target.fa"

# Read single baseline sequence from its file
baseline_record = Bio.SeqIO.read(baseline_file, "fasta")
baseline_seq = str(baseline_record.seq)

print(f"Baseline sequence '{baseline_record.id}' loaded successfully.")

### **Step 4: Process Each Target Variant**
# Loop through target file, generating an embedding for each variant and comparing it to baseline.

import torch

### Step 4: Calculating and Displaying Grammaticality

import torch
import numpy as np

def get_sequence_embedding(sequence, tokenizer, model, device):
    """
    Generates single vector embedding for protein sequence.
    This version is compatible with AutoModelForMaskedLM.
    """
    inputs = tokenizer(sequence, return_tensors="pt", truncation=True, max_length=1024).to(device)
    with torch.no_grad():
        # Ask for the hidden states explicitly.
        outputs = model(**inputs, output_hidden_states=True)

    # Hidden states are in a different attribute.
    # Use last hidden state from the tuple of all hidden states.
    embedding = outputs.hidden_states[-1].mean(dim=1)
    return embedding

def get_sequence_grammaticality(sequence, tokenizer, model, device):
    """Calculates the average log-likelihood of a sequence, a proxy for 'grammaticality'."""
    input_ids = tokenizer.encode(sequence, return_tensors="pt").to(device)
    with torch.no_grad():
        outputs = model(input_ids)
        logits = outputs.logits

    log_probs = torch.nn.functional.log_softmax(logits, dim=-1)
    sequence_log_probs = log_probs[0, :-1].gather(1, input_ids[0, 1:].unsqueeze(-1)).squeeze()
    return sequence_log_probs.mean().item()


# --- Main Analysis Loop ---

# Generate the embedding for the baseline sequence
baseline_embedding = get_sequence_embedding(baseline_seq, tokenizer, model, device)
print(f"Generated embedding for baseline '{baseline_record.id}'.")

# Calculate grammaticality of the baseline sequence
baseline_grammaticality = get_sequence_grammaticality(baseline_seq, tokenizer, model, device)
print(f"Baseline Grammaticality Score: {baseline_grammaticality:.4f}\n")


print("--- Comparing Target Variants to Baseline ---")
# Use Bio.SeqIO.parse again to reset iterator
target_variants = Bio.SeqIO.parse(target_file, "fasta")

# Loop through each variant in the target file
for target_record in target_variants:
    target_seq = str(target_record.seq)

    # 1. Calculate Semantic Change
    target_embedding = get_sequence_embedding(target_seq, tokenizer, model, device)
    cosine_similarity = torch.nn.functional.cosine_similarity(baseline_embedding, target_embedding).item()
    semantic_change = 1 - cosine_similarity

    # 2. Calculate Grammaticality
    target_grammaticality = get_sequence_grammaticality(target_seq, tokenizer, model, device)

    # Print results for this specific variant
    print(f"Variant ID: {target_record.id}")
    print(f"  Semantic Change:   {semantic_change:.6f}")
    print(f"  Grammaticality:    {target_grammaticality:.4f}") # Higher (less negative) is better
    print("-" * 20)

### **Step 1: Define "Past" and "Future" Sequences**

# First, isolate key sequences from the data file: the Delta variant ("past") and
# the Omicron variant ("future" test case).

# Create dictionary to hold key sequences
key_sequences = {}

# Re-parse the target file to find Delta and Omicron
for record in Bio.SeqIO.parse(target_file, "fasta"):
    if record.id == 'spike_delta':
        key_sequences['delta'] = {'id': record.id, 'seq': str(record.seq)}
    elif record.id == 'spike_omicron':
        key_sequences['omicron'] = {'id': record.id, 'seq': str(record.seq)}

if 'delta' in key_sequences and 'omicron' in key_sequences:
    print("Successfully loaded Delta and Omicron sequences.")
    delta_seq = key_sequences['delta']['seq']
    omicron_seq = key_sequences['omicron']['seq']
else:
    print("Error: Could not find Delta or Omicron in the target file.")

### **Step 2: Generate Hypothetical Mutations**

# Create set of random, hypothetical variants based on our Delta sequence. This will be the control group.

import random

def generate_mutants(base_seq, num_mutants, max_mutations):
    """Generates a list of random mutants from a base sequence."""
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    mutants = []

    for i in range(num_mutants):
        mutant_seq = list(base_seq)
        num_to_mutate = random.randint(1, max_mutations)

        for _ in range(num_to_mutate):
            # Pick random position to mutate
            pos = random.randint(0, len(mutant_seq) - 1)
            original_aa = mutant_seq[pos]

            # Pick new amino acid that is different from the original
            new_aa = original_aa
            while new_aa == original_aa:
                new_aa = random.choice(amino_acids)

            mutant_seq[pos] = new_aa

        mutants.append({'id': f'hypothetical_{i}', 'seq': "".join(mutant_seq)})

    return mutants

# Generate 500 hypothetical mutants, each with 1 to 30 random mutations
hypothetical_variants = generate_mutants(delta_seq, num_mutants=500, max_mutations=30)

print(f"Generated {len(hypothetical_variants)} hypothetical variants based on Delta.")

### **Step 3: Analyze All Sequences (Hypothesis Test)**

# Core of the experiment. Calculate two metrics for Omicron and for all 500 hypothetical
# variants, all relative to Delta baseline.

from tqdm.notebook import tqdm # progress bar
import pandas as pd

def get_sequence_embedding(sequence, tokenizer, model, device):
    inputs = tokenizer(sequence, return_tensors="pt", truncation=True, max_length=1024).to(device)
    with torch.no_grad():
        outputs = model(**inputs, output_hidden_states=True)
    return outputs.hidden_states[-1].mean(dim=1)

def get_sequence_grammaticality(sequence, tokenizer, model, device):
    input_ids = tokenizer.encode(sequence, return_tensors="pt").to(device)
    with torch.no_grad():
        outputs = model(input_ids)
        logits = outputs.logits
    log_probs = torch.nn.functional.log_softmax(logits, dim=-1)
    sequence_log_probs = log_probs[0, :-1].gather(1, input_ids[0, 1:].unsqueeze(-1)).squeeze()
    return sequence_log_probs.mean().item()

# --- Run analysis ---
results = []
all_variants_to_test = [key_sequences['omicron']] + hypothetical_variants

# Calculate baseline metrics for Delta
delta_embedding = get_sequence_embedding(delta_seq, tokenizer, model, device)
delta_grammaticality = get_sequence_grammaticality(delta_seq, tokenizer, model, device)

# Add Delta to results for plotting
results.append({
    'id': 'delta',
    'semantic_change': 0.0, # Change from itself is 0
    'grammaticality': delta_grammaticality,
    'type': 'Baseline'
})

# Loop through test set
for variant in tqdm(all_variants_to_test, desc="Analyzing Variants"):
    seq = variant['seq']

    # Calculate metrics
    embedding = get_sequence_embedding(seq, tokenizer, model, device)
    grammaticality = get_sequence_grammaticality(seq, tokenizer, model, device)
    cosine_similarity = torch.nn.functional.cosine_similarity(delta_embedding, embedding).item()
    semantic_change = 1 - cosine_similarity

    results.append({
        'id': variant['id'],
        'semantic_change': semantic_change,
        'grammaticality': grammaticality,
        'type': 'Omicron' if variant['id'] == 'spike_omicron' else 'Hypothetical'
    })

# Convert results to pandas DataFrame for easy plotting
results_df = pd.DataFrame(results)
print("\nAnalysis complete!")
results_df.head()

### **Step 4: Visualize the Results**

# Create 2D scatter plot to verify if Omicron is an outlier.

import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use('seaborn-v0_8-whitegrid')
fig, ax = plt.subplots(figsize=(12, 8))

# Plot hypothetical mutations
sns.scatterplot(
    data=results_df[results_df['type'] == 'Hypothetical'],
    x='semantic_change',
    y='grammaticality',
    alpha=0.5,
    label='Hypothetical Variants',
    ax=ax
)

# Plot Delta baseline
sns.scatterplot(
    data=results_df[results_df['type'] == 'Baseline'],
    x='semantic_change',
    y='grammaticality',
    color='green',
    marker='P', # 'P' for plus sign
    s=200,
    label='Delta (Baseline)',
    ax=ax
)

# Plot Omicron variant
sns.scatterplot(
    data=results_df[results_df['type'] == 'Omicron'],
    x='semantic_change',
    y='grammaticality',
    color='red',
    marker='*',
    s=400,
    label='Omicron (Future Variant)',
    ax=ax
)

ax.set_title('Predictive Analysis: Can We Identify Omicron?', fontsize=16)
ax.set_xlabel('Semantic Change (from Delta)', fontsize=12)
ax.set_ylabel('Grammaticality (Fitness Score)', fontsize=12)
ax.legend()
plt.show()

### **Fitness Heatmap**

```python
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# --- Step 1: Data Preparation ---
# This code assumes data is already in a pandas DataFrame named 'df'.
# If DataFrame has different column names, please update them in the line below.
# For example, if fitness column is named 'score', change 'fitness' to 'score'.

# Pivot the DataFrame to create a matrix suitable for a heatmap.
# Index: Amino Acids ('mutation')
# Columns: Sequence Positions ('position')
# Values: Fitness scores ('fitness')
try:
    heatmap_data = df.pivot(index='mutation', columns='position', values='fitness')

    # Optional: Define a specific order for amino acids on the y-axis if desired.
    # amino_acid_order = list('ARNDCEQGHILKMFPSTWYV')
    # heatmap_data = heatmap_data.reindex(amino_acid_order)

    # --- Step 2: Plotting ---
    plt.figure(figsize=(20, 8)) # Adjust figure size as needed

    # Create heatmap using seaborn
    sns.heatmap(
        heatmap_data,
        cmap='viridis',    # A perceptually uniform colormap is the choice. 'coolwarm' is for divergent data.
        cbar_kws={'label': 'Experimental Fitness Score'} # Label for the color bar
    )

    # --- Step 3: Heatmap ---
    plt.title('Heatmap of Mutation Fitness Landscape', fontsize=16)
    plt.xlabel('Sequence Position', fontsize=12)
    plt.ylabel('Amino Acid Substitution', fontsize=12)
    plt.xticks(rotation=45) # Rotate x-axis labels if overlap
    plt.tight_layout()      # Adjust plot to ensure fit without overlapping
    plt.show()

except KeyError as e:
    print(f"Error: A column was not found. Please check your DataFrame's column names. Missing column: {e}")
except Exception as e:
    print(f"An unexpected error occurred: {e}")