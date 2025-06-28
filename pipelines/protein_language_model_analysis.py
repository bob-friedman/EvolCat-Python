"""
Protein Language Model Analysis Script

This script analyzes protein sequences using a pre-trained protein language model (ESM-2).
It performs the following main tasks:
1. Loads a baseline protein sequence and target variant sequences from FASTA files.
2. Calculates semantic change and grammaticality for each target variant relative to the baseline.
3. Identifies specific "past" (e.g., Delta) and "future" (e.g., Omicron) variants.
4. Generates a set of hypothetical mutants based on the "past" variant.
5. Compares the "future" variant against the "past" variant and the hypothetical mutants
   based on semantic change and grammaticality.
6. Visualizes these comparisons using a scatter plot.
7. Includes a section for generating a fitness heatmap (requires a DataFrame 'df' with
   mutation, position, and fitness data).

Workflow:
- Initialize Model and Tokenizer: Loads the ESM-2 model and tokenizer.
- Load Data: Reads baseline and target sequences.
- Calculate Metrics:
    - `get_sequence_embedding`: Generates a vector embedding for a sequence.
    - `get_sequence_grammaticality`: Calculates the log-likelihood of a sequence.
- Compare Variants: Computes semantic change and grammaticality for target variants.
- Hypothesis Test:
    - Defines past/future sequences (e.g., Delta/Omicron).
    - `generate_mutants`: Creates random mutations from the past sequence.
    - Analyzes future and hypothetical variants against the past sequence.
- Visualize: Plots semantic change vs. grammaticality.
- Fitness Heatmap: Pivots fitness data and plots a heatmap.

Input:
- Baseline FASTA file (e.g., 'example_wt.fa')
- Target variants FASTA file (e.g., 'example_target.fa')
- For heatmap: A pandas DataFrame `df` with 'mutation', 'position', 'fitness' columns.

Output:
- Console output: Metrics for each variant, analysis progress.
- Scatter plot: Visualizing Omicron vs. hypothetical variants.
- Heatmap plot (if fitness data is provided).

Dependencies:
- torch
- transformers
- biopython
- numpy
- pandas
- matplotlib
- seaborn
- tqdm (for progress bars)
- GitPython (if cloning repository from script)

Data Source and Setup:
The example data and some conceptual parts of this script are inspired by or may use resources
from the `viral-mutation` project by Brian Hie: https://github.com/brianhie/viral-mutation.
This repository contains data and models relevant to viral evolution and protein engineering.

To use data from this repository, you would typically clone it first.
You can clone it manually using git:
  git clone https://github.com/brianhie/viral-mutation.git

Or, you can clone it from within a Python script using the `GitPython` library:
  import git
  try:
    git.Repo.clone_from("https://github.com/brianhie/viral-mutation.git", "viral-mutation")
    print("Repository cloned successfully.")
  except git.exc.GitCommandError as e:
    print(f"Error cloning repository: {e}")
    print("It might already exist, or there might be other issues.")

The FASTA file paths (e.g., `DEFAULT_BASELINE_FILE`, `DEFAULT_TARGET_FILE`) might need
to be adjusted to point to the correct locations within the cloned `viral-mutation`
repository, or your own data files. For example, if `viral-mutation` is cloned into
the same directory as this script, paths might be like:
`os.path.join("viral-mutation", "examples", "example_wt.fa")`

Note: The original script included shell commands for `git clone` and `pip install`.
These should be handled outside of this script (e.g., in a setup script or manually
if not using the Python cloning method described above).
File paths for FASTA files are placeholders and may need to be adjusted.
The heatmap generation part assumes a DataFrame `df` is already loaded with specific columns.
"""

import torch
from transformers import AutoTokenizer, AutoModelForMaskedLM
import Bio.SeqIO
import numpy as np
import random
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm # Use standard tqdm, not tqdm.notebook
import argparse
import os

# --- Configuration ---
DEFAULT_MODEL_NAME = "facebook/esm2_t6_8M_UR50D"
# Placeholder file paths - these should be configurable or determined dynamically
# DEFAULT_BASELINE_FILE = "/content/viral-mutation/examples/example_wt.fa"
# DEFAULT_TARGET_FILE = "/content/viral-mutation/examples/example_target.fa"
# For the script to run without the /content/viral-mutation path,
# we assume the example files are in a subdirectory 'examples' relative to this script,
# or specified via command-line arguments.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_BASELINE_FILE = os.path.join(SCRIPT_DIR, "examples", "example_wt.fa")
DEFAULT_TARGET_FILE = os.path.join(SCRIPT_DIR, "examples", "example_target.fa")


# --- Core Functions ---

def load_protein_model(model_name: str = DEFAULT_MODEL_NAME) -> tuple:
    """
    Loads a pre-trained protein language model and tokenizer from Hugging Face.
    Moves the model to GPU if available.
    """
    print(f"Loading model '{model_name}'...")
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModelForMaskedLM.from_pretrained(model_name)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)
    print(f"Model '{model_name}' with Language Modeling head loaded successfully on {device}.")
    return tokenizer, model, device

def load_fasta_sequence(file_path: str, is_single_record: bool = True) -> any:
    """
    Loads sequence(s) from a FASTA file.
    If is_single_record is True, expects one sequence and returns a SeqRecord.
    Otherwise, returns an iterator of SeqRecords.
    """
    if not os.path.exists(file_path):
        print(f"Error: FASTA file not found at {file_path}")
        return None
    if is_single_record:
        try:
            record = Bio.SeqIO.read(file_path, "fasta")
            print(f"Sequence '{record.id}' loaded successfully from {file_path}.")
            return record
        except ValueError:
            print(f"Error: More than one sequence found in {file_path}, expected single record.")
            return None
    else:
        return Bio.SeqIO.parse(file_path, "fasta")


def get_sequence_embedding(sequence: str, tokenizer, model, device) -> torch.Tensor:
    """
    Generates a single vector embedding for a protein sequence using the specified model.
    Compatible with AutoModelForMaskedLM.
    """
    inputs = tokenizer(sequence, return_tensors="pt", truncation=True, max_length=1024).to(device)
    with torch.no_grad():
        outputs = model(**inputs, output_hidden_states=True)
    # Use the last hidden state, averaged across sequence length
    embedding = outputs.hidden_states[-1].mean(dim=1)
    return embedding

def get_sequence_grammaticality(sequence: str, tokenizer, model, device) -> float:
    """
    Calculates the average log-likelihood of a sequence, a proxy for 'grammaticality'.
    Higher (less negative) scores indicate better grammaticality.
    """
    input_ids = tokenizer.encode(sequence, return_tensors="pt").to(device)
    if input_ids.shape[1] == 0: # Handle empty sequences if they occur
        return -float('inf')
    with torch.no_grad():
        outputs = model(input_ids, labels=input_ids.clone()) # Provide labels for MLM
        # Logits are for MLM, loss is already calculated if labels are provided
        # For grammaticality, we often use pseudo-log-likelihoods if not directly using MLM loss
        # The original script's approach:
        logits = outputs.logits

    log_probs = torch.nn.functional.log_softmax(logits, dim=-1)
    # Exclude start/end tokens if tokenizer adds them; ensure input_ids correspond to actual sequence tokens
    # The original code slices input_ids[0, 1:], assuming a start token.
    # For ESM, this might depend on tokenizer settings.
    # Let's ensure the slicing is safe.
    if input_ids.shape[1] <= 1: # Sequence too short
        return -float('inf')

    sequence_log_probs = log_probs[0, :-1].gather(1, input_ids[0, 1:].unsqueeze(-1)).squeeze()
    return sequence_log_probs.mean().item()

def compare_variants_to_baseline(baseline_seq_record, target_file_path: str, tokenizer, model, device):
    """
    Processes each target variant, calculating and printing its semantic change
    and grammaticality compared to the baseline sequence.
    """
    if baseline_seq_record is None:
        print("Baseline sequence not loaded. Skipping variant comparison.")
        return

    baseline_seq = str(baseline_seq_record.seq)
    baseline_embedding = get_sequence_embedding(baseline_seq, tokenizer, model, device)
    baseline_grammaticality = get_sequence_grammaticality(baseline_seq, tokenizer, model, device)

    print(f"\nBaseline ID: {baseline_seq_record.id}")
    print(f"Baseline Grammaticality Score: {baseline_grammaticality:.4f}\n")
    print("--- Comparing Target Variants to Baseline ---")

    target_variants_iterator = load_fasta_sequence(target_file_path, is_single_record=False)
    if target_variants_iterator is None:
        print(f"Could not load target variants from {target_file_path}. Skipping comparison.")
        return

    for target_record in target_variants_iterator:
        target_seq = str(target_record.seq)
        if not target_seq:
            print(f"Warning: Empty sequence for variant ID: {target_record.id}. Skipping.")
            continue

        # 1. Calculate Semantic Change
        target_embedding = get_sequence_embedding(target_seq, tokenizer, model, device)
        cosine_similarity = torch.nn.functional.cosine_similarity(baseline_embedding, target_embedding).item()
        semantic_change = 1 - cosine_similarity

        # 2. Calculate Grammaticality
        target_grammaticality = get_sequence_grammaticality(target_seq, tokenizer, model, device)

        print(f"Variant ID: {target_record.id}")
        print(f"  Semantic Change (from baseline): {semantic_change:.6f}")
        print(f"  Grammaticality:                  {target_grammaticality:.4f}")
        print("-" * 30)

def define_past_future_sequences(target_file_path: str, past_id: str = 'spike_delta', future_id: str = 'spike_omicron') -> dict:
    """
    Isolates key sequences (e.g., Delta and Omicron) from the target FASTA file.
    Returns a dictionary containing these sequences.
    """
    key_sequences = {}
    target_variants_iterator = load_fasta_sequence(target_file_path, is_single_record=False)
    if target_variants_iterator is None:
        print(f"Could not load target variants from {target_file_path} for past/future definition.")
        return key_sequences # Return empty dict

    for record in target_variants_iterator:
        if record.id == past_id:
            key_sequences['past'] = {'id': record.id, 'seq': str(record.seq)}
        elif record.id == future_id:
            key_sequences['future'] = {'id': record.id, 'seq': str(record.seq)}

    if 'past' in key_sequences:
        print(f"Successfully loaded '{past_id}' as 'past' sequence.")
    else:
        print(f"Error: Could not find '{past_id}' in {target_file_path}.")

    if 'future' in key_sequences:
        print(f"Successfully loaded '{future_id}' as 'future' sequence.")
    else:
        print(f"Error: Could not find '{future_id}' in {target_file_path}.")

    return key_sequences

def generate_mutants(base_seq: str, num_mutants: int, max_mutations: int) -> list:
    """
    Generates a list of random mutants from a base sequence.
    Each mutant has a random number of mutations up to max_mutations.
    """
    if not base_seq:
        print("Error: Base sequence for mutant generation is empty.")
        return []
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY' # Standard amino acids
    mutants = []

    for i in range(num_mutants):
        mutant_seq_list = list(base_seq)
        num_to_mutate = random.randint(1, max_mutations)

        positions_to_mutate = random.sample(range(len(mutant_seq_list)), min(num_to_mutate, len(mutant_seq_list)))

        for pos in positions_to_mutate:
            original_aa = mutant_seq_list[pos]
            new_aa = random.choice([aa for aa in amino_acids if aa != original_aa]) # Ensure different AA
            mutant_seq_list[pos] = new_aa

        mutants.append({'id': f'hypothetical_{i}', 'seq': "".join(mutant_seq_list)})
    print(f"Generated {len(mutants)} hypothetical variants based on the 'past' sequence.")
    return mutants

def analyze_hypothesis_sequences(key_sequences: dict, hypothetical_variants: list, tokenizer, model, device) -> pd.DataFrame:
    """
    Analyzes the 'future' variant and hypothetical variants relative to the 'past' variant.
    Calculates semantic change and grammaticality for all.
    Returns a pandas DataFrame with the results.
    """
    results = []
    if 'past' not in key_sequences or not key_sequences['past']['seq']:
        print("Error: 'Past' sequence not available for hypothesis analysis.")
        return pd.DataFrame()

    past_seq = key_sequences['past']['seq']
    past_id = key_sequences['past']['id']

    # Calculate baseline metrics for the 'past' sequence
    past_embedding = get_sequence_embedding(past_seq, tokenizer, model, device)
    past_grammaticality = get_sequence_grammaticality(past_seq, tokenizer, model, device)

    results.append({
        'id': past_id,
        'semantic_change': 0.0, # Change from itself is 0
        'grammaticality': past_grammaticality,
        'type': 'Baseline (Past)'
    })

    all_variants_to_test = []
    if 'future' in key_sequences and key_sequences['future']['seq']:
        all_variants_to_test.append(key_sequences['future'])
    else:
        print("Warning: 'Future' sequence not available for hypothesis analysis.")

    all_variants_to_test.extend(hypothetical_variants)

    if not all_variants_to_test:
        print("No variants to test in hypothesis analysis.")
        return pd.DataFrame(results) # Return df with only past/baseline

    print(f"\n--- Analyzing 'Future' and Hypothetical Variants (relative to '{past_id}') ---")
    for variant in tqdm(all_variants_to_test, desc="Analyzing Variants for Hypothesis Test"):
        seq = variant['seq']
        if not seq:
            print(f"Warning: Empty sequence for variant ID: {variant['id']}. Skipping.")
            continue

        embedding = get_sequence_embedding(seq, tokenizer, model, device)
        grammaticality = get_sequence_grammaticality(seq, tokenizer, model, device)
        cosine_similarity = torch.nn.functional.cosine_similarity(past_embedding, embedding).item()
        semantic_change = 1 - cosine_similarity

        variant_type = 'Future' if 'future' in key_sequences and variant['id'] == key_sequences['future']['id'] else 'Hypothetical'

        results.append({
            'id': variant['id'],
            'semantic_change': semantic_change,
            'grammaticality': grammaticality,
            'type': variant_type
        })

    results_df = pd.DataFrame(results)
    print("\nHypothesis analysis complete!")
    print(results_df.head())
    return results_df

def visualize_hypothesis_results(results_df: pd.DataFrame, output_plot_file: str = "predictive_analysis_plot.png"):
    """
    Creates and saves a 2D scatter plot of semantic change vs. grammaticality.
    Highlights Baseline (Past), Future, and Hypothetical variants.
    """
    if results_df.empty:
        print("No data to visualize for hypothesis results.")
        return

    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(12, 8))

    # Plot hypothetical mutations
    hypothetical_data = results_df[results_df['type'] == 'Hypothetical']
    if not hypothetical_data.empty:
        sns.scatterplot(
            data=hypothetical_data,
            x='semantic_change',
            y='grammaticality',
            alpha=0.5,
            label='Hypothetical Variants',
            ax=ax
        )

    # Plot Baseline (Past)
    baseline_data = results_df[results_df['type'] == 'Baseline (Past)']
    if not baseline_data.empty:
        sns.scatterplot(
            data=baseline_data,
            x='semantic_change',
            y='grammaticality',
            color='green',
            marker='P', # 'P' for plus sign
            s=200,
            label='Baseline (Past)',
            ax=ax
        )

    # Plot Future variant
    future_data = results_df[results_df['type'] == 'Future']
    if not future_data.empty:
        sns.scatterplot(
            data=future_data,
            x='semantic_change',
            y='grammaticality',
            color='red',
            marker='*',
            s=400,
            label='Future Variant',
            ax=ax
        )
    else:
        print("Note: 'Future' variant data not found for plotting.")


    ax.set_title('Predictive Analysis: Semantic Change vs. Grammaticality', fontsize=16)
    ax.set_xlabel('Semantic Change (from Baseline/Past)', fontsize=12)
    ax.set_ylabel('Grammaticality (Fitness Score)', fontsize=12)
    ax.legend()
    plt.savefig(output_plot_file)
    print(f"Predictive analysis plot saved to {output_plot_file}")
    # plt.show() # Comment out if running in a non-interactive environment

def plot_fitness_heatmap(df: pd.DataFrame, output_heatmap_file: str = "fitness_heatmap.png"):
    """
    Generates and saves a heatmap of mutation fitness landscape.
    Assumes 'df' has 'mutation', 'position', and 'fitness' columns.
    """
    print("\n--- Generating Fitness Heatmap ---")
    if df is None or df.empty:
        print("DataFrame for heatmap is empty or not provided. Skipping heatmap generation.")
        return
    if not all(col in df.columns for col in ['mutation', 'position', 'fitness']):
        print("Error: DataFrame for heatmap must contain 'mutation', 'position', and 'fitness' columns.")
        print(f"DataFrame columns found: {df.columns.tolist()}")
        return

    try:
        heatmap_data = df.pivot(index='mutation', columns='position', values='fitness')

        # Optional: Define a specific order for amino acids on the y-axis if desired.
        # amino_acid_order = list('ACDEFGHIKLMNPQRSTVWY') # Standard order
        # heatmap_data = heatmap_data.reindex(amino_acid_order).dropna(how='all')


        plt.figure(figsize=(20, 8))
        sns.heatmap(
            heatmap_data,
            cmap='viridis',
            cbar_kws={'label': 'Experimental Fitness Score'},
            annot=False # Annotations can make it too busy for large heatmaps
        )
        plt.title('Heatmap of Mutation Fitness Landscape', fontsize=16)
        plt.xlabel('Sequence Position', fontsize=12)
        plt.ylabel('Amino Acid Substitution', fontsize=12)
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        plt.savefig(output_heatmap_file)
        print(f"Fitness heatmap saved to {output_heatmap_file}")
        # plt.show() # Comment out if running in a non-interactive environment

    except KeyError as e:
        print(f"Error creating heatmap: A required column might be missing or misnamed. Missing column: {e}")
    except Exception as e:
        print(f"An unexpected error occurred during heatmap generation: {e}")


# --- Main Execution ---

def main():
    parser = argparse.ArgumentParser(description="Protein Language Model Analysis Script")
    parser.add_argument("--model_name", type=str, default=DEFAULT_MODEL_NAME,
                        help="Name of the pre-trained model from Hugging Face (e.g., 'facebook/esm2_t6_8M_UR50D')")
    parser.add_argument("--baseline_file", type=str, default=DEFAULT_BASELINE_FILE,
                        help="Path to the FASTA file containing the baseline sequence.")
    parser.add_argument("--target_file", type=str, default=DEFAULT_TARGET_FILE,
                        help="Path to the FASTA file containing target variant sequences.")
    parser.add_argument("--past_variant_id", type=str, default="spike_delta",
                        help="ID of the 'past' variant in the target file (e.g., 'spike_delta').")
    parser.add_argument("--future_variant_id", type=str, default="spike_omicron",
                        help="ID of the 'future' variant in the target file (e.g., 'spike_omicron').")
    parser.add_argument("--num_hypothetical_mutants", type=int, default=500,
                        help="Number of hypothetical mutants to generate.")
    parser.add_argument("--max_mutations_per_mutant", type=int, default=30,
                        help="Maximum number of mutations in each hypothetical mutant.")
    parser.add_argument("--skip_variant_comparison", action="store_true",
                        help="Skip the initial comparison of all target variants to baseline.")
    parser.add_argument("--skip_hypothesis_test", action="store_true",
                        help="Skip the hypothesis test involving past, future, and hypothetical variants.")
    parser.add_argument("--skip_heatmap", action="store_true",
                        help="Skip generating the fitness heatmap (requires a specific 'df').")
    parser.add_argument("--output_plot_file", type=str, default="predictive_analysis_plot.png",
                        help="Filename for the output predictive analysis scatter plot.")
    parser.add_argument("--output_heatmap_file", type=str, default="fitness_heatmap.png",
                        help="Filename for the output fitness heatmap.")

    args = parser.parse_args()

    # --- Step 1: Load Pre-trained Protein Model ---
    tokenizer, model, device = load_protein_model(args.model_name)

    # --- Step 2: Load Baseline Sequence ---
    baseline_seq_record = load_fasta_sequence(args.baseline_file, is_single_record=True)
    if baseline_seq_record is None:
        print(f"Critical error: Baseline sequence from {args.baseline_file} could not be loaded. Exiting.")
        return

    # --- Step 3: Process Each Target Variant (Optional) ---
    if not args.skip_variant_comparison:
        print("\n### Part 1: Comparing all target variants to baseline ###")
        compare_variants_to_baseline(baseline_seq_record, args.target_file, tokenizer, model, device)

    if not args.skip_hypothesis_test:
        print("\n### Part 2: Hypothesis Test (Past vs. Future vs. Hypothetical) ###")
        # --- Step 3a: Define "Past" and "Future" Sequences ---
        key_sequences = define_past_future_sequences(args.target_file, args.past_variant_id, args.future_variant_id)

        if 'past' not in key_sequences or not key_sequences['past']['seq']:
            print("Critical error: 'Past' sequence could not be loaded or is empty. Skipping hypothesis test.")
        else:
            past_seq = key_sequences['past']['seq']

            # --- Step 3b: Generate Hypothetical Mutations ---
            hypothetical_variants = generate_mutants(past_seq, args.num_hypothetical_mutants, args.max_mutations_per_mutant)

            # --- Step 3c: Analyze All Sequences for Hypothesis Test ---
            results_df = analyze_hypothesis_sequences(key_sequences, hypothetical_variants, tokenizer, model, device)

            # --- Step 3d: Visualize the Hypothesis Results ---
            if not results_df.empty:
                visualize_hypothesis_results(results_df, args.output_plot_file)
            else:
                print("Skipping visualization as hypothesis analysis yielded no results.")
    else:
        print("\nSkipping hypothesis test as per user request.")

    # --- Step 4: Fitness Heatmap (Example, requires specific DataFrame 'df') ---
    if not args.skip_heatmap:
        print("\n### Part 3: Fitness Heatmap ###")
        # This part is illustrative as the original script implies 'df' comes from elsewhere.
        # For demonstration, we create a dummy df. Replace with actual data loading if needed.
        print("Note: The fitness heatmap requires a DataFrame 'df' with 'mutation', 'position', and 'fitness' columns.")
        print("Attempting to generate heatmap with a dummy DataFrame if no actual df is loaded.")

        # Example of how 'df' might be structured:
        example_fitness_data = {
            'mutation': ['A1C', 'A1D', 'A1E', 'C2A', 'C2G', 'Y171N', 'Y171C'],
            'position': [1, 1, 1, 2, 2, 171, 171],
            'fitness': [0.1, 0.5, 0.2, 0.8, 0.3, 0.9, 0.4],
            # The 'mutation' column in the original script's heatmap seems to be the target amino acid.
            # Let's adjust to make it compatible with the pivot:
            'amino_acid_substitution': ['C', 'D', 'E', 'A', 'G', 'N', 'C']

        }
        # df_heatmap = pd.DataFrame(example_fitness_data)
        # df_heatmap['mutation'] = df_heatmap['amino_acid_substitution'] # Use this for index

        # A more direct interpretation of the pivot df.pivot(index='mutation', columns='position', values='fitness')
        # is that 'mutation' is the *new amino acid* at a given position.
        # Let's create a dummy df that fits this structure for the heatmap function to run.
        positions = list(range(1, 5)) # Example positions
        amino_acids = list('ACDE')    # Example amino acids
        example_data_for_heatmap = []
        for pos in positions:
            for aa in amino_acids:
                example_data_for_heatmap.append({
                    'position': pos,
                    'mutation': aa, # This is the new amino acid
                    'fitness': random.uniform(0,1)
                })
        df_for_heatmap = pd.DataFrame(example_data_for_heatmap)

        plot_fitness_heatmap(df_for_heatmap, args.output_heatmap_file)
        print("If you have actual fitness data in a CSV, you might load it like: ")
        print("# df_actual = pd.read_csv('your_fitness_data.csv')")
        print("# plot_fitness_heatmap(df_actual, args.output_heatmap_file)")
    else:
        print("\nSkipping fitness heatmap generation as per user request.")

    print("\nScript execution finished.")

if __name__ == "__main__":
    # Create dummy example files if they don't exist, so the script can run.
    # This is for demonstration purposes. In a real scenario, these files would be inputs.
    if not os.path.exists(DEFAULT_BASELINE_FILE):
        os.makedirs(os.path.dirname(DEFAULT_BASELINE_FILE), exist_ok=True)
        with open(DEFAULT_BASELINE_FILE, "w") as f:
            f.write(">example_wt\n")
            f.write("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPTSCCSWEQRMTVVCYSCGTEYPTMTSETVEFAYGTADKPYQGRT\n") # Short example sequence
    if not os.path.exists(DEFAULT_TARGET_FILE):
        os.makedirs(os.path.dirname(DEFAULT_TARGET_FILE), exist_ok=True)
        with open(DEFAULT_TARGET_FILE, "w") as f:
            f.write(">spike_delta\n")
            f.write("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPTSCCSWEQRMTVVCYSCGTEYPTMTSETVEFAYGTADKPYQGRT\n") # Same as wt for simplicity
            f.write(">spike_omicron\n")
            f.write("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPTSCCSWEQRMTVVCYSCGTEYPTMTSETVEFAYGTADKPYQGRS\n") # One mutation
            f.write(">another_variant\n")
            f.write("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPTSCCSWEQRMTVVCYSCGTEYPTMTSETVEFAYGTTDKPYQGRT\n") # Different mutation

    main()
