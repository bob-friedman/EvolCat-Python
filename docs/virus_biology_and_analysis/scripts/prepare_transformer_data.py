# prepare_transformer_data.py
#
# This script prepares data for training a transformer model to predict viral evolution.
# It can also generate a multi-FASTA file of descendant sequences.
#
# The script will perform the following steps:
# 1. Parse FASTA files: Read reference genome sequences.
# 2. Parse mutation files: Read mutation data from specified file formats.
# 3. Apply mutations to reference sequence: Generate mutated sequences based on a reference genome and mutation data.
# (Optional) Write descendant sequences to a multi-FASTA file.
# 4. Tokenize sequences: Convert DNA sequences into a numerical format suitable for the model.
# 5. Create input-output pairs: Generate pairs of (original sequence, mutated sequence) for training, now tokenized.
# 6. Pad sequences: Ensure all tokenized sequences have the same length by padding shorter sequences.
# 7. Split data: Divide the dataset into training, validation, and test sets.

import warnings # Keep for now, might be used by dependencies or other parts not touched.
import re
import random
import os
import argparse

# --- Vocabulary Definition ---
NUCLEOTIDE_VOCAB = {'<PAD>': 0, 'A': 1, 'C': 2, 'G': 3, 'T': 4, 'N': 5, '<UNK>': 6}
REV_NUCLEOTIDE_VOCAB = {v: k for k, v in NUCLEOTIDE_VOCAB.items()}

def parse_fasta_file(file_path):
    """
    Parses a FASTA file and returns the genome sequence as a single string.
    Ignores headers and concatenates sequences if multiple entries are present.
    Converts sequence to uppercase.
    """
    sequence = []
    try:
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    continue
                sequence.append(line.strip().upper())
        return "".join(sequence)
    except FileNotFoundError:
        print(f"Error: FASTA file not found at {file_path}")
        return ""

def parse_mutation_file(file_path, mutations_file_delimiter=":"):
    """
    Parses a mutation file to extract identifier and mutation string tuples.
    Returns a list of (identifier, MUTATIONS_string) tuples.
    """
    parsed_data_list = []
    try:
        with open(file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                if mutations_file_delimiter in line:
                    parts = line.split(mutations_file_delimiter, 1)
                    identifier = parts[0].strip()
                    mutations_part = parts[1].strip()

                    parsed_data_list.append((identifier, mutations_part))
                    if not mutations_part:
                        print(f"Info: Empty mutation string (but identifier present) for '{identifier}' on line {line_num} in {file_path}: '{line}'")
                else:
                    print(f"Info: Line {line_num} in {file_path} does not contain '{mutations_file_delimiter}' and was skipped: '{line}'")
        return parsed_data_list
    except FileNotFoundError:
        print(f"Error: Mutations file not found at {file_path}")
        return []

def apply_mutations(ancestral_sequence_str, mutations_str):
    """
    Applies mutations to an ancestral sequence string.
    Validates mutations against the ancestral sequence.
    """
    if not ancestral_sequence_str:
        print("Info: Ancestral sequence is empty. No mutations applied.")
        return ""
    if not mutations_str:
        return ancestral_sequence_str

    sequence_list = list(ancestral_sequence_str)
    mutation_operations = mutations_str.split(',')
    mutation_pattern = re.compile(r"([ACGTNacgtn])(\d+)([ACGTNacgtn])")

    for mut_op in mutation_operations:
        mut_op = mut_op.strip()
        if not mut_op:
            continue

        match = mutation_pattern.fullmatch(mut_op)
        if not match:
            print(f"Info: Invalid mutation format '{mut_op}' in '{mutations_str}'. Skipped.")
            continue

        original_base = match.group(1).upper()
        position_1_indexed = int(match.group(2))
        new_base = match.group(3).upper()
        position_0_indexed = position_1_indexed - 1

        if not (0 <= position_0_indexed < len(sequence_list)):
            print(f"Info: Position {position_1_indexed} in mutation '{mut_op}' (from '{mutations_str}') is out of bounds. Skipped.")
            continue

        if sequence_list[position_0_indexed].upper() != original_base:
            print(f"Info: Mismatch for mutation '{mut_op}' (from '{mutations_str}'). Expected '{original_base}' (at 0-indexed {position_0_indexed}), found '{sequence_list[position_0_indexed]}'. Mutation skipped.")
            continue

        sequence_list[position_0_indexed] = new_base
    return "".join(sequence_list)

def tokenize_sequence(sequence_str, vocab):
    """
    Tokenizes a DNA sequence string using the provided vocabulary.
    Maps unknown characters to a special <UNK> token.
    """
    tokens = []
    unknown_token = vocab.get('<UNK>')

    for nucleotide in sequence_str.upper():
        token = vocab.get(nucleotide, unknown_token)
        if token is unknown_token and nucleotide not in vocab : # Ensure we only warn for truly unknown, not if <UNK> itself is a char
             print(f"Info: Unknown nucleotide '{nucleotide}' found in sequence. Mapped to <UNK>.")
        tokens.append(token)
    return tokens

def create_input_output_pairs(main_reference_sequence_str, identifier_mutation_list, vocab):
    """
    Creates tokenized input-output pairs from a list of (identifier, mutation_str) tuples.
    Returns a list of (tokenized_ancestral, tokenized_descendant, identifier, mutations_str) tuples.
    """
    pairs = []
    if not main_reference_sequence_str:
        print("Info: Main reference sequence is empty for create_input_output_pairs. No pairs created.")
        return pairs

    tokenized_main_reference = tokenize_sequence(main_reference_sequence_str, vocab)

    for i, (identifier, mutations_str) in enumerate(identifier_mutation_list):
        descendant_sequence_str = apply_mutations(main_reference_sequence_str, mutations_str)
        tokenized_descendant = tokenize_sequence(descendant_sequence_str, vocab)
        pairs.append((tokenized_main_reference, tokenized_descendant, identifier, mutations_str))
    return pairs

def pad_sequences(tokenized_pairs, max_len, padding_token=0):
    """
    Pads or truncates tokenized sequences in pairs to max_len.
    Input: list of (tokenized_ancestral, tokenized_descendant, identifier, mutations_str)
    Output: list of (padded_ancestral, padded_descendant, identifier, mutations_str)
    """
    padded_pairs = []
    for tok_anc, tok_desc, identifier, mut_str in tokenized_pairs:
        padded_anc = tok_anc[:max_len] + [padding_token] * max(0, max_len - len(tok_anc))
        padded_desc = tok_desc[:max_len] + [padding_token] * max(0, max_len - len(tok_desc))
        padded_pairs.append((padded_anc, padded_desc, identifier, mut_str))
    return padded_pairs

def split_data(all_padded_pairs, train_ratio=0.8, val_ratio=0.1, random_seed=None):
    """
    Splits a list of padded pairs into training, validation, and test sets.
    """
    if not (0 <= train_ratio <= 1 and 0 <= val_ratio <= 1 and (train_ratio + val_ratio) <= 1):
        raise ValueError("Train and validation ratios must be between 0 and 1, and their sum must be <= 1.")

    shuffled_pairs = list(all_padded_pairs)
    if random_seed is not None:
        random.seed(random_seed)
    random.shuffle(shuffled_pairs)

    n_total = len(shuffled_pairs)
    n_train = int(n_total * train_ratio)
    n_val = int(n_total * val_ratio)

    train_set = shuffled_pairs[:n_train]
    val_set = shuffled_pairs[n_train : n_train + n_val]
    test_set = shuffled_pairs[n_train + n_val:]

    return train_set, val_set, test_set

def write_multifasta(filepath, sequences_data):
    """
    Writes a list of sequences to a single multi-FASTA file.

    Args:
        filepath (str): The path to the output multi-FASTA file.
        sequences_data (list): A list of tuples, where each tuple is
                               (header_string, sequence_string).
    """
    try:
        with open(filepath, 'w') as f:
            for header, sequence in sequences_data:
                f.write(f">{header}\n")
                f.write(f"{sequence}\n")
        print(f"Info: Successfully wrote {len(sequences_data)} sequences to {filepath}")
    except IOError as e:
        print(f"Error: Could not write to file {filepath}. Reason: {e}")

def main():
    parser = argparse.ArgumentParser(description="Prepare data for Transformer model and optionally generate FASTA files.")
    parser.add_argument('--reference_fasta', type=str, default="dummy_reference.fasta", help='Path to the reference FASTA file.')
    parser.add_argument('--mutations_file', type=str, default="dummy_mutations.txt", help='Path to the mutations file.')
    parser.add_argument('--mutations_delimiter', type=str, default=":", help='Delimiter used in the mutations file.')
    parser.add_argument('--output_multifasta_file', type=str, default="output_fasta_sequences/descendant_sequences.fasta", help='Path to save the output multi-FASTA file. If empty, FASTA generation is skipped.')

    args = parser.parse_args()

    REFERENCE_FASTA_PATH = args.reference_fasta
    MUTATIONS_FILE_PATH = args.mutations_file
    mutations_file_delimiter = args.mutations_delimiter
    output_fasta_filepath = args.output_multifasta_file

    print("Starting data preparation script with arguments:")
    print(f"  Reference FASTA: {REFERENCE_FASTA_PATH}")
    print(f"  Mutations File: {MUTATIONS_FILE_PATH}")
    print(f"  Mutations Delimiter: '{mutations_file_delimiter}'")
    print(f"  Output Multi-FASTA: {output_fasta_filepath if output_fasta_filepath else 'Not generating'}")

    # Create dummy files if defaults are used and files don't exist
    if REFERENCE_FASTA_PATH == "dummy_reference.fasta" and not os.path.exists(REFERENCE_FASTA_PATH):
        with open(REFERENCE_FASTA_PATH, 'w') as f:
            f.write(">dummy_ref_genome\nAAAAAGGGGGTTTTTCCCCCNNNNNX") # Length 26
        print(f"Info: Created dummy reference FASTA: {REFERENCE_FASTA_PATH}")

    if MUTATIONS_FILE_PATH == "dummy_mutations.txt" and not os.path.exists(MUTATIONS_FILE_PATH):
        with open(MUTATIONS_FILE_PATH, 'w') as f:
            f.write("node1:G6T,C18G\nnode2:A1C,G10N,N22A\nnode3:X5Y,G7A\nnode4:C6T\nnode5:N21A,N22T,N23G,N24C,N25A\nnode6:G26X\nnode7:\nnode8:A30T\nnode9:A2T\nnode10:C3G")
        print(f"Info: Created dummy mutations file: {MUTATIONS_FILE_PATH}")

    # Create output directory for FASTA if specified and doesn't exist
    if output_fasta_filepath:
        output_fasta_dir = os.path.dirname(output_fasta_filepath)
        if output_fasta_dir and not os.path.exists(output_fasta_dir):
            os.makedirs(output_fasta_dir, exist_ok=True)
            print(f"Info: Created output directory for FASTA: {output_fasta_dir}")

    print(f"\n--- Parsing FASTA file: {REFERENCE_FASTA_PATH} ---")
    reference_genome = parse_fasta_file(REFERENCE_FASTA_PATH)
    if not reference_genome:
        print(f"Error: Could not load reference genome from {REFERENCE_FASTA_PATH}. Exiting.")
        return # Changed from exit(1) to return for easier testing if needed

    print(f"Loaded reference genome. Length: {len(reference_genome)}")
    print(f"Sequence (first 60 chars): {reference_genome[:60]}\n")

    max_len = len(reference_genome)
    print(f"Using max_len for padding/analysis: {max_len}")

    print(f"\n--- Parsing Mutation file: {MUTATIONS_FILE_PATH} ---")
    identifier_mutation_list = parse_mutation_file(MUTATIONS_FILE_PATH, mutations_file_delimiter=mutations_file_delimiter)

    if not identifier_mutation_list:
        print("Info: No mutation data loaded or file problem. FASTA generation will only contain reference if specified. Transformer data pipeline might be empty.")
        # Allow to proceed to write reference if specified, or generate empty transformer data
    else:
        print(f"Loaded {len(identifier_mutation_list)} mutation entries (identifier, mutations_str).")
        if identifier_mutation_list:
            print(f"  Example entry: Identifier='{identifier_mutation_list[0][0]}', Mutations='{identifier_mutation_list[0][1]}'")

    # --- Generate and Write Descendant Sequences for FASTA output ---
    if output_fasta_filepath: # Only proceed if an output path is given
        descendant_sequences_for_fasta = []
        if reference_genome: # Need reference_genome to proceed
            print(f"\n--- Generating descendant sequences for FASTA output... ---")
            # Add reference itself to the list
            descendant_sequences_for_fasta.append((os.path.basename(REFERENCE_FASTA_PATH) + "_reference", reference_genome))

            if identifier_mutation_list:
                for identifier, mutation_str in identifier_mutation_list:
                    if not mutation_str:
                        descendant_seq = reference_genome
                        print(f"Info: No mutations for identifier '{identifier}', using reference sequence for FASTA entry.")
                    else:
                        descendant_seq = apply_mutations(reference_genome, mutation_str)
                    descendant_sequences_for_fasta.append((identifier, descendant_seq))

            if descendant_sequences_for_fasta:
                write_multifasta(output_fasta_filepath, descendant_sequences_for_fasta)
            else: # Should not happen if reference_genome exists
                print(f"Info: No descendant sequences (including reference) to write. FASTA file '{output_fasta_filepath}' will not be created or will be empty.")
        else:
            print("Info: Reference genome not loaded. Cannot generate descendant sequences for FASTA.")


    # --- Data Preparation for Transformer Model ---
    print(f"\n--- Creating Tokenized Input-Output Pairs for Transformer ---")
    if not identifier_mutation_list and reference_genome:
         # If no mutations, create one pair: (ref, ref) for model training to have some data
        print("Info: No mutations found, creating a self-reference pair for transformer data.")
        identifier_mutation_list_for_transformer = [("reference_self_pair", "")]
    elif not reference_genome:
        print("Error: Reference genome not available. Cannot create pairs for transformer. Exiting data prep for transformer.")
        return
    else:
        identifier_mutation_list_for_transformer = identifier_mutation_list

    tokenized_input_output_pairs = create_input_output_pairs(reference_genome, identifier_mutation_list_for_transformer, NUCLEOTIDE_VOCAB)

    if not tokenized_input_output_pairs:
        print("No tokenized pairs created for transformer. Exiting data prep for transformer.")
        return # Stop if no pairs, e.g. if reference_genome was also empty

    print(f"Created {len(tokenized_input_output_pairs)} tokenized input-output pairs for Transformer.")
    if tokenized_input_output_pairs:
         print(f"  Example pair content (first pair): Ancestral tokens ({len(tokenized_input_output_pairs[0][0])}), Descendant tokens ({len(tokenized_input_output_pairs[0][1])}), Identifier ('{tokenized_input_output_pairs[0][2]}'), Mutations ('{tokenized_input_output_pairs[0][3]}')")


    print(f"\n--- Padding Sequences for Transformer ---")
    padded_pairs = pad_sequences(tokenized_input_output_pairs, max_len, NUCLEOTIDE_VOCAB['<PAD>'])
    print(f"Created {len(padded_pairs)} padded pairs for Transformer.")

    print(f"\n--- Splitting Data for Transformer ---")
    if padded_pairs:
        train_data, val_data, test_data = split_data(padded_pairs, train_ratio=0.8, val_ratio=0.1, random_seed=42)
        print(f"Split data into:")
        print(f"  Training set: {len(train_data)} samples")
        print(f"  Validation set: {len(val_data)} samples")
        print(f"  Test set: {len(test_data)} samples")

        if train_data:
            print(f"\nExample: First item in training_data (identifier: '{train_data[0][2]}', mutations: '{train_data[0][3]}')")
    else:
        print("Skipping data splitting as no padded pairs were available.")

    print("\nData preparation script finished.")

if __name__ == "__main__":
    main()
```
