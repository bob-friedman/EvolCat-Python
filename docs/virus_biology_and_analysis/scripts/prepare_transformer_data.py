# prepare_transformer_data.py
#
# This script prepares data for training a transformer model to predict viral evolution.
#
# The script will perform the following steps:
# 1. Parse FASTA files: Read reference genome sequences.
# 2. Parse mutation files: Read mutation data from specified file formats.
# 3. Apply mutations to reference sequence: Generate mutated sequences based on a reference genome and mutation data.
# 4. Tokenize sequences: Convert DNA sequences into a numerical format suitable for the model.
# 5. Create input-output pairs: Generate pairs of (original sequence, mutated sequence) for training, now tokenized.
# 6. Pad sequences: Ensure all tokenized sequences have the same length by padding shorter sequences.
# 7. Split data: Divide the dataset into training, validation, and test sets.

import warnings
import re
import random

# --- Vocabulary Definition ---
NUCLEOTIDE_VOCAB = {'<PAD>': 0, 'A': 1, 'C': 2, 'G': 3, 'T': 4, 'N': 5, '<UNK>': 6} # Added <UNK> for unknown chars
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

def parse_mutation_file(file_path):
    """
    Parses a mutation file to extract mutation strings.
    Returns a list of MUTATIONS strings.
    """
    mutation_strings = []
    try:
        with open(file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                if ':' in line:
                    parts = line.split(':', 1)
                    mutations = parts[1].strip()
                    # Add mutation string even if it's empty, as long as identifier was present.
                    # This allows create_input_output_pairs to correctly create (ref,ref) for such cases.
                    mutation_strings.append(mutations)
                    if not mutations:
                        warnings.warn(f"Warning: Empty mutation string (but identifier present) on line {line_num} in {file_path}: '{line}'")
                else:
                    warnings.warn(f"Warning: Line {line_num} in {file_path} does not contain ':' and was skipped: '{line}'")
        return mutation_strings
    except FileNotFoundError:
        print(f"Error: Mutations file not found at {file_path}")
        return []

def apply_mutations(ancestral_sequence_str, mutations_str):
    """
    Applies mutations to an ancestral sequence string.
    Validates mutations against the ancestral sequence.
    """
    if not ancestral_sequence_str:
        warnings.warn("Warning: Ancestral sequence is empty. No mutations applied.")
        return ""
    if not mutations_str: # If mutations_str is empty, return original sequence
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
            warnings.warn(f"Warning: Invalid mutation format '{mut_op}' in '{mutations_str}'. Skipped.")
            continue

        original_base = match.group(1).upper()
        position_1_indexed = int(match.group(2))
        new_base = match.group(3).upper()
        position_0_indexed = position_1_indexed - 1

        if not (0 <= position_0_indexed < len(sequence_list)):
            warnings.warn(f"Warning: Position {position_1_indexed} in mutation '{mut_op}' (from '{mutations_str}') is out of bounds. Skipped.")
            continue

        if sequence_list[position_0_indexed].upper() != original_base:
            warnings.warn(f"Warning: Mismatch for mutation '{mut_op}' (from '{mutations_str}'). Expected '{original_base}' found '{sequence_list[position_0_indexed]}'. Skipped.")
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
        if token is unknown_token: # More robust check if get() returns default
             warnings.warn(f"Warning: Unknown nucleotide '{nucleotide}' found in sequence. Mapped to <UNK>.")
        tokens.append(token)
    return tokens

def create_input_output_pairs(main_reference_sequence_str, all_mutation_strings_list, vocab):
    """
    Creates tokenized input-output pairs.
    Returns a list of (tokenized_ancestral, tokenized_descendant, mutations_str) tuples.
    """
    pairs = []
    if not main_reference_sequence_str:
        warnings.warn("Warning: Main reference sequence is empty for create_input_output_pairs. No pairs created.")
        return pairs

    tokenized_main_reference = tokenize_sequence(main_reference_sequence_str, vocab)

    for i, mutations_str in enumerate(all_mutation_strings_list):
        descendant_sequence_str = apply_mutations(main_reference_sequence_str, mutations_str)
        tokenized_descendant = tokenize_sequence(descendant_sequence_str, vocab)
        pairs.append((tokenized_main_reference, tokenized_descendant, mutations_str))
    return pairs

def pad_sequences(tokenized_pairs, max_len, padding_token=0):
    """
    Pads or truncates tokenized sequences in pairs to max_len.
    Input: list of (tokenized_ancestral, tokenized_descendant, mutations_str)
    Output: list of (padded_ancestral, padded_descendant, mutations_str)
    """
    padded_pairs = []
    for tok_anc, tok_desc, mut_str in tokenized_pairs:
        padded_anc = tok_anc[:max_len] + [padding_token] * max(0, max_len - len(tok_anc))
        padded_desc = tok_desc[:max_len] + [padding_token] * max(0, max_len - len(tok_desc))
        padded_pairs.append((padded_anc, padded_desc, mut_str))
    return padded_pairs

def split_data(all_padded_pairs, train_ratio=0.8, val_ratio=0.1, random_seed=None):
    """
    Splits a list of padded pairs into training, validation, and test sets.
    """
    if not (0 <= train_ratio <= 1 and 0 <= val_ratio <= 1 and (train_ratio + val_ratio) <= 1):
        raise ValueError("Train and validation ratios must be between 0 and 1, and their sum must be <= 1.")

    shuffled_pairs = list(all_padded_pairs) # Make a copy
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

if __name__ == "__main__":
    print("Starting data preparation script for transformer model.")

    reference_fasta_path = "reference.fasta"
    mutations_file_path = "mutations.txt"

    with open(reference_fasta_path, 'w') as f:
        f.write(">test_ref\n")
        f.write("AAAAAGGGGGTTTTTCCCCCNNNNNX\n") # Length 26

    # Create enough entries for a reasonable split. 10 entries for 80/10/10.
    with open(mutations_file_path, 'w') as f:
        f.write("node1: G6T,C18G\n")
        f.write("node2: A1C,G10N,N22A\n")
        f.write("node3: X5Y,G7A\n")
        f.write("node4: C6T\n")
        f.write("node5: N21A,N22T,N23G,N24C,N25A\n")
        f.write("node6: G26X\n")
        f.write("node7: \n")      # Empty mutation string after colon
        f.write("node8: A30T\n") # Out of bounds
        f.write("node9: A2T\n")  # Additional valid mutation
        f.write("node10: C3G\n") # Additional valid mutation


    print(f"\n--- Parsing FASTA file: {reference_fasta_path} ---")
    reference_genome = parse_fasta_file(reference_fasta_path)
    if not reference_genome:
        print("Failed to load reference genome. Exiting.\n")
        exit(1)
    print(f"Loaded reference genome. Length: {len(reference_genome)}")
    print(f"Sequence: {reference_genome}\n")

    max_len = len(reference_genome)
    print(f"Using max_len: {max_len}")

    print(f"\n--- Parsing Mutation file: {mutations_file_path} ---")
    mutation_data_list = parse_mutation_file(mutations_file_path)
    if not mutation_data_list and any(m_str for m_str in mutation_data_list): # Check if list is not empty AND contains non-empty strings
        print("No mutation data loaded or file problem (and not just empty mutation strings). Exiting.\n")
        exit(1)
    print(f"Loaded {len(mutation_data_list)} mutation entries.")

    print(f"\n--- Creating Tokenized Input-Output Pairs ---")
    tokenized_input_output_pairs = create_input_output_pairs(reference_genome, mutation_data_list, NUCLEOTIDE_VOCAB)
    print(f"Created {len(tokenized_input_output_pairs)} tokenized input-output pairs.")

    if not tokenized_input_output_pairs:
        print("No tokenized pairs created. Exiting.\n")
        exit(1)

    # Print details for the first tokenized pair
    # print("\nExample of first tokenized pair (before padding):")
    # tok_anc_ex, tok_desc_ex, mut_str_ex = tokenized_input_output_pairs[0]
    # print(f"  Mutations Applied: {mut_str_ex}")
    # print(f"  Tokenized Ancestral (first 10): {tok_anc_ex[:10]}... Length: {len(tok_anc_ex)}")
    # print(f"  Tokenized Descendant (first 10): {tok_desc_ex[:10]}... Length: {len(tok_desc_ex)}")

    print(f"\n--- Padding Sequences ---")
    padded_pairs = pad_sequences(tokenized_input_output_pairs, max_len, NUCLEOTIDE_VOCAB['<PAD>'])
    print(f"Created {len(padded_pairs)} padded pairs.")

    # if padded_pairs:
    #     print("\nExample of first padded pair:")
    #     pad_anc_ex, pad_desc_ex, mut_str_pad_ex = padded_pairs[0]
    #     print(f"  Mutations Applied: {mut_str_pad_ex}")
    #     print(f"  Padded Ancestral (first 10): {pad_anc_ex[:10]}... Length: {len(pad_anc_ex)}")
    #     print(f"  Padded Descendant (first 10): {pad_desc_ex[:10]}... Length: {len(pad_desc_ex)}")

    #     if len(reference_genome) > 20:
    #         shorter_max_len = 15
    #         print(f"\nTesting padding/truncation with max_len = {shorter_max_len}")
    #         temp_padded_pairs = pad_sequences(tokenized_input_output_pairs[:1], shorter_max_len, NUCLEOTIDE_VOCAB['<PAD>'])
    #         if temp_padded_pairs:
    #             s_pad_anc, s_pad_desc, s_mut_str = temp_padded_pairs[0]
    #             print(f"  Truncated/Padded Ancestral (len {len(s_pad_anc)}): {s_pad_anc}")
    #             print(f"  Truncated/Padded Descendant (len {len(s_pad_desc)}): {s_pad_desc}")

    print(f"\n--- Splitting Data ---")
    if padded_pairs:
        train_data, val_data, test_data = split_data(padded_pairs, train_ratio=0.8, val_ratio=0.1, random_seed=42)
        print(f"Split data into:")
        print(f"  Training set: {len(train_data)} samples")
        print(f"  Validation set: {len(val_data)} samples")
        print(f"  Test set: {len(test_data)} samples")

        if train_data:
            # The third element in the tuple is the mutation string
            print(f"\nExample: First item in training_data (mutations: '{train_data[0][2]}')")
            # print(f"  Padded Ancestral (first 10): {train_data[0][0][:10]}")
            # print(f"  Padded Descendant (first 10): {train_data[0][1][:10]}")
    else:
        print("Skipping data splitting as no padded pairs were available.")

    print("\nData preparation pipeline finished.")
