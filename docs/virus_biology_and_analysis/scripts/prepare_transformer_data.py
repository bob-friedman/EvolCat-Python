"""
Prepare data for viral evolution analysis and Transformer model training.

This script processes a reference FASTA file and a corresponding mutations file
to generate various data formats:
1.  Descendant Sequences: Applies mutations to the reference sequence to generate
    a set of descendant sequences. These can be optionally saved to a multi-FASTA
    file.
2.  Transformer Model Data: Prepares tokenized and padded input-output pairs
    from the ancestral (reference) and descendant sequences. This data is then
    split into training, validation, and test sets suitable for training a
    Transformer model (Note: Transformer model scripts are now in the
    `proof_of_concept` directory).

Key Command-Line Arguments for FASTA Generation:
  --reference_fasta: Path to the input reference FASTA file.
  --mutations_file: Path to the file containing mutation information.
  --mutations_delimiter: Delimiter used in the mutations file (default ':').
  --output_multifasta_file: Path to save the generated descendant sequences
                            in multi-FASTA format. If not specified, this step is skipped.

The script also includes internal functions for parsing, mutation application,
tokenization, padding, and data splitting, which are used to prepare data
for machine learning models.
"""
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

    It reads the specified FASTA file, ignores any header lines (lines starting
    with '>'), and concatenates all other lines to form a single nucleotide
    sequence string. The sequence is converted to uppercase. If multiple
    sequence entries exist in the file, they are all concatenated together.

    Args:
        file_path (str): The path to the FASTA file.

    Returns:
        str: The concatenated nucleotide sequence from the FASTA file.
             Returns an empty string if the file is not found or is empty.
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

    Each line in the mutation file is expected to be in the format:
    `identifier<delimiter>MUTATIONS` (e.g., "node_1:A123G,C456T").
    The function splits each line by the specified `mutations_file_delimiter`.
    The part before the delimiter is taken as the identifier, and the part
    after is taken as the raw mutation string. Both are stripped of leading/trailing
    whitespace.

    Args:
        file_path (str): The path to the mutations file.
        mutations_file_delimiter (str, optional): The delimiter separating
            the identifier from the mutation string in each line. Defaults to ":".

    Returns:
        list: A list of tuples, where each tuple is `(identifier, mutation_string)`.
              If a line does not contain the delimiter, it's skipped with an info message.
              If a line has an identifier but an empty mutation string, the tuple will
              contain an empty string for `mutation_string`.
              Returns an empty list if the file is not found.
    """
    parsed_data_list = []
    try:
        with open(file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line: # Skip empty lines
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
    Applies a comma-separated string of mutations to an ancestral sequence.

    The function takes an ancestral nucleotide sequence string and a string
    containing one or more mutations separated by commas (e.g., "A123G,C456T").
    Each individual mutation operation (e.g., "A123G") is parsed to get the
    original nucleotide, the 1-indexed position, and the new nucleotide.

    Validation checks performed for each mutation:
    -   Format check: Ensures the mutation string matches the pattern (e.g., "A123G").
    -   Out-of-bounds check: Verifies the position is within the sequence length.
    -   Mismatch check: Confirms the original nucleotide in the mutation string
        matches the nucleotide at the specified position in the ancestral sequence.

    If a mutation is invalid, an info message is printed, and that specific
    mutation is skipped. Valid mutations are applied to the sequence.

    Args:
        ancestral_sequence_str (str): The ancestral nucleotide sequence.
        mutations_str (str): A string of comma-separated mutations (e.g., "C1191T,C11674T").
                             Can be empty, in which case the original sequence is returned.

    Returns:
        str: The descendant sequence string after applying valid mutations.
             Returns the original ancestral sequence if `mutations_str` is empty
             or if all mutations within it are invalid.
             Returns an empty string if `ancestral_sequence_str` is empty.
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
    Converts a nucleotide sequence string into a list of integer tokens.

    Each nucleotide in the input `sequence_str` (converted to uppercase) is
    mapped to its corresponding integer token using the provided `vocab`
    dictionary. If a nucleotide is encountered that is not present in the
    `vocab` and `vocab` contains an '<UNK>' key, that nucleotide is mapped
    to the '<UNK>' token's value; an info message is printed for this case.
    If an unknown nucleotide is found and no '<UNK>' token is in vocab,
    it would ideally raise an error (current implementation appends None if '<UNK>' not in vocab).

    Args:
        sequence_str (str): The nucleotide sequence string (e.g., "ACGTN").
        vocab (dict): A dictionary mapping nucleotide characters (and special
                      tokens like '<PAD>', '<UNK>') to integer IDs.

    Returns:
        list: A list of integer tokens representing the input sequence.
    """
    tokens = []
    unknown_token = vocab.get('<UNK>') # Get the token for unknown characters

    for nucleotide in sequence_str.upper(): # Ensure uppercase for consistency with vocab
        token = vocab.get(nucleotide, unknown_token) # Default to unknown_token if not found
        if token is unknown_token and nucleotide not in vocab: # Print info only if truly unknown and not '<UNK>' itself
             print(f"Info: Unknown nucleotide '{nucleotide}' found in sequence. Mapped to <UNK>.")
        tokens.append(token)
    return tokens

def create_input_output_pairs(main_reference_sequence_str, identifier_mutation_list, vocab):
    """
    Generates tokenized (ancestral, descendant) pairs for Transformer model training.

    For each entry in `identifier_mutation_list` (which contains tuples of
    (identifier, mutation_string)), this function:
    1.  Applies the `mutation_string` to the `main_reference_sequence_str` to get a
        `descendant_sequence_str`.
    2.  Tokenizes both the `main_reference_sequence_str` (ancestral) and the
        `descendant_sequence_str` using `tokenize_sequence` and the provided `vocab`.
    3.  Stores the result as a tuple:
        `(tokenized_ancestral, tokenized_descendant, identifier, mutation_string)`.
        The identifier and original mutation string are carried along for context.

    Args:
        main_reference_sequence_str (str): The single ancestral sequence string.
        identifier_mutation_list (list): A list of tuples, where each tuple is
                                         `(identifier, mutation_string)`.
        vocab (dict): The nucleotide-to-token vocabulary.

    Returns:
        list: A list of tuples, each being
              `(tokenized_ancestral, tokenized_descendant, identifier, mutation_string)`.
              Returns an empty list if `main_reference_sequence_str` is empty.
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
    Pads or truncates tokenized ancestral and descendant sequences in pairs to a fixed length.

    Each element in `tokenized_pairs` is expected to be a tuple:
    `(tokenized_ancestral, tokenized_descendant, identifier, mutations_str)`.
    Both `tokenized_ancestral` and `tokenized_descendant` sequences are processed:
    - If a sequence is longer than `max_len`, it is truncated to `max_len`.
    - If a sequence is shorter than `max_len`, it is padded with `padding_token`
      (default 0, corresponding to '<PAD>') until it reaches `max_len`.

    Args:
        tokenized_pairs (list): A list of tuples, each containing tokenized
                                ancestral and descendant sequences, an identifier,
                                and a mutation string.
        max_len (int): The target length for all sequences after padding/truncation.
        padding_token (int, optional): The token ID used for padding. Defaults to 0.

    Returns:
        list: A list of tuples, where tokenized sequences are replaced by their
              padded/truncated versions:
              `(padded_ancestral, padded_descendant, identifier, mutations_str)`.
    """
    padded_pairs = []
    for tok_anc, tok_desc, identifier, mut_str in tokenized_pairs:
        padded_anc = tok_anc[:max_len] + [padding_token] * max(0, max_len - len(tok_anc))
        padded_desc = tok_desc[:max_len] + [padding_token] * max(0, max_len - len(tok_desc))
        padded_pairs.append((padded_anc, padded_desc, identifier, mut_str))
    return padded_pairs

def split_data(all_padded_pairs, train_ratio=0.8, val_ratio=0.1, random_seed=None):
    """
    Splits a list of data items (e.g., padded pairs) into training, validation, and test sets.

    The data is first shuffled (optionally with a `random_seed` for reproducibility).
    Then, it's split based on the provided `train_ratio` and `val_ratio`.
    The test set size is implicitly `1.0 - train_ratio - val_ratio`.

    Args:
        all_padded_pairs (list): The list of data items to split.
        train_ratio (float, optional): Proportion of data for the training set. Defaults to 0.8.
        val_ratio (float, optional): Proportion of data for the validation set. Defaults to 0.1.
        random_seed (int, optional): Seed for `random.shuffle` for reproducible splits. Defaults to None.

    Returns:
        tuple: A tuple containing three lists: `(train_set, val_set, test_set)`.

    Raises:
        ValueError: If ratios are invalid (e.g., not between 0 and 1, or sum > 1).
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

    Each sequence is written with its header, followed by the sequence itself
    on a new line. This implementation writes the entire sequence on a single line
    without wrapping, which is simpler but might not be ideal for extremely long
    sequences viewed in some FASTA tools.

    Args:
        filepath (str): The path to the output multi-FASTA file.
        sequences_data (list): A list of tuples, where each tuple is
                               `(header_string, sequence_string)`.

    Handles `IOError` if the file cannot be written and prints an error message.
    """
    try:
        with open(filepath, 'w') as f:
            for header, sequence in sequences_data:
                f.write(f">{header}\n")
                f.write(f"{sequence}\n")
        print(f"Info: Successfully wrote {len(sequences_data)} sequences to {filepath}")
    except IOError as e:
        print(f"Error: Could not write to file {filepath}. Reason: {e}")

# This function is not used by the main script but is kept for potential direct use or testing.
def to_tf_dataset(pairs_list, batch_size, max_len_for_padding): # Added max_len for empty dataset shape
    """
    Converts a list of processed pairs into a TensorFlow Dataset.
    This version is designed for data structured as (ancestral, descendant, identifier, mutation_str).
    It prepares data for a sequence-to-sequence model where:
    - Encoder input: padded_tokenized_ancestral
    - Decoder input (teacher forcing): padded_tokenized_descendant (shifted right, e.g. [START] + seq[:-1])
    - Target for loss: padded_tokenized_descendant (shifted left, e.g. seq[1:] + [END])

    Note: The current implementation uses p[1][:-1] for decoder input and p[1][1:] for target,
    implying a simple shift without explicit START/END tokens for this version.

    Args:
        pairs_list (list): List of tuples, each being (padded_ancestral, padded_descendant, identifier, mut_str).
        batch_size (int): The batch size for the dataset.
        max_len_for_padding (int): The sequence length used for padding, needed for empty dataset signature.


    Returns:
        tf.data.Dataset: A TensorFlow Dataset object, batched and prefetched.
    """
    if not pairs_list:
        # Create an empty dataset with the expected structure and shapes
        # This helps prevent errors if, for example, val_pairs or test_pairs is empty.
        # Shapes are (batch_size, sequence_length). Decoder input and target are one shorter due to shifting.
        # The sequence length for inputs_dec and targets should be max_len_for_padding - 1.
        return tf.data.Dataset.from_tensor_slices(
            ((tf.zeros((0, max_len_for_padding), dtype=tf.int32), tf.zeros((0, max_len_for_padding - 1), dtype=tf.int32) ),
             tf.zeros((0, max_len_for_padding - 1), dtype=tf.int32)) # Target shape
        ).batch(batch_size).prefetch(tf.data.AUTOTUNE)

    # Extract sequences for model input/output
    # p[0] is padded_ancestral, p[1] is padded_descendant
    inputs_enc = np.array([p[0] for p in pairs_list], dtype=np.int32)
    inputs_dec = np.array([p[1][:-1] for p in pairs_list], dtype=np.int32) # Decoder input: seq[:-1]
    targets = np.array([p[1][1:] for p in pairs_list], dtype=np.int32)    # Target output: seq[1:]

    dataset = tf.data.Dataset.from_tensor_slices(((inputs_enc, inputs_dec), targets))
    dataset = dataset.shuffle(len(pairs_list), reshuffle_each_iteration=True) # Shuffle before batching
    dataset = dataset.batch(batch_size)
    dataset = dataset.prefetch(tf.data.AUTOTUNE) # Optimize performance
    return dataset


def main():
    # --- Argument Parsing ---
    # Sets up command-line arguments for configuring the script's behavior.
    parser = argparse.ArgumentParser(description="Prepare data for Transformer model and optionally generate FASTA files.")
    parser.add_argument('--reference_fasta', type=str, default="dummy_reference.fasta", help='Path to the reference FASTA file.')
    parser.add_argument('--mutations_file', type=str, default="dummy_mutations.txt", help='Path to the mutations file.')
    parser.add_argument('--mutations_delimiter', type=str, default=":", help='Delimiter used in the mutations file.')
    parser.add_argument('--output_multifasta_file', type=str, default="output_fasta_sequences/descendant_sequences.fasta", help='Path to save the output multi-FASTA file. If empty or not provided, FASTA generation is skipped.')

    args = parser.parse_args()

    # Use parsed arguments for file paths and parameters
    REFERENCE_FASTA_PATH = args.reference_fasta
    MUTATIONS_FILE_PATH = args.mutations_file
    mutations_file_delimiter = args.mutations_delimiter
    output_fasta_filepath = args.output_multifasta_file

    print("Starting data preparation script with arguments:")
    print(f"  Reference FASTA: {REFERENCE_FASTA_PATH}")
    print(f"  Mutations File: {MUTATIONS_FILE_PATH}")
    print(f"  Mutations Delimiter: '{mutations_file_delimiter}'")
    print(f"  Output Multi-FASTA: {output_fasta_filepath if output_fasta_filepath else 'Not generating'}")

    # --- Dummy File Creation (if defaults are used and files don't exist) ---
    # This allows the script to run with example data for testing purposes.
    if REFERENCE_FASTA_PATH == "dummy_reference.fasta" and not os.path.exists(REFERENCE_FASTA_PATH):
        with open(REFERENCE_FASTA_PATH, 'w') as f:
            f.write(">dummy_ref_genome\nAAAAAGGGGGTTTTTCCCCCNNNNNX") # Length 26
        print(f"Info: Created dummy reference FASTA: {REFERENCE_FASTA_PATH}")

    if MUTATIONS_FILE_PATH == "dummy_mutations.txt" and not os.path.exists(MUTATIONS_FILE_PATH):
        with open(MUTATIONS_FILE_PATH, 'w') as f:
            # Example mutations for dummy data
            f.write("node1:G6T,C18G\nnode2:A1C,G10N,N22A\nnode3:X5Y,G7A\nnode4:C6T\nnode5:N21A,N22T,N23G,N24C,N25A\nnode6:G26X\nnode7:\nnode8:A30T\nnode9:A2T\nnode10:C3G")
        print(f"Info: Created dummy mutations file: {MUTATIONS_FILE_PATH}")

    # --- Output Directory Creation for FASTA ---
    # Ensures the directory for the output FASTA file exists before writing.
    if output_fasta_filepath:
        output_fasta_dir = os.path.dirname(output_fasta_filepath)
        if output_fasta_dir and not os.path.exists(output_fasta_dir): # Check if output_fasta_dir is not empty
            os.makedirs(output_fasta_dir, exist_ok=True)
            print(f"Info: Created output directory for FASTA: {output_fasta_dir}")

    # --- Core Data Processing ---
    print(f"\n--- Parsing FASTA file: {REFERENCE_FASTA_PATH} ---")
    reference_genome = parse_fasta_file(REFERENCE_FASTA_PATH)
    if not reference_genome: # Critical error if reference cannot be loaded
        print(f"Error: Could not load reference genome from {REFERENCE_FASTA_PATH}. Exiting.")
        return

    print(f"Loaded reference genome. Length: {len(reference_genome)}")
    print(f"Sequence (first 60 chars): {reference_genome[:60]}\n")

    # max_len is determined by the actual reference genome length for this script's run
    max_len = len(reference_genome)
    print(f"Using max_len for padding/analysis: {max_len}")

    print(f"\n--- Parsing Mutation file: {MUTATIONS_FILE_PATH} ---")
    identifier_mutation_list = parse_mutation_file(MUTATIONS_FILE_PATH, mutations_file_delimiter=mutations_file_delimiter)

    if not identifier_mutation_list:
        print("Info: No mutation data loaded or file problem. FASTA generation will only contain reference if specified. Transformer data pipeline might be empty.")
    else:
        print(f"Loaded {len(identifier_mutation_list)} mutation entries (identifier, mutations_str).")
        if identifier_mutation_list: # Print example only if list is not empty
            print(f"  Example entry: Identifier='{identifier_mutation_list[0][0]}', Mutations='{identifier_mutation_list[0][1]}'")

    # --- FASTA File Generation ---
    # This section generates descendant sequences from the reference and mutations,
    # then writes them to a multi-FASTA file if an output path is specified.
    if output_fasta_filepath:
        descendant_sequences_for_fasta = []
        if reference_genome:
            print(f"\n--- Generating descendant sequences for FASTA output... ---")
            # Always include the reference sequence itself in the FASTA output.
            ref_fasta_header = os.path.basename(REFERENCE_FASTA_PATH).replace(".fasta", "").replace(".fa", "") + "_reference"
            descendant_sequences_for_fasta.append((ref_fasta_header, reference_genome))

            # Generate descendant sequences from mutations
            if identifier_mutation_list:
                for identifier, mutation_str in identifier_mutation_list:
                    if not mutation_str: # If mutation string is empty, descendant is same as reference
                        descendant_seq = reference_genome
                        print(f"Info: No mutations for identifier '{identifier}', using reference sequence for its FASTA entry.")
                    else:
                        descendant_seq = apply_mutations(reference_genome, mutation_str)
                    descendant_sequences_for_fasta.append((identifier, descendant_seq))

            # Write all collected sequences (reference + descendants) to the multi-FASTA file.
            if descendant_sequences_for_fasta:
                write_multifasta(output_fasta_filepath, descendant_sequences_for_fasta)
            # This else case should ideally not be reached if reference_genome is present.
            # else:
            #     print(f"Info: No descendant sequences (including reference) to write. FASTA file '{output_fasta_filepath}' will not be created or will be empty.")
        else:
            # This case should not be reached due to earlier exit if reference_genome is not loaded.
            print("Info: Reference genome not loaded. Cannot generate descendant sequences for FASTA.")


    # --- Data Preparation for Transformer Model ---
    # This section prepares data specifically for training a Transformer model.
    # It involves tokenization, creating input-output pairs, padding, and splitting.
    print(f"\n--- Creating Tokenized Input-Output Pairs for Transformer ---")

    # Use identifier_mutation_list for Transformer data; if it's empty, create a self-pair.
    # This ensures the Transformer pipeline can run even if no external mutations are provided.
    identifier_mutation_list_for_transformer = identifier_mutation_list
    if not identifier_mutation_list and reference_genome:
        print("Info: No mutations found from file, creating a self-reference pair for transformer data.")
        identifier_mutation_list_for_transformer = [("reference_self_pair", "")]
    elif not reference_genome: # Should have exited earlier if no reference.
        print("Error: Reference genome not available. Cannot create pairs for transformer.")
        return # Stop if no reference genome for Transformer data prep.

    tokenized_input_output_pairs = create_input_output_pairs(reference_genome, identifier_mutation_list_for_transformer, NUCLEOTIDE_VOCAB)

    if not tokenized_input_output_pairs: # If list is empty (e.g. reference was empty)
        print("No tokenized pairs created for transformer. Ending Transformer data preparation.")
        return

    print(f"Created {len(tokenized_input_output_pairs)} tokenized input-output pairs for Transformer.")
    if tokenized_input_output_pairs: # Print example if pairs exist
         print(f"  Example pair content (first pair): Ancestral tokens ({len(tokenized_input_output_pairs[0][0])}), Descendant tokens ({len(tokenized_input_output_pairs[0][1])}), Identifier ('{tokenized_input_output_pairs[0][2]}'), Mutations ('{tokenized_input_output_pairs[0][3]}')")

    print(f"\n--- Padding Sequences for Transformer ---")
    padded_pairs = pad_sequences(tokenized_input_output_pairs, max_len, NUCLEOTIDE_VOCAB['<PAD>'])
    print(f"Created {len(padded_pairs)} padded pairs for Transformer.")

    print(f"\n--- Splitting Data for Transformer ---")
    if padded_pairs:
        train_data, val_data, test_data = split_data(padded_pairs, train_ratio=0.8, val_ratio=0.1, random_seed=42)
        print(f"Split data into (for Transformer):")
        print(f"  Training set: {len(train_data)} samples")
        print(f"  Validation set: {len(val_data)} samples")
        print(f"  Test set: {len(test_data)} samples")

        if train_data: # Print example if training data exists
            print(f"\nExample: First item in training_data (identifier: '{train_data[0][2]}', mutations: '{train_data[0][3]}')")
    else:
        print("Skipping data splitting as no padded pairs were available for Transformer.")

    print("\nData preparation script finished.")

if __name__ == "__main__":
    main()
```
