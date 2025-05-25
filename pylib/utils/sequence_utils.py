import re

def generate_kmers(sequence_string, kmer_len, step=1):
    """
    Generates k-mers from a sequence string.

    Args:
        sequence_string (str): The input sequence.
        kmer_len (int): The length of the k-mers to generate.
        step (int): The step size to move along the sequence. Default is 1 (overlapping k-mers).

    Yields:
        str: The k-mers found in the sequence.
    """
    if not isinstance(sequence_string, str):
        raise TypeError("Input sequence must be a string.")
    if not isinstance(kmer_len, int) or kmer_len <= 0:
        raise ValueError("K-mer length must be a positive integer.")
    if not isinstance(step, int) or step <= 0:
        raise ValueError("Step must be a positive integer.")

    if len(sequence_string) < kmer_len:
        return # Or raise ValueError("Sequence length is less than k-mer length.")

    for i in range(0, len(sequence_string) - kmer_len + 1, step):
        yield sequence_string[i : i + kmer_len]

def is_valid_dna(sequence_string, alphabet={'A', 'T', 'C', 'G'}):
    """
    Validates if a sequence string contains only characters from the given DNA alphabet.
    Assumes sequence_string is already uppercase.

    Args:
        sequence_string (str): The input sequence (expected to be uppercase).
        alphabet (set): A set of allowed characters. Default is {'A', 'T', 'C', 'G'}.

    Returns:
        bool: True if the sequence is valid DNA, False otherwise.
    """
    if not isinstance(sequence_string, str):
        return False # Or raise TypeError
    
    # Using a set for efficient lookup of unique characters in the sequence
    # against the provided alphabet.
    return set(sequence_string).issubset(alphabet)

def is_valid_dna_strict_regex(sequence_string):
    """
    Validates if a sequence string contains only A, T, C, G using regex.
    Assumes sequence_string is already uppercase.
    This is similar to what was implemented in count_kmers.py initially.
    """
    if not isinstance(sequence_string, str):
        return False
    return bool(re.fullmatch(r"^[ATCG]+$", sequence_string))
