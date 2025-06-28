"""
Evolutionary Data Retrieval and Analysis Pipeline

This script automates several steps in a common bioinformatics workflow:
1.  **BLAST Search:** Takes a protein sequence (e.g., SARS-CoV-2 RBD) and searches for
    homologous sequences in the NCBI 'nr' database using BLASTp.
2.  **Parse BLAST Results:** Extracts homologous sequences from the BLAST XML output.
3.  **Multiple Sequence Alignment (MSA):** Aligns the query sequence and its homologs
    using MAFFT to create an MSA.
4.  **Conservation Analysis:** Calculates per-position conservation scores based on the MSA.
5.  **Structural Analysis (Burial Score):**
    -   Downloads a specified PDB file.
    -   Runs DSSP (Define Secondary Structure of Proteins) to calculate Relative Solvent
        Accessibility (RSA) for each residue.
    -   Calculates a "Burial Score" (1 - RSA) and adds it to the conservation data.
6.  **Physicochemical Analysis:**
    -   Calculates a Physicochemical Score (P-score) for potential amino acid substitutions
        based on the BLOSUM62 matrix. This score indicates how radical a mutation is.

Workflow:
- Define Query Sequence: A specific protein sequence (e.g., SARS-CoV-2 RBD) is defined.
- Run BLAST: `run_blast_search` function sends the query to NCBI BLAST.
- Parse Results: `parse_blast_xml` function reads the XML output and extracts sequences.
- Generate MSA: `generate_msa_with_mafft` function uses MAFFT to align sequences.
- Calculate Conservation: `calculate_conservation_scores` function analyzes the MSA.
- Calculate Burial Score:
    - `fetch_pdb_file` (helper, not explicitly in original script but implied for PDB)
    - `inspect_pdb_for_chain_residues` function helps identify correct residue numbering.
    - `calculate_burial_scores` function uses DSSP and PDB data.
- Calculate Physicochemical Score: `calculate_p_score` function uses BLOSUM62.
- Combine Data: Results from conservation and burial scores are combined into a DataFrame.

Inputs:
-   Target protein sequence string (hardcoded or configurable).
-   PDB ID (e.g., "7VYR") for structural analysis.
-   Chain ID within the PDB file (e.g., 'C').
-   Residue number range for DSSP analysis (e.g., 336-521).

Outputs:
-   `blast_results.xml`: Raw XML output from NCBI BLAST.
-   `unaligned_sequences.fasta`: Homologous sequences before alignment.
-   `aligned_sequences.fasta`: Multiple Sequence Alignment in FASTA format.
-   Console output: Progress messages, summaries of results, snippets of data.
-   A pandas DataFrame (`filter_df`) containing position, reference AA, conservation score,
    RSA, and burial score.
-   Example P-scores for sample mutations.

Dependencies:
-   BioPython (`biopython`)
-   pandas
-   MAFFT (must be installed separately, e.g., `sudo apt-get install mafft`)
-   DSSP (mkdssp) (must be installed separately, e.g., `sudo apt-get install dssp`)

Notes:
-   The script assumes MAFFT and DSSP (mkdssp) are in the system PATH.
-   NCBI BLAST queries can be time-consuming.
-   File paths for outputs are relative to the script's execution directory.
-   Error handling is included for common issues like file not found or API errors.
-   The final part of the original script ("loop through all 187 positions...") for P-scores
    is described but not fully implemented as a loop saving results; this refactoring
    maintains the P-score calculation function and its test.
"""

import os
import argparse
from io import StringIO
from collections import Counter

# BioPython modules
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MafftCommandline
from Bio.PDB import PDBParser, PDBList
from Bio.PDB.DSSP import DSSP
from Bio.Align import substitution_matrices

# Other libraries
import pandas as pd

# --- Configuration & Constants ---
# Default sequence: SARS-CoV-2 Spike RBD (Wuhan-Hu-1, GenBank: MN908947.3, residues 319-541 of Spike)
# The provided sequence is 201 AA. RBD is often cited as ~331-531.
# Let's use the provided sequence and slicing.
DEFAULT_RBD_SEQUENCE_FULL = "NIDTFLDSVYKGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVN"
DEFAULT_RBD_SUBSEQUENCE_SLICE = (330, 531) # Python slice indices (0-indexed start, exclusive end) for 331-531

DEFAULT_PDB_ID = "7VYR" # Example PDB ID for Spike RBD complex
DEFAULT_PDB_CHAIN_ID = 'C' # Example chain
DEFAULT_DSSP_START_RES = 336 # Example start residue for DSSP analysis on PDB 7VYR chain C
DEFAULT_DSSP_END_RES = 521   # Example end residue

OUTPUT_BLAST_XML_FILE = "blast_results.xml"
OUTPUT_UNALIGNED_FASTA_FILE = "unaligned_sequences.fasta"
OUTPUT_ALIGNED_FASTA_FILE = "aligned_sequences.fasta"
OUTPUT_PDB_FILE_PREFIX = DEFAULT_PDB_ID.lower() # PDB files are often lowercase, e.g., 7vyr.pdb

# --- Helper Functions ---
def get_rbd_subsequence(full_sequence_str: str, slice_tuple: tuple) -> str:
    """Extracts a subsequence using Python slicing."""
    return full_sequence_str[slice_tuple[0]:slice_tuple[1]]

def create_seq_record(sequence_str: str, record_id: str, description: str) -> SeqRecord:
    """Creates a Bio.SeqRecord object."""
    return SeqRecord(Seq(sequence_str), id=record_id, description=description)

# --- Core Workflow Functions ---

def run_blast_search(query_record: SeqRecord, output_xml_file: str = OUTPUT_BLAST_XML_FILE, blast_program: str = "blastp", database: str = "nr") -> bool:
    """
    Performs a BLAST search against the NCBI database and saves the results.
    Returns True on success, False on failure.
    """
    print(f"\n--- Running BLAST Search ({blast_program} against {database}) ---")
    print(f"Query: {query_record.id} ({len(query_record.seq)} aa)")
    print("This may take several minutes depending on NCBI server load. Please be patient.")
    try:
        result_handle = NCBIWWW.qblast(blast_program, database, query_record.seq)
        with open(output_xml_file, "w") as out_handle:
            out_handle.write(result_handle.read())
        print(f"BLAST search complete. Results saved to '{output_xml_file}'.")
        return True
    except Exception as e:
        print(f"An error occurred during the BLAST search: {e}")
        print("This can be due to network issues or NCBI server problems. Consider trying again later.")
        return False

def parse_blast_xml(blast_xml_file: str, query_record: SeqRecord) -> list[SeqRecord]:
    """
    Parses a BLAST XML file and extracts homologous sequences as SeqRecord objects.
    The original query sequence is added first to the list.
    """
    print(f"\n--- Parsing BLAST XML: {blast_xml_file} ---")
    homologous_sequences = [query_record] # Add query sequence first

    if not os.path.exists(blast_xml_file):
        print(f"ERROR: BLAST XML file '{blast_xml_file}' not found.")
        return homologous_sequences # Return list with only query

    try:
        with open(blast_xml_file, "r") as result_handle:
            blast_records = NCBIXML.parse(result_handle) # Iterator
            for blast_record in blast_records: # Should be only one for a single query
                print(f"Processing hits for query: {blast_record.query_id} ({blast_record.query_letters} letters)")
                print(f"Found {len(blast_record.alignments)} alignments (homologous sequences).")
                for alignment in blast_record.alignments:
                    # Each alignment can have multiple HSPs. Typically, the first HSP is the most significant.
                    if not alignment.hsps:
                        continue
                    hsp = alignment.hsps[0] # Use the first HSP

                    hit_id = alignment.hit_id
                    # Clean up description: often looks like "pdb|6M0J|A Chain A, Spike..."
                    hit_def = alignment.hit_def.split(']')[0] + ']' if ']' in alignment.hit_def else alignment.hit_def
                    hit_sequence_str = hsp.sbjct # The aligned part of the subject sequence (hit)

                    hit_record = create_seq_record(hit_sequence_str, hit_id, hit_def)
                    homologous_sequences.append(hit_record)
        print(f"Extracted {len(homologous_sequences) - 1} homologous sequences from BLAST results.")
        print("First 5 extracted sequence IDs (includes query):")
        for seq_rec in homologous_sequences[:5]:
            print(f"- {seq_rec.id}: {seq_rec.description[:70]}...")
    except Exception as e:
        print(f"An error occurred during BLAST XML parsing: {e}")
    return homologous_sequences


def generate_msa_with_mafft(sequences: list[SeqRecord],
                              unaligned_file: str = OUTPUT_UNALIGNED_FASTA_FILE,
                              aligned_file: str = OUTPUT_ALIGNED_FASTA_FILE) -> AlignIO.MultipleSeqAlignment | None:
    """
    Generates a Multiple Sequence Alignment (MSA) using MAFFT.
    Requires MAFFT to be installed and in the system PATH.
    Saves unaligned and aligned sequences to FASTA files.
    Returns the MSA object or None on failure.
    """
    print("\n--- Generating Multiple Sequence Alignment (MSA) with MAFFT ---")
    if not sequences:
        print("No sequences provided for MSA.")
        return None

    SeqIO.write(sequences, unaligned_file, "fasta")
    print(f"Saved {len(sequences)} unaligned sequences to '{unaligned_file}'.")
    print("Running MAFFT for alignment. This may take some time for many/long sequences...")

    try:
        # Check if mafft is available
        mafft_cline_test = MafftCommandline(clwstrict=True) # A simple way to test

        mafft_cline = MafftCommandline(input=unaligned_file)
        # Common MAFFT options:
        # --auto: Automatically selects appropriate strategy
        # --thread N: Use N threads
        # mafft_cline.auto = True
        # mafft_cline.thread = 4 # Example: use 4 threads if available

        stdout, stderr = mafft_cline() # Execute MAFFT

        if stderr:
            print(f"MAFFT STDERR:\n{stderr}")

        msa = AlignIO.read(StringIO(stdout), "fasta")
        AlignIO.write(msa, aligned_file, "fasta")

        print("MSA created successfully.")
        print(f"The MSA has {len(msa)} sequences and an alignment length of {msa.get_alignment_length()} columns.")
        print(f"Alignment saved to '{aligned_file}'.")
        print("\nSnippet of the final MSA (first 60 columns):")
        print(msa[:, :60])
        return msa
    except FileNotFoundError:
        print("ERROR: MAFFT command not found. Please ensure MAFFT is installed and in your system PATH.")
        return None
    except Exception as e:
        print(f"An error occurred during MAFFT alignment: {e}")
        return None

def calculate_conservation_scores(msa: AlignIO.MultipleSeqAlignment) -> pd.DataFrame:
    """
    Calculates per-position conservation scores from an MSA.
    Conservation is defined as the frequency of the most common amino acid at that position.
    The reference sequence for 'Reference_AA' is taken as the first sequence in the MSA.
    Returns a pandas DataFrame with scores.
    """
    print("\n--- Calculating Conservation Scores ---")
    if msa is None or len(msa) == 0:
        print("MSA object is empty or not provided. Cannot calculate conservation.")
        return pd.DataFrame()

    conservation_data = []
    num_sequences = len(msa)
    alignment_length = msa.get_alignment_length()
    reference_sequence = msa[0].seq # Query sequence is assumed to be the first

    for i in range(alignment_length): # Iterate through each column of the MSA
        column_residues = msa[:, i] # All residues in the current column
        counts = Counter(column_residues)

        most_common_count = 0
        if counts:
            # Exclude gaps from being the "most common" if desired, though simple frequency is used here.
            # most_common_aa, most_common_count = counts.most_common(1)[0]

            # If you want to ignore gaps in the count for the most common residue:
            # non_gap_counts = Counter(res for res in column_residues if res != '-')
            # if non_gap_counts:
            #     most_common_count = non_gap_counts.most_common(1)[0][1]
            # else: # Column is all gaps
            #     most_common_count = 0
            most_common_count = counts.most_common(1)[0][1]


        score = most_common_count / num_sequences if num_sequences > 0 else 0
        reference_aa = reference_sequence[i]

        conservation_data.append({
            "Position": i + 1,  # 1-based numbering for biological convention
            "Reference_AA": reference_aa,
            "Conservation_Score": score
        })

    conservation_df = pd.DataFrame(conservation_data)
    print("Conservation scores calculated.")
    print("\nFirst 5 positions (Conservation):")
    print(conservation_df.head(5))
    # print("\nScores for specific positions (e.g., N501Y location if applicable to your query length):")
    # Example: If N501Y is at position 171 in your ~187aa query (501 - 330 = 171)
    # if 171 <= len(conservation_df):
    #     print(conservation_df[conservation_df['Position'] == 171])
    return conservation_df

def fetch_pdb_file(pdb_id: str, download_dir: str = ".") -> str | None:
    """
    Downloads a PDB file using Bio.PDB.PDBList.
    Returns the path to the downloaded file or None on failure.
    """
    print(f"\n--- Fetching PDB File: {pdb_id} ---")
    pdbl = PDBList()
    # Structure of PDBList().retrieve_pdb_file:
    # It downloads to a directory like <download_dir>/pdb_id[1:3]/pdb<pdb_id>.ent
    # We want it as <pdb_id>.pdb or <pdb_id>.cif in the current directory for simplicity.
    # It returns the path to the downloaded file.
    try:
        # This will download to a subdirectory named after middle chars of PDB ID, e.g. "vy" for 7vyr
        # And filename like pdb7vyr.ent
        # For DSSP, we often need .pdb or .cif format. PDBList usually gets .ent (PDB format) or .cif
        filepath = pdbl.retrieve_pdb_file(pdb_id, pdir=download_dir, file_format='pdb') # 'pdb', 'cif', 'xml', 'mmtf'

        if not os.path.exists(filepath): # Should not happen if retrieve_pdb_file succeeded
            print(f"Error: PDB file for {pdb_id} downloaded by PDBList but not found at expected path: {filepath}")
            # Try to find it if naming is different
            expected_path_pdb = os.path.join(download_dir, f"{pdb_id.lower()}.pdb")
            expected_path_ent = os.path.join(download_dir, f"pdb{pdb_id.lower()}.ent")
            if os.path.exists(expected_path_pdb): filepath = expected_path_pdb
            elif os.path.exists(expected_path_ent): filepath = expected_path_ent
            else: return None

        # Rename to <pdb_id>.pdb for consistency if needed
        base_dir = os.path.dirname(filepath) # Might be a subdirectory
        final_pdb_path = os.path.join(download_dir, f"{pdb_id.lower()}.pdb")

        if filepath != final_pdb_path :
            if os.path.exists(final_pdb_path): # If it already exists from a previous run
                 os.remove(final_pdb_path) # Remove to avoid error on rename
            os.rename(filepath, final_pdb_path)
            # If filepath was in a subdirectory created by PDBList, remove that dir if empty
            if base_dir != download_dir and os.path.exists(base_dir) and not os.listdir(base_dir):
                os.rmdir(base_dir)
            print(f"PDB file saved as '{final_pdb_path}'.")
            return final_pdb_path
        else:
            print(f"PDB file already at '{final_pdb_path}'.")
            return filepath

    except Exception as e:
        print(f"Error downloading PDB file {pdb_id}: {e}")
        return None

def inspect_pdb_for_chain_residues(pdb_filepath: str, chain_id: str, dssp_path: str = 'mkdssp'):
    """
    Inspects a PDB file using DSSP to find actual residue numbering for a given chain.
    This is a diagnostic step.
    """
    print(f"\n--- Inspecting PDB File for Chain '{chain_id}' Residue Numbering ---")
    if not os.path.exists(pdb_filepath):
        print(f"ERROR: PDB file '{pdb_filepath}' not found for inspection.")
        return

    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("PDB_structure", pdb_filepath)
        model = structure[0] # Assuming single model PDB

        # Check if DSSP can run before proceeding
        try:
            dssp = DSSP(model, pdb_filepath, dssp=dssp_path)
        except FileNotFoundError:
            print(f"ERROR: DSSP executable ('{dssp_path}') not found. Please install DSSP and ensure it's in PATH.")
            print("Skipping PDB inspection and subsequent burial score calculation.")
            return
        except Exception as e: # Other DSSP errors
             print(f"An error occurred while trying to run DSSP: {e}")
             print("Skipping PDB inspection and subsequent burial score calculation.")
             return


        present_keys = list(dssp.keys()) # List of (chain_id, residue_id_tuple)
        chain_residue_keys = [key for key in present_keys if key[0] == chain_id]

        if not chain_residue_keys:
            print(f"No residues found for Chain '{chain_id}' in PDB file '{pdb_filepath}' via DSSP.")
            print("Available chain IDs in DSSP output:", sorted(list(set(key[0] for key in present_keys))))
        else:
            # residue_id_tuple is like (' ', res_num, ' ') for standard residues
            chain_residue_numbers = sorted([key[1][1] for key in chain_residue_keys if key[1][0] == ' ']) # Filter for standard residues
            if chain_residue_numbers:
                print(f"Found {len(chain_residue_numbers)} residues for Chain '{chain_id}'.")
                print(f"Residue numbering in PDB for Chain '{chain_id}' starts at: {chain_residue_numbers[0]}")
                print(f"Residue numbering in PDB for Chain '{chain_id}' ends at:   {chain_residue_numbers[-1]}")
                # print("First 10 residue numbers found:", chain_residue_numbers[:10])
            else:
                print(f"No standard residues (hetflag=' ') found for Chain '{chain_id}'.")
    except Exception as e:
        print(f"An error occurred during PDB inspection: {e}")


def calculate_burial_scores(conservation_df: pd.DataFrame, pdb_id: str,
                             pdb_filepath: str, chain_id: str,
                             dssp_start_res: int, dssp_end_res: int,
                             dssp_path: str = 'mkdssp') -> pd.DataFrame:
    """
    Calculates Relative Solvent Accessibility (RSA) using DSSP and PDB data,
    then computes a Burial Score (1 - RSA). Adds these to the conservation DataFrame.
    Requires DSSP (mkdssp) to be installed.
    The length of RSA values extracted must match the length of the conservation_df
    (which is based on MSA length). Adjustments (padding/truncating) might be needed if they differ.
    """
    print("\n--- Calculating Burial Scores (RSA from DSSP) ---")
    if conservation_df.empty:
        print("Conservation DataFrame is empty. Cannot add burial scores.")
        return pd.DataFrame()
    if not os.path.exists(pdb_filepath):
        print(f"ERROR: PDB file '{pdb_filepath}' not found for DSSP.")
        return conservation_df # Return original df

    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_id, pdb_filepath)
        model = structure[0]

        try:
            dssp = DSSP(model, pdb_filepath, dssp=dssp_path)
        except FileNotFoundError:
            print(f"ERROR: DSSP executable ('{dssp_path}') not found. Please install DSSP and ensure it's in PATH.")
            return conservation_df
        except Exception as e: # Other DSSP errors
             print(f"An error occurred while trying to run DSSP: {e}")
             return conservation_df


        rsa_values = []
        print(f"Extracting RSA for Chain '{chain_id}', PDB residues {dssp_start_res} to {dssp_end_res}.")
        for res_num_pdb in range(dssp_start_res, dssp_end_res + 1):
            dssp_key = (chain_id, (' ', res_num_pdb, ' ')) # (' ', res_num, insertion_code)
            try:
                # DSSP output: (dssp_index, aa, sec_struct, rsa, phi, psi, ..., x_ca, y_ca, z_ca)
                rsa = dssp[dssp_key][3]
                rsa_values.append(rsa)
            except KeyError:
                # Residue might be missing in DSSP output (e.g., disorder, not resolved in structure)
                # Or numbering mismatch if PDB doesn't have this exact residue number for the chain.
                rsa_values.append(float('nan')) # Use NaN for missing residues

        print(f"Extracted {len(rsa_values)} RSA values from PDB {pdb_id}, Chain {chain_id}.")

        # --- Alignment of RSA values with MSA ---
        # The number of RSA values must match the number of positions in the MSA (len(conservation_df)).
        # This is a critical step and highly dependent on how the query sequence maps to the PDB structure.
        # The original script had a simple adjustment for a specific length mismatch.
        # A robust solution needs a proper sequence alignment between MSA reference and PDB sequence.
        # For now, we'll use a simplified approach: if lengths mismatch, pad/truncate.

        msa_length = len(conservation_df)
        if len(rsa_values) != msa_length:
            print(f"Warning: RSA value count ({len(rsa_values)}) differs from MSA length ({msa_length}).")
            print("This may indicate a mismatch between your query sequence and the PDB structure's residue range.")
            # Simple adjustment: pad with NaN or truncate if necessary.
            if len(rsa_values) < msa_length:
                # Assuming missing values are at the end (could be start or internal)
                rsa_values.extend([float('nan')] * (msa_length - len(rsa_values)))
                print(f"Padded RSA values to {msa_length} with NaNs.")
            else: # rsa_values > msa_length
                rsa_values = rsa_values[:msa_length]
                print(f"Truncated RSA values to {msa_length}.")

        # Add to DataFrame
        conservation_df['RSA'] = rsa_values
        # Burial score: 1 - RSA. Fill NaN RSA with 0.5 (neutral) before calculating burial.
        # A common RSA threshold for buried is < 0.2 or 0.25. Max RSA is typically > 1 for some residues like Gly.
        # Normalizing RSA to 0-1 before calculating burial score might be needed if values exceed 1.
        # For simplicity, we use raw RSA. If RSA can be > 1, burial score can be < 0.
        conservation_df['Burial_Score'] = 1 - conservation_df['RSA'].fillna(0.5)

        print("Added RSA and Burial Score to the DataFrame.")
        print("\nDataFrame with Conservation and Burial Scores (first 5 rows):")
        print(conservation_df.head(5))

    except FileNotFoundError: # For PDBParser
        print(f"ERROR: PDB file '{pdb_filepath}' could not be parsed (might be corrupted or wrong format).")
    except Exception as e:
        print(f"An unexpected error occurred during burial score calculation: {e}")
    return conservation_df


def calculate_p_score(original_aa: str, new_aa: str, matrix_name: str = "BLOSUM62") -> float | None:
    """
    Calculates a normalized physicochemical score (P-score) for an amino acid mutation.
    The score is normalized to a 0-1 scale, where 1 represents the most radical change
    (most negative or smallest BLOSUM score) and 0 represents the least radical
    (most positive BLOSUM score or self-mutation).
    Returns None if matrix cannot be loaded or AAs are invalid.
    """
    try:
        subs_matrix = substitution_matrices.load(matrix_name)
    except FileNotFoundError:
        print(f"Error: Substitution matrix '{matrix_name}' not found by BioPython.")
        return None
    except Exception as e:
        print(f"Error loading substitution matrix '{matrix_name}': {e}")
        return None


    if original_aa == new_aa or original_aa == '*' or new_aa == '*': # '*' for stop codon
        return 0.0  # No physicochemical cost for self-mutation or involving stop codon here

    # Get min/max scores from the matrix for normalization
    # Some matrices like BLOSUM62 store values as (AA1, AA2) and (AA2, AA1)
    # We need to iterate through its actual values to find true min/max
    all_scores = []
    for aa1_key in subs_matrix.alphabet: # Iterate through unique amino acids in matrix
        for aa2_key in subs_matrix.alphabet:
            try:
                all_scores.append(subs_matrix[aa1_key, aa2_key])
            except KeyError: # Should not happen if iterating alphabet
                pass

    if not all_scores:
        print(f"Warning: No scores found in matrix {matrix_name}. Cannot normalize.")
        return 0.5 # Return neutral if matrix is empty or unusual

    min_score = min(all_scores)
    max_score = max(all_scores)
    score_range = max_score - min_score

    if score_range == 0: # All scores in matrix are the same (unlikely for BLOSUM)
        return 0.5 # Neutral score

    try:
        # BLOSUM matrices are symmetrical; (A,R) score is same as (R,A)
        # substitution_matrices handles this: subs_matrix['A','R'] or subs_matrix[('A','R')]
        score = subs_matrix[original_aa, new_aa]
    except KeyError:
        # If AA pair not in matrix (e.g., non-standard AA, or if matrix doesn't include one)
        # Try reverse pair, though Bio.Align.substitution_matrices usually handles this.
        try:
            score = subs_matrix[new_aa, original_aa]
        except KeyError:
            print(f"Warning: Amino acid pair ({original_aa}, {new_aa}) not found in {matrix_name}. Returning neutral P-score.")
            return 0.5 # Neutral score for unlisted pairs

    # Normalize: (max_score - current_score) / range
    # This makes highly negative (radical) scores approach 1,
    # and highly positive (conservative) scores approach 0.
    normalized_p_score = (max_score - score) / score_range
    return normalized_p_score

def test_p_score_calculation():
    """Tests the P-score calculation with examples."""
    print("\n--- Testing Physicochemical Score (P-score) Calculation (BLOSUM62) ---")
    # Example 1: Conservative mutation (Leucine to Isoleucine) - expect low P-score
    original_L, new_I = 'L', 'I'
    p_score_LI = calculate_p_score(original_L, new_I)
    if p_score_LI is not None:
        print(f"P-score for conservative change ({original_L} -> {new_I}): {p_score_LI:.3f} (expect low)")

    # Example 2: Radical mutation (Tryptophan to Glycine) - expect high P-score
    original_W, new_G = 'W', 'G'
    p_score_WG = calculate_p_score(original_W, new_G)
    if p_score_WG is not None:
        print(f"P-score for radical change ({original_W} -> {new_G}): {p_score_WG:.3f} (expect high)")

    # Example 3: Self-mutation
    original_A, new_A = 'A', 'A'
    p_score_AA = calculate_p_score(original_A, new_A)
    if p_score_AA is not None:
        print(f"P-score for self-mutation ({original_A} -> {new_A}): {p_score_AA:.3f} (expect 0.0)")

    # Example 4: Invalid AA
    original_X, new_Y = 'X', 'A' # X is often unknown
    p_score_XA = calculate_p_score(original_X, new_Y)
    if p_score_XA is not None: # Will likely print warning and return 0.5
        print(f"P-score for mutation involving unknown AA ({original_X} -> {new_Y}): {p_score_XA:.3f} (expect 0.5 or error)")


# --- Main Execution ---
def main():
    parser = argparse.ArgumentParser(description="Evolutionary Data Retrieval and Analysis Pipeline.")
    parser.add_argument("--query_sequence", type=str, default=None,
                        help="Custom query protein sequence string. If None, uses default SARS-CoV-2 RBD subsequence.")
    parser.add_argument("--query_id", type=str, default="SARS-CoV-2_RBD_Query",
                        help="ID for the query sequence.")
    parser.add_argument("--query_desc", type=str, default="Partial Spike protein receptor-binding domain",
                        help="Description for the query sequence.")

    parser.add_argument("--skip_blast", action="store_true", help="Skip BLAST search and try to use existing blast_results.xml.")
    parser.add_argument("--blast_xml_input", type=str, default=OUTPUT_BLAST_XML_FILE,
                        help="Path to BLAST XML input file if skipping BLAST search.")

    parser.add_argument("--skip_mafft", action="store_true", help="Skip MAFFT alignment and try to use existing aligned_sequences.fasta.")
    parser.add_argument("--aligned_fasta_input", type=str, default=OUTPUT_ALIGNED_FASTA_FILE,
                        help="Path to aligned FASTA input file if skipping MAFFT.")

    parser.add_argument("--pdb_id", type=str, default=DEFAULT_PDB_ID, help="PDB ID for structural analysis.")
    parser.add_argument("--pdb_file_path", type=str, default=None,
                        help=f"Direct path to PDB file. If provided, fetching by PDB ID is skipped. Default: <pdb_id>.pdb")
    parser.add_argument("--pdb_chain_id", type=str, default=DEFAULT_PDB_CHAIN_ID, help="Chain ID in PDB for DSSP.")
    parser.add_argument("--dssp_start_res", type=int, default=DEFAULT_DSSP_START_RES, help="Start residue number in PDB chain for DSSP.")
    parser.add_argument("--dssp_end_res", type=int, default=DEFAULT_DSSP_END_RES, help="End residue number in PDB chain for DSSP.")
    parser.add_argument("--dssp_path", type=str, default="mkdssp", help="Path to mkdssp executable.")

    parser.add_argument("--skip_conservation", action="store_true", help="Skip conservation score calculation.")
    parser.add_argument("--skip_burial", action="store_true", help="Skip burial score calculation.")
    parser.add_argument("--skip_pscore_test", action="store_true", help="Skip P-score calculation tests.")


    args = parser.parse_args()

    # 1. Define Query Sequence
    print("--- Initializing Query Sequence ---")
    query_seq_string = args.query_sequence
    if query_seq_string is None:
        query_seq_string = get_rbd_subsequence(DEFAULT_RBD_SEQUENCE_FULL, DEFAULT_RBD_SUBSEQUENCE_SLICE)
        print(f"Using default SARS-CoV-2 RBD subsequence (residues {DEFAULT_RBD_SUBSEQUENCE_SLICE[0]+1}-{DEFAULT_RBD_SUBSEQUENCE_SLICE[1]}).")

    query_record = create_seq_record(query_seq_string, args.query_id, args.query_desc)
    print(f"Query: {query_record.id}, Length: {len(query_record.seq)} aa")
    print(f"Sequence: {query_record.seq[:30]}...")

    # 2. Run BLAST Search (or use existing results)
    blast_xml_to_parse = args.blast_xml_input
    if not args.skip_blast:
        success = run_blast_search(query_record, output_xml_file=OUTPUT_BLAST_XML_FILE)
        if not success:
            print("BLAST search failed. Attempting to use existing XML if available, or exiting.")
            if not os.path.exists(OUTPUT_BLAST_XML_FILE):
                print(f"No existing '{OUTPUT_BLAST_XML_FILE}' found. Cannot proceed without BLAST results.")
                return
        blast_xml_to_parse = OUTPUT_BLAST_XML_FILE
    else:
        print(f"Skipping BLAST search. Using specified XML: '{blast_xml_to_parse}'")
        if not os.path.exists(blast_xml_to_parse):
            print(f"ERROR: Specified BLAST XML file '{blast_xml_to_parse}' not found. Cannot proceed.")
            return

    # 3. Parse BLAST XML
    homologous_sequences = parse_blast_xml(blast_xml_to_parse, query_record)
    if len(homologous_sequences) <= 1 and not args.skip_mafft: # Only query sequence, no homologs
        print("No homologous sequences found or extracted. MSA and further analysis might be limited.")
        # Depending on desired behavior, one might choose to stop or continue with just the query.
        # For now, we'll allow it to proceed to MAFFT if not skipped (MAFFT with 1 seq is trivial).

    # 4. Generate MSA with MAFFT (or use existing alignment)
    msa = None
    aligned_fasta_to_use = args.aligned_fasta_input
    if not args.skip_mafft:
        if not homologous_sequences: # Should not happen if query_record was added
            print("Critical error: No sequences available for MAFFT. Skipping MSA.")
        else:
            msa = generate_msa_with_mafft(homologous_sequences,
                                          unaligned_file=OUTPUT_UNALIGNED_FASTA_FILE,
                                          aligned_file=OUTPUT_ALIGNED_FASTA_FILE)
            if msa is None:
                print("MAFFT alignment failed. Attempting to use existing aligned FASTA if available.")
                if not os.path.exists(OUTPUT_ALIGNED_FASTA_FILE):
                    print(f"No existing '{OUTPUT_ALIGNED_FASTA_FILE}'. Cannot proceed with conservation/burial scores.")
                    return
            aligned_fasta_to_use = OUTPUT_ALIGNED_FASTA_FILE # Use the one just created (or attempted)
    else:
        print(f"Skipping MAFFT. Using specified aligned FASTA: '{aligned_fasta_to_use}'")

    if msa is None: # If MAFFT was skipped or failed, try to load from specified/default file
        if os.path.exists(aligned_fasta_to_use):
            try:
                msa = AlignIO.read(aligned_fasta_to_use, "fasta")
                print(f"Successfully loaded MSA from '{aligned_fasta_to_use}'.")
                print(f"MSA has {len(msa)} sequences, length {msa.get_alignment_length()}.")
            except Exception as e:
                print(f"Error loading MSA from '{aligned_fasta_to_use}': {e}")
                print("Cannot proceed without a valid MSA for conservation/burial scores.")
                return
        else:
            print(f"ERROR: Aligned FASTA file '{aligned_fasta_to_use}' not found. Cannot proceed.")
            return

    # 5. Calculate Conservation Scores
    conservation_df = pd.DataFrame()
    if not args.skip_conservation:
        if msa:
            conservation_df = calculate_conservation_scores(msa)
        else:
            print("No MSA available, skipping conservation score calculation.")
    else:
        print("Skipping conservation score calculation.")

    # 6. Calculate Burial Scores (requires PDB and DSSP)
    final_df = conservation_df.copy() # Start with conservation data
    if not args.skip_burial:
        if conservation_df.empty and not args.skip_conservation : # If conservation was meant to run but failed
            print("Conservation scores are needed for adding burial scores. Skipping burial score calculation.")
        else:
            # Fetch PDB file
            pdb_filepath_to_use = args.pdb_file_path
            if not pdb_filepath_to_use: # If no direct path given, fetch by ID
                 pdb_filepath_to_use = fetch_pdb_file(args.pdb_id, download_dir=".") # Download to current dir

            if pdb_filepath_to_use and os.path.exists(pdb_filepath_to_use):
                # Diagnostic: Inspect PDB for chain residue info
                inspect_pdb_for_chain_residues(pdb_filepath_to_use, args.pdb_chain_id, args.dssp_path)

                # Calculate burial scores
                final_df = calculate_burial_scores(conservation_df, args.pdb_id,
                                                   pdb_filepath_to_use, args.pdb_chain_id,
                                                   args.dssp_start_res, args.dssp_end_res,
                                                   args.dssp_path)
            else:
                print(f"PDB file for {args.pdb_id} not available (path: {pdb_filepath_to_use}). Skipping burial score calculation.")
    else:
        print("Skipping burial score calculation.")

    # Display final combined DataFrame (if any analysis was done)
    if not final_df.empty:
        print("\n--- Final Combined DataFrame (first 10 rows) ---")
        print(final_df.head(10))
        # Example of saving:
        # final_df.to_csv("protein_analysis_scores.csv", index=False)
        # print("\nSaved combined scores to 'protein_analysis_scores.csv'")

    # 7. Test Physicochemical Score Calculation
    if not args.skip_pscore_test:
        test_p_score_calculation()
        print("\nNote: The P-score calculation for all mutations as described in the original")
        print("script's last comment ('loop through all 187 positions...') is not implemented")
        print("as a full loop here, but the `calculate_p_score` function is available.")

    print("\n--- Evolutionary Data Retrieval Pipeline Finished ---")

if __name__ == "__main__":
    # Note: To run this script effectively:
    # 1. Ensure BioPython and pandas are installed (`pip install biopython pandas`).
    # 2. Ensure MAFFT is installed and in PATH (`sudo apt-get install mafft`).
    # 3. Ensure DSSP (mkdssp) is installed and in PATH (`sudo apt-get install dssp`).
    # 4. Adjust file paths and parameters as needed, especially for PDB/DSSP sections
    #    to match your specific query sequence and chosen PDB structure.
    main()
