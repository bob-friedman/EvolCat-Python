### **Evolutionary Data Retrieval**

# --- Step 1: Install necessary library ---
# BioPython is the standard library for computational biology in Python.
print("Installing BioPython...")
!pip install biopython
print("Installation complete.")

# --- Step 2: Import necessary modules and define our target sequence ---
from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Define the sequence for the SARS-CoV-2 Spike RBD (Wuhan-Hu-1, residues 331-531)
rbd_sequence_str = "NIDTFLDSVYKGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVN"
rbd_subsequence = rbd_sequence_str[330:531] # Python is 0-indexed, so this selects residues 331 to 531

# Create a SeqRecord object, which is the standard input format for BioPython tools
rbd_record = SeqRecord(
    Seq(rbd_subsequence),
    id="SARS-CoV-2_RBD",
    description="Spike protein receptor-binding domain from Wuhan-Hu-1"
)

# --- Step 3: Perform BLAST Search ---
# Send sequence to NCBI servers for comparison in the 'nr' (non-redundant) database.
# This step is computationally intensive and will take a significant amount of time.
print("\nStarting BLAST search against the NCBI 'nr' database.")
print("This may take 10-30 minutes or more, depending on server load. Please be patient.")

try:
    # `blastp` is the program for protein-protein comparison
    result_handle = NCBIWWW.qblast("blastp", "nr", rbd_record.seq)

    # --- Step 4: Save the Results ---
    # Save raw output as an XML file. This file contains all the information
    # about the sequences that were found to be similar to our target.
    with open("blast_results.xml", "w") as out_handle:
        out_handle.write(result_handle.read())

    print("\n--- Success! ---")
    print("BLAST search complete and results have been saved to 'blast_results.xml'.")
    print("You will see this file appear in the Colab file browser on the left.")

except Exception as e:
    print(f"\nAn error occurred during the BLAST search: {e}")
    print("This can sometimes happen due to network issues or NCBI server load. Please try running the script again.")

### **Parse BLAST XML and Extract Homologous Sequences**

# BioPython is the standard library for computational biology in Python
print("Installing BioPython...")
!pip install biopython
print("Installation complete.")

# --- Step 1: Import necessary modules ---
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO # To treat strings as files, useful for alignment

# --- Step 2: Define file and initialize a list for our sequences ---
blast_xml_file = "blast_results.xml"
homologous_sequences = []

# Define the sequence for the SARS-CoV-2 Spike RBD (Wuhan-Hu-1, residues 331-531)
rbd_sequence_str = "NIDTFLDSVYKGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVN"
rbd_subsequence = rbd_sequence_str[330:531] # Python is 0-indexed, so this selects residues 331 to 531

# Create a SeqRecord object, which is the standard input format for BioPython tools
rbd_record = SeqRecord(
    Seq(rbd_subsequence),
    id="SARS-CoV-2_RBD",
    description="Spike protein receptor-binding domain from Wuhan-Hu-1"
)

# Add the original query sequence to the list first
# This ensures it is part of the alignment for downstream use.
homologous_sequences.append(rbd_record)

print(f"Starting to parse '{blast_xml_file}'...")

# --- Step 3: Parse BLAST results ---
try:
    with open(blast_xml_file, "r") as result_handle:
        # The NCBIXML.parse() function is an iterator, which is memory-efficient
        # for large BLAST files. It yields one record at a time.
        blast_records = NCBIXML.parse(result_handle)

        # Since only one query, work with the first and only record.
        for blast_record in blast_records:
            print(f"\nProcessing hits for query: {blast_record.query[:30]}...")
            print(f"Found {len(blast_record.alignments)} alignments.")

            # Each 'alignment' in the record corresponds to a hit from the database.
            for alignment in blast_record.alignments:
                # An alignment can have multiple 'HSPs' (High-scoring Segment Pairs),
                # which are the specific regions that align. Use the first one.
                for hsp in alignment.hsps:
                    # Indentify 'hit' or 'subject'
                    hit_sequence_str = hsp.sbjct

                    # Aligned sequence may contain gaps ('-'). Retain.
                    # for now as they inform the alignment.
                    hit_id = alignment.hit_id
                    hit_def = alignment.hit_def.split(']')[0] + ']' # Clean up description

                    # Create SeqRecord for this hit
                    hit_record = SeqRecord(
                        Seq(hit_sequence_str),
                        id=hit_id,
                        description=hit_def
                    )
                    homologous_sequences.append(hit_record)

    print("\n--- Success! ---")
    print(f"Extracted a total of {len(homologous_sequences)} sequences (including the original query).")
    print("These sequences are now stored in the 'homologous_sequences' list.")
    print("\nFirst 5 extracted sequence IDs:")
    for seq in homologous_sequences[:5]:
        print(f"- {seq.id}: {seq.description[:70]}...")


except FileNotFoundError:
    print(f"ERROR: The file '{blast_xml_file}' was not found.")
    print("Please make sure you have successfully run the previous BLAST search script.")
except Exception as e:
    print(f"An error occurred during parsing: {e}")

### **Generate the MSA**

# --- Step 3a: Install MAFFT alignment program ---
# This uses the system package manager, not pip as MAFFT is a standalone program.
print("Installing MAFFT...")
# The '-y' flag auto-confirms the installation, and '>/dev/null' hides the long output.
!sudo apt-get -y -qq install mafft > /dev/null
print("MAFFT installation complete.")

# --- Step 3b: Import necessary modules for alignment ---
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
from io import StringIO # To treat strings as if they were files

# Define name for a temporary file to hold unaligned sequences
unaligned_fasta_file = "unaligned_sequences.fasta"

# Write list of SeqRecord objects to the FASTA file
SeqIO.write(homologous_sequences, unaligned_fasta_file, "fasta")

print(f"\nSaved {len(homologous_sequences)} unaligned sequences to '{unaligned_fasta_file}'.")
print("Running MAFFT to create the alignment. This may take a minute...")

# --- Step 3c: Run MAFFT ---
try:
    # The MafftCommandline object builds the command string for us.
    mafft_cline = MafftCommandline(input=unaligned_fasta_file)

    # Execute the command. stdout and stderr capture the output from the program.
    stdout, stderr = mafft_cline()

    # The 'stdout' variable now holds alignment in FASTA format as a single large string.
    # Parse this string directly using AlignIO.
    # The StringIO wrapper treats the string as a file, which AlignIO expects.
    msa = AlignIO.read(StringIO(stdout), "fasta")

    print("\n--- Success! ---")
    print("Multiple Sequence Alignment (MSA) created successfully.")
    print(f"The MSA has {len(msa)} sequences and an alignment length of {msa.get_alignment_length()} columns.")

    # Optional: Save the alignment to a file for inspection
    aligned_fasta_file = "aligned_sequences.fasta"
    AlignIO.write(msa, aligned_fasta_file, "fasta")
    print(f"Alignment saved to '{aligned_fasta_file}'. You can view it in the file browser.")

    # Print a small snippet of the MSA
    print("\nSnippet of the final MSA:")
    print(msa[:, :60]) # Print all rows, but only the first 60 columns

except Exception as e:
    print(f"\nAn error occurred during alignment: {e}")

### **Code to Calculate Conservation Scores**

# --- Step 4a: Import libraries ---
import pandas as pd
from collections import Counter

# Use the 'msa' object created in the previous step.
if 'msa' not in globals():
    print("MSA object not found. Please run the previous MAFFT alignment step first.")
else:
    # --- Step 4b: Calculate conservation for each column ---
    print("Calculating conservation scores for each position in the MSA...")

    conservation_scores = []
    num_sequences = len(msa)
    alignment_length = msa.get_alignment_length()

    # Get the sequence of reference protein (the first one in the alignment)
    reference_sequence = msa[0].seq

    for i in range(alignment_length):
        # Get all the amino acids in the current column (position i)
        column_residues = msa[:, i]

        # Use Counter to get the counts of each amino acid in the column
        # e.g., Counter({'N': 51})
        counts = Counter(column_residues)

        # Find the count of the most frequent amino acid
        if not counts: # Handle empty columns, though unlikely here
             most_common_count = 0
        else:
             most_common_count = counts.most_common(1)[0][1]

        # Calculate conservation score
        score = most_common_count / num_sequences

        # Retrieve amino acid from our original query at this position
        reference_aa = reference_sequence[i]

        conservation_scores.append({
            "Position": i + 1,  # Use 1-based numbering for biological convention
            "Reference_AA": reference_aa,
            "Conservation_Score": score
        })

    # --- Step 4c: Create a pandas DataFrame for easy viewing ---
    conservation_df = pd.DataFrame(conservation_scores)

    print("\n--- Success ---")
    print("Conservation scores calculated.")

    # Display the first 10 and last 10 rows of the results
    print("\nFirst 10 positions:")
    print(conservation_df.head(10))

    print("\nLast 10 positions:")
    print(conservation_df.tail(10))

    print("\nScores for specific positions (e.g., N501Y location):")
    # Position 501 in full spike is ~171 in our 187aa query
    # N501Y is N -> Y. The original is N.
    # Note: Our sequence starts at 331, but the PDB/literature often uses full spike numbering.
    # Our "Position" is relative to the start of our 187aa sequence.
    # The 'N' of N501Y is at position 171 in our sequence (501 - 330 = 171).
    print(conservation_df[conservation_df['Position'] == 171])

### **Diagnostic Code to Find Correct PDB Residue Numbers**

# Re-import and parse since the runtime was restarted.
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

print("Inspecting the PDB file to find the correct residue numbering...")

try:
    p = PDBParser(QUIET=True)
    structure = p.get_structure("7VYR", "7VYR.pdb")
    model = structure[0]
    dssp = DSSP(model, "7VYR.pdb", dssp='mkdssp')

    # Get a list of all residue keys that are actually in the DSSP results
    present_keys = list(dssp.keys())

    # Filter for keys belonging to Chain 'C'
    chain_c_keys = [key for key in present_keys if key[0] == 'C']

    if not chain_c_keys:
        print("\nERROR: No residues found for Chain C. The chain ID might be different.")
    else:
        # Extract just the residue numbers
        chain_c_residue_numbers = sorted([key[1][1] for key in chain_c_keys])

        print(f"\nFound {len(chain_c_residue_numbers)} residues for Chain C.")
        print(f"The residue numbering in the PDB file for Chain C starts at: {chain_c_residue_numbers[0]}")
        print(f"The residue numbering in the PDB file for Chain C ends at:   {chain_c_residue_numbers[-1]}")
        print("\nFirst 10 residue numbers found:", chain_c_residue_numbers[:10])

except Exception as e:
    print(f"An error occurred during inspection: {e}")

### **Burial Score (Step 5)**

print("\n--- Starting Final Step 5: Calculating Burial Score ---")
# The PDB file should already be downloaded from the previous step.

# Import necessary libraries
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import pandas as pd # Make sure pandas is available

# --- Run DSSP and parse results ---
p = PDBParser(QUIET=True)
structure = p.get_structure("7VYR", "7VYR.pdb")
model = structure[0]

# Run DSSP on the model with the correct PDB file path
dssp = DSSP(model, "7VYR.pdb", dssp='mkdssp')

# Iterate through correct residue numbers for Chain C.
start_res = 336
end_res = 521
chain_id = 'C'

rsa_values = []
for res_num in range(start_res, end_res + 1):
    try:
        a_key = (chain_id, (' ', res_num, ' '))
        rsa = dssp[a_key][3]
        rsa_values.append(rsa)
    except KeyError:
        # This handles cases where a residue within range is missing
        rsa_values.append(float('nan'))

print(f"\nExtracted {len(rsa_values)} RSA values using the correct numbering.")

# We have 186 RSA values for our 187-long alignment.
# The most likely scenario is that the first residue is missing. Add a NaN at the start.
if len(rsa_values) == 186:
    rsa_values.insert(0, float('nan')) # Add a placeholder for the first residue
    print("Adjusted RSA list to match alignment length of 187.")


# --- Add the Burial Score to DataFrame ---
if len(rsa_values) == len(conservation_df):
    conservation_df['RSA'] = rsa_values
    conservation_df['Burial_Score'] = 1 - conservation_df['RSA'].fillna(0.5)
    filter_df = conservation_df # Update the final DataFrame name

    print("\n--- Success! ---")
    print("Added RSA and Burial Score to the DataFrame using correct PDB numbering.")

    print("\nDataFrame with Conservation and Burial Scores (first 10 rows):")
    print(filter_df.head(10))

    print("\nScores for specific positions of interest:")
    # Check a known surface residue (e.g., K417, which is pos 87 in seq)
    # and a known buried residue (e.g., V367, which is pos 37 in seq)
    print(filter_df[filter_df['Position'].isin([37, 87])])

else:
    print(f"\nERROR: Mismatch in length between MSA ({len(conservation_df)}) and PDB residues ({len(rsa_values)}).")

### **Step 6: Calculating the Physicochemical Score (P)**

'''
The Physicochemical Score (`P`) quantifies how "radical" a proposed mutation is. Changing a small,
nonpolar Alanine (A) to another small, nonpolar Valine (V) is a minor change. Changing it to a large,
charged Arginine (R) is a radical change.

Use the industry-standard **BLOSUM62 matrix** to retrieve this score. This matrix contains pre-calculated scores
for every possible amino acid substitution, based on how often they are observed in conserved blocks of alignments.
*   A high positive score (e.g., +5) means the substitution is common and conservative.
*   A high negative score (e.g., -4) means the substitution is rare and radical.

Our plan:
1.  Load the BLOSUM62 matrix using BioPython.
2.  Create a function that takes an original amino acid and a new one.
3.  Inside the function, look up the score and **normalize it** to our 0-1 scale, where `1` represents the most
    radical change.
'''

# --- Step 6a: Import the substitution matrix tools (UPDATED) ---
from Bio.Align import substitution_matrices

# Load the BLOSUM62 matrix from the new location
blosum62 = substitution_matrices.load("BLOSUM62")

# --- Step 6b: Define a function to calculate the normalized P-score ---

# For normalization, use the min and max values in the matrix.
# The matrix object has convenient methods for this.
min_score = blosum62.min()
max_score = blosum62.max()
score_range = max_score - min_score

def calculate_p_score(original_aa, new_aa):
    """
    Calculates a normalized physicochemical score (0 to 1) for a mutation,
    where 1 represents the most radical change.
    """
    if original_aa == new_aa or original_aa == '*' or new_aa == '*':
        return 0.0 # A mutation to itself or involving a stop codon has no cost

    # The new matrix object handles pairs automatically.
    try:
        score = blosum62[original_aa, new_aa]
    except KeyError:
        # This can happen if one of the pair is not a standard AA.
        try:
            score = blosum62[new_aa, original_aa]
        except KeyError:
            return 0.5 # Return a neutral score if the pair is not found

    # Normalize score. Radical changes (high negative scores)
    # are close to 1, and conservative changes (high positive scores)
    # are close to 0.
    normalized_score = (max_score - score) / score_range
    return normalized_score

# --- Step 6c: Test the function with examples ---
print("Testing the Physicochemical Score (P) function:")

# Example 1: A conservative mutation (Leucine to Isoleucine)
original = 'L'
new = 'I'
p_score_conservative = calculate_p_score(original, new)
print(f"P-score for a conservative change ({original} -> {new}): {p_score_conservative:.3f}")

# Example 2: A radical mutation (Tryptophan to Glycine)
original = 'W'
new = 'G'
p_score_radical = calculate_p_score(original, new)
print(f"P-score for a radical change ({original} -> {new}): {p_score_radical:.3f}")