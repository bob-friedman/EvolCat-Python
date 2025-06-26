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

# Install necessary bioinformatics tools
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

# --- Fetch Accession Numbers ---
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

# --- Fetch GenBank Records & Extract Spike Proteins ---
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

# --- Step 3: Run IQ-TREE for phylogeny and ASR ---
print(f"\nRunning IQ-TREE for phylogenetic inference and ASR...")
# Step 1: Build the tree and find the best model, but do NOT run ASR yet.
!iqtree -s sequence_aligned.fas -m MFP -T AUTO -redo -safe
# Step 2: Now, run ASR using the tree and model that were just created.
# This explicitly tells iqtree which tree to use for the reconstruction.
!iqtree -s sequence_aligned.fas -t sequence_aligned.fas.treefile -asr --redo
# Above two line are a problem: mismatch between the node names in the phylogenetic tree file (`.treefile`) and
# the sequence names in the ancestral state file (`.state`).
# !iqtree -s sequence_aligned.fas -m MFP -asr -T AUTO -redo -safe

print("\nIQ-TREE run complete.")
!ls -lh sequence_aligned.fas*

# 5. Create Ancestor-Descendant Pairs (Definitive Method using Pandas)
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

# 2. Get the sequences for the tips (leaves) of the tree.
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

### **SASA vs. Attention Score Analysis**

# Hypothesis: The more a residue is exposed on the surface (higher SASA), the more attention the model will pay to it.

# This code block will:
# 1.  Install `FreeSASA`, a tool to calculate solvent accessibility.
# 2.  Run `FreeSASA` on our PDB file (`6VXX.pdb`).
# 3.  Parse the output to get per-residue SASA scores.
# 4.  Align these scores with the per-residue attention scores we already calculated.
 #5.  Generate the scatter plot and calculate the statistical correlation.

# --- New Cell (Corrected Version) ---

# 0. Set up Colab with MyDrive
from google.colab import drive
import os

# Mount Google Drive and set the working directory
drive.mount('/content/drive')
os.chdir('/content/drive/MyDrive')
print(f"Current working directory: {os.getcwd()}")

# 1. Install all necessary Python libraries with pip
print("Installing dependencies (freesasa, biopython, scipy)...")
!pip install -q freesasa biopython scipy

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser
from scipy.stats import pearsonr
import os
import freesasa

print("Dependencies installed successfully.")

# 2. Define filenames and load existing attention data
PDB_FILE = "6VXX.pdb"

# --- This is a critical check ---
if 'per_residue_attention' not in locals():
    print("Warning: 'per_residue_attention' not found in memory.")
    print("This is expected if running this cell standalone.")
    print("For demonstration, a random array will be used. In a real run, load this from the model.")
    per_residue_attention = np.random.rand(1300) * 10

# 3. Calculate Solvent Accessible Surface Area (SASA) using the Python API
if not os.path.exists(PDB_FILE):
    raise FileNotFoundError(f"{PDB_FILE} not found. Please ensure it is in the correct directory.")

print(f"Calculating SASA for {PDB_FILE} using the freesasa Python API...")
structure = freesasa.Structure(PDB_FILE)
result = freesasa.calc(structure)

# 4. Extract per-residue SASA scores with robust error handling
print("Extracting per-residue SASA scores...")
sasa_scores = {}
residue_areas = result.residueAreas()

if 'A' in residue_areas:
    chain_a_areas = residue_areas['A']
    for res_string, area_obj in chain_a_areas.items():
        try:
            parts = res_string.split()
            if not parts: continue
            res_id_int = int(parts[-1])
            sasa_scores[res_id_int] = area_obj.total
        except (IndexError, ValueError):
            print(f"  - Warning: Could not parse residue ID from key '{res_string}'. Skipping.")
            continue
else:
    print("Warning: Chain 'A' was not found in the SASA calculation results.")

print(f"Extracted and summed SASA scores for {len(sasa_scores)} residues from Chain A.")

# 5. Align SASA and Attention scores into a single DataFrame for analysis
data_points = []
for res_id, sasa_value in sasa_scores.items():
    score_index = res_id
    if score_index < len(per_residue_attention):
        attention_value = per_residue_attention[score_index]
        data_points.append({'ResidueID': res_id, 'SASA': sasa_value, 'Attention': attention_value})

df_analysis = pd.DataFrame(data_points)
print(f"Aligned {len(df_analysis)} data points for statistical analysis.")

# 6. Perform Correlation Analysis and Visualization
if not df_analysis.empty:
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(10, 7))

    # The keyword for scatter point edge width is 'linewidths' (plural)
    sns.regplot(data=df_analysis, x='SASA', y='Attention', ax=ax,
                scatter_kws={'alpha':0.6, 's':25, 'edgecolor':'w', 'linewidths':0.5},
                line_kws={'color':'#e41a1c', 'linestyle':'--'})

    correlation, p_value = pearsonr(df_analysis['SASA'], df_analysis['Attention'])

    ax.set_title('Correlation between Solvent Accessibility and Model Attention', fontsize=16, pad=20)
    ax.set_xlabel('Total Solvent Accessible Surface Area per Residue (SASA, in Å²)', fontsize=12)
    ax.set_ylabel('Transformer Attention Score (arbitrary units)', fontsize=12)

    stats_text = f'Pearson r = {correlation:.3f}\np-value = {p_value:.3e}'
    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', bbox=dict(boxstyle='round,pad=0.5', fc='aliceblue', alpha=0.8))

    plt.show()
else:
    print("Could not generate plot because no data points were aligned.")

'''
# Result:
*   Pearson r = -0.029: This value is close to zero, indicating that there is virtually no
    linear relationship between the Solvent Accessible Surface Area (SASA) and the attention score the model
    assigns to a residue.
*   p-value = 0.3727: This value is very high (much greater than the standard significance threshold of 0.05).
      This tells us that the tiny negative correlation observed is not statistically significant and is very likely due to random chance.
''';

### Retrieving Data from the IEDB Database

'''
1.  Go to the [IEDB Database Export page](https://www.iedb.org/database_export_v3.php).
2.  Find the "B cell epitopes" section.
3.  Click the link to download `bcell_full_v3.zip`.
4.  Unzip the file to get `bcell_full_v3.csv`.
5.  Upload this CSV file to your Colab environment.

Once the file is uploaded, we can use the following code to load and process it.
'''

import os
from google.colab import drive

# Mount Google Drive to save our work
drive.mount('/content/drive')
os.chdir('/content/drive/MyDrive')
print(f"Current working directory: {os.getcwd()}")

!unzip -x bcell_full_v3_single_file.zip

# Diagnostic #1 for Fields in the IEDB Database

import pandas as pd

file_path = 'bcell_full_v3.csv'

try:
    # Load the file
    df_iedb = pd.read_csv(file_path, header=1, low_memory=False)
    print(f"Successfully loaded '{file_path}'.")

    # --- Diagnostic Section ---
    print("\n--- Investigating Filter Values ---")

    # 1. Check values in the 'Source Organism' column for anything containing 'coronavirus'
    organism_options = df_iedb[df_iedb['Source Organism'].str.contains('coronavirus', case=False, na=False)]['Source Organism'].unique()
    print("\nPotential 'Source Organism' values:")
    print(organism_options)

    # 2. Check values in the 'Name' column for anything containing 'Spike'
    antigen_options = df_iedb[df_iedb['Name'].str.contains('Spike', case=False, na=False)]['Name'].unique()
    print("\nPotential 'Name' (Antigen) values:")
    print(antigen_options)

    # 3. Check all unique values in the 'Type' column
    type_options = df_iedb['Type'].unique()
    print("\nUnique 'Type' values:")
    print(type_options)
    print("--- End of Investigation ---")

except FileNotFoundError:
    print(f"File '{file_path}' not found. Please upload it to the Colab environment.")
except Exception as e:
    print(f"An error occurred: {e}")

# Diagnostic #2 for Fields in the IEDB Database

import pandas as pd

file_path = 'bcell_full_v3.csv'

try:
    # Load the file
    df_iedb = pd.read_csv(file_path, header=1, low_memory=False)
    print(f"Successfully loaded '{file_path}'.")

    # --- Targeted Diagnostic ---
    # Check values in the 'Source Molecule' column for anything containing 'Spike'
    print("\n--- Investigating 'Source Molecule' column ---")
    antigen_options = df_iedb[df_iedb['Source Molecule'].str.contains('Spike', case=False, na=False)]['Source Molecule'].unique()
    print("\nPotential 'Source Molecule' values for Spike:")
    print(antigen_options)
    print("--- End of Investigation ---")

except FileNotFoundError:
    print(f"File '{file_path}' not found. Please upload it to the Colab environment.")
except Exception as e:
    print(f"An error occurred: {e}")

# Diagnostic #3 for Fields in the IEDB Database

import pandas as pd

file_path = 'bcell_full_v3.csv'

try:
    # Load the file
    df_iedb = pd.read_csv(file_path, header=1, low_memory=False)
    print(f"Successfully loaded '{file_path}'.")

    # --- Final Diagnostic Step ---
    print("\n--- Finding Molecule Name for the Correct Organism ---")

    # First, filter for ONLY the organism
    df_sars2_organism = df_iedb[df_iedb['Source Organism'] == 'Severe acute respiratory syndrome coronavirus 2'].copy()

    # Now, from that subset, find the unique molecule names used
    molecule_names = df_sars2_organism['Source Molecule'].unique()

    print(f"Found {len(df_sars2_organism)} entries for the organism 'Severe acute respiratory syndrome coronavirus 2'.")
    print("\nMolecule names associated with this organism:")
    print(molecule_names)
    print("\n--- End of Investigation ---")


except FileNotFoundError:
    print(f"File '{file_path}' not found. Please upload it to your environment.")
except Exception as e:
    print(f"An error occurred: {e}")

# Epitope Analysis of Surface Residues

import pandas as pd
import numpy as np
from Bio.PDB import PDBParser, PDBList
from Bio.PDB.SASA import ShrakeRupley
from scipy.stats import mannwhitneyu

print("--- Starting Full B-Cell Epitope Analysis ---")

try:
    # --- PART 1: Load and filter the IEDB epitope data ---
    print("\n[1/4] Loading and filtering B-cell epitope data...")
    df_iedb = pd.read_csv('bcell_full_v3.csv', header=1, low_memory=False)
    df_sars2_organism = df_iedb[df_iedb['Source Organism'] == 'Severe acute respiratory syndrome coronavirus 2'].copy()
    spike_keywords = 'Spike|surface glycoprotein'
    df_sars2_spike = df_sars2_organism[df_sars2_organism['Source Molecule'].str.contains(spike_keywords, case=False, na=False)].copy()
    df_sars2_spike['Starting Position'] = pd.to_numeric(df_sars2_spike['Starting Position'], errors='coerce')
    df_sars2_spike['Ending Position'] = pd.to_numeric(df_sars2_spike['Ending Position'], errors='coerce')
    df_sars2_spike.dropna(subset=['Starting Position', 'Ending Position'], inplace=True)
    df_sars2_spike[['Starting Position', 'Ending Position']] = df_sars2_spike[['Starting Position', 'Ending Position']].astype(int)
    print(f"Loaded {len(df_sars2_spike)} epitope entries.")

    # --- PART 2: Define Residue Sets A and B ---
    print("\n[2/4] Defining epitope and surface residue sets...")
    # Create Set A (all unique residues in B-cell epitopes)
    b_cell_epitope_residues = set()
    for index, row in df_sars2_spike.iterrows():
        for residue_pos in range(row['Starting Position'], row['Ending Position'] + 1):
            b_cell_epitope_residues.add(residue_pos)
    Set_A = b_cell_epitope_residues

    # Calculate surface residues using BioPython
    pdbl = PDBList()
    pdbl.retrieve_pdb_file('6VXX', pdir='.', file_format='pdb', overwrite=True)
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('6VXX', 'pdb6vxx.ent')
    sasa_calc = ShrakeRupley()
    sasa_calc.compute(structure[0], level="R")

    # Define all surface residues, ensuring they are within our model's range (1-1300)
    surface_residues = set()
    sasa_cutoff = 1.0
    for residue in structure[0]['A']:
        res_id = residue.get_id()[1]
        # This is the crucial fix to filter out-of-range PDB residues
        if res_id <= 1300 and residue.sasa > sasa_cutoff:
            surface_residues.add(res_id)

    print(f"--- Diagnostic: Total surface residues found: {len(surface_residues)}")
    intersection_size = len(surface_residues.intersection(b_cell_epitope_residues))
    print(f"--- Diagnostic: Surface residues that are also B-cell epitopes: {intersection_size}")

    # --- Analysis: Compare Surface vs. Buried Residues ---

    # Set_A is now all surface residues
    Set_A_surface = surface_residues

    # Set_B is now all buried residues (all model positions minus surface positions)
    all_model_residues = set(attention_map.keys())
    Set_B_buried = all_model_residues.difference(Set_A_surface)

    print(f"Set A (Surface Residues): {len(Set_A_surface)} residues")
    print(f"Set B (Buried Residues): {len(Set_B_buried)} residues")

    # --- PART 4: Perform the final analysis ---
    print("\n[4/4] Performing statistical analysis (Surface vs. Buried)...")
    attention_A = [attention_map[pos] for pos in Set_A_surface if pos in attention_map]
    attention_B = [attention_map[pos] for pos in Set_B_buried if pos in attention_map]

    if not attention_A or not attention_B:
         raise ValueError("One of the new residue sets (Surface or Buried) is empty.")

    mean_attention_A = np.mean(attention_A)
    mean_attention_B = np.mean(attention_B)
    statistic, p_value = mannwhitneyu(attention_A, attention_B, alternative='greater')

    # --- Report Final Results ---
    print("\n--- ANALYSIS COMPLETE (Surface vs. Buried) ---")
    print(f"Mean Attention for Set A (Surface): {mean_attention_A:.6f}")
    print(f"Mean Attention for Set B (Buried): {mean_attention_B:.6f}")
    print(f"\nMann-Whitney U test p-value: {p_value:.6f}")

    if p_value < 0.05:
        print("\nConclusion: The result is statistically significant (p < 0.05).")
        print("The mean attention on SURFACE residues is significantly higher than on BURIED residues.")
    else:
        print("\nConclusion: The result is not statistically significant (p >= 0.05).")
        print("We cannot conclude that attention on SURFACE residues is higher than on BURIED residues.")

except FileNotFoundError as fnf:
    print(f"\nERROR: A required file was not found. {fnf}")
except Exception as e:
    print(f"\nAn unexpected error occurred: {e}")