# Guide: Accessing MHC Sequence Databases

This guide provides a conceptual overview and a practical Python script for programmatically accessing Major Histocompatibility Complex (MHC) sequence data from primary biological databases.

---

## Table of Contents
1.  [Primary Data Sources](#1-primary-data-sources)
2.  [How to Use the Python Script](#2-how-to-use-the-python-script)
3.  [Script Methodology Explained](#3-script-methodology-explained)
    *   [NCBI Entrez API](#ncbi-entrez-api)
    *   [IMGT/IPD FTP Download](#imgt-ipd-ftp-download)
    *   [IEDB API Interaction](#iedb-api-interaction)
4.  [The Complete Python Script](#4-the-complete-python-script)
5.  [Next Steps: Further Analysis](#5-next-steps-further-analysis)

---

## 1. Primary Data Sources

This guide and script focus on accessing data from the world's most prominent MHC and immunology databases:

*   **[IMGT/HLA Database](https://www.ebi.ac.uk/ipd/imgt/hla/)**: The gold standard for Human Leukocyte Antigen (HLA) sequences and nomenclature.
*   **[IPD-MHC Database](https://www.ebi.ac.uk/ipd/mhc/)**: The official repository for non-human MHC sequences.
*   **[NCBI (GenBank/RefSeq)](https://www.ncbi.nlm.nih.gov/)**: A vast, general-purpose sequence database containing a large number of MHC sequences.
*   **[IEDB (Immune Epitope Database)](https://www.iedb.org/)**: Primarily for immune epitopes, but offers powerful tools for MHC binding prediction and analysis.

---

## 2. How to Use the Python Script

### Prerequisites
1.  Install Python 3.
2.  Install the necessary libraries:
    ```bash
    pip install biopython requests pandas
    ```

### Quickstart
1.  Save the code from the [complete script section](#4-the-complete-python-script) below as a Python file (e.g., `fetch_mhc.py`).
2.  **Crucially, edit the script** and replace `"your_email@example.com"` with your actual email address. NCBI's Entrez API requires this.
3.  Run the script from your terminal:
    ```bash
    python fetch_mhc.py
    ```
4.  The script will execute the default examples, querying NCBI for human and bovine MHC sequences and saving the results to `.fasta` files in the current directory.
5.  To customize, edit the `if __name__ == "__main__":` block at the bottom of the script. You can change the search terms, organisms, or uncomment the FTP download sections.

---

## 3. Script Methodology Explained

The script is divided into functions for interacting with different databases.

### NCBI Entrez API
*   **Method:** Uses the `Bio.Entrez` module from Biopython to programmatically query and fetch data. This is ideal for targeted, flexible searches.
*   **Functions:**
    *   `query_ncbi_mhc()`: Constructs a search term to find relevant MHC sequences for a given gene and organism.
    *   `fetch_sequences_from_ncbi()`: Takes a list of IDs from the search and downloads the corresponding sequence data in a specified format (e.g., FASTA).

### IMGT IPD FTP Download
*   **Method:** Connects to the EBI's public FTP server to download entire curated datasets. This is the best approach for getting comprehensive, gold-standard allele sets.
*   **Function:**
    *   `download_imgt_data()`: Uses Python's built-in `ftplib` to navigate to the correct directory on the FTP server and download specific files (e.g., `hla_nuc.fasta`).
> **⚠️ Important:** FTP server paths can change. Always verify the correct download paths on the official [IMGT/HLA](https://www.ebi.ac.uk/ipd/imgt/hla/download/) and [IPD-MHC](https://www.ebi.ac.uk/ipd/mhc/download/) websites before running the FTP functions.

### IEDB API Interaction
*   **Method:** Demonstrates how to send a request to one of IEDB's web tools APIs.
*   **Function:**
    *   `query_iedb_mhc_allele_info()`: Shows an example of submitting peptides to the MHC-I binding prediction tool. This highlights that IEDB is best used for *analysis and prediction* after you have obtained MHC sequences from a primary source like IMGT or NCBI.

---

## 4. The Complete Python Script

```python
#!/usr/bin/env python3

import os
import ftplib
import requests
import pandas as pd
from Bio import Entrez, SeqIO
from io import StringIO

# --- Configuration ---
# ALWAYS provide your email to NCBI
Entrez.email = "your_email@example.com"
# Optional: Get an API key from your NCBI account for higher request rates
Entrez.api_key = None # "YOUR_API_KEY_HERE"

# --- Helper Functions ---
def safe_filename(name):
    """Creates a safe filename by replacing non-alphanumeric characters."""
    return "".join(c if c.isalnum() or c in ('.', '_') else '_' for c in name)

# --- 1. NCBI (GenBank/RefSeq) ---
def query_ncbi_mhc(gene_name, organism, database="nucleotide", retmax=10):
    """Queries NCBI for MHC sequences and returns a list of Entrez IDs."""
    print(f"\n[NCBI] Searching for '{gene_name}' in '{organism}' ({database})...")
    search_term = f'"{gene_name}"[Gene Name] OR "{gene_name}"[All Fields] AND "{organism}"[Organism] AND (MHC OR "major histocompatibility complex" OR HLA)[All Fields]'
    
    try:
        handle = Entrez.esearch(db=database, term=search_term, retmax=retmax, idtype="acc")
        record = Entrez.read(handle)
        handle.close()
        print(f"[NCBI] Found {record['Count']} records. Returning IDs for up to {retmax}.")
        return record["IdList"]
    except Exception as e:
        print(f"[NCBI] Error during search: {e}")
        return []

def fetch_sequences_from_ncbi(ids, database="nucleotide", output_format="fasta", filename_prefix="ncbi_mhc"):
    """Fetches sequences from NCBI and saves them to a file."""
    if not ids:
        print("[NCBI] No IDs provided to fetch.")
        return None

    print(f"[NCBI] Fetching {len(ids)} sequences in {output_format} format...")
    try:
        handle = Entrez.efetch(db=database, id=ids, rettype=output_format, retmode="text")
        content = handle.read()
        handle.close()
        
        if not content.strip():
            print("[NCBI] No sequence data returned.")
            return None

        output_filename = f"{filename_prefix}_{database}.{output_format}"
        with open(output_filename, "w") as f:
            f.write(content)
        print(f"[NCBI] Sequences saved to {output_filename}")
        return output_filename
    except Exception as e:
        print(f"[NCBI] Error fetching sequences: {e}")
        return None

# --- 2. IMGT/IPD (FTP Download) ---
def download_imgt_data(ftp_host, ftp_path, local_dir="imgt_data", filename_pattern=None):
    """Downloads curated datasets from an IMGT FTP site."""
    print(f"\n[IMGT/IPD] Connecting to FTP: {ftp_host}{ftp_path}")
    os.makedirs(local_dir, exist_ok=True)
    
    try:
        with ftplib.FTP(ftp_host) as ftp:
            ftp.login()  # Anonymous login
            ftp.cwd(ftp_path)
            
            filenames = ftp.nlst()
            for filename in filenames:
                if filename_pattern and not filename.endswith(filename_pattern.strip("*")):
                    continue
                
                local_filepath = os.path.join(local_dir, filename)
                print(f"[IMGT/IPD] Downloading {filename}...")
                with open(local_filepath, 'wb') as f_local:
                    ftp.retrbinary(f'RETR {filename}', f_local.write)
            print(f"[IMGT/IPD] Download complete. Files saved in '{local_dir}'.")
    except ftplib.all_errors as e:
        print(f"[IMGT/IPD] FTP Error: {e}")

# --- 3. IEDB (API Interaction Example) ---
def predict_mhc_binding_iedb(peptides, mhc_allele_name="HLA-A*02:01"):
    """Submits peptides to the IEDB MHC-I binding prediction tool."""
    print(f"\n[IEDB] Submitting peptide binding prediction for allele {mhc_allele_name}...")
    iedb_api_url = "https://tools.iedb.org/tools_api/mhci/"
    sequence_text = "\\n".join(peptides)
    
    payload = {
        'method': 'recommended',
        'sequence_text': sequence_text,
        'allele': mhc_allele_name,
    }
    
    try:
        response = requests.post(iedb_api_url, data=payload)
        response.raise_for_status()
        df = pd.read_csv(StringIO(response.text), sep='\\t')
        print("[IEDB] Prediction successful.")
        return df
    except requests.exceptions.RequestException as e:
        print(f"[IEDB] API request error: {e}")
        return None

# --- 4. Basic Analysis ---
def analyze_fasta_file(fasta_file):
    """Performs basic analysis on a FASTA file."""
    if not fasta_file or not os.path.exists(fasta_file):
        print(f"\n[Analysis] File not found: {fasta_file}")
        return

    print(f"\n[Analysis] Analyzing: {fasta_file}")
    count = sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))
    print(f"[Analysis] Found {count} sequences.")

# --- Main Execution ---
if __name__ == "__main__":
    # --- NCBI Example ---
    hla_ids = query_ncbi_mhc(gene_name="HLA-A", organism="Homo sapiens", retmax=5)
    if hla_ids:
        hla_fasta_file = fetch_sequences_from_ncbi(hla_ids, filename_prefix="human_hla_a")
        analyze_fasta_file(hla_fasta_file)

    # --- IMGT/IPD FTP Download Example (Commented Out) ---
    print("\n--- IMGT/IPD FTP Download Note ---")
    print("The FTP download examples are commented out by default.")
    print("Please uncomment and VERIFY THE FTP PATHS on the official EBI/IPD websites before running.")
    # download_imgt_data(
    #     ftp_host="ftp.ebi.ac.uk",
    #     ftp_path="/pub/databases/ipd/imgt/hla/fasta/", # Verify this path!
    #     local_dir="imgt_hla_data",
    #     filename_pattern="hla_nuc.fasta"
    # )
    # analyze_fasta_file("imgt_hla_data/hla_nuc.fasta")
    
    # --- IEDB API Example ---
    example_peptides = ['SLYNTVATL', 'NLVPMVATV']
    binding_predictions = predict_mhc_binding_iedb(peptides=example_peptides, mhc_allele_name="HLA-A*02:01")
    if binding_predictions is not None:
        print("[IEDB] Binding Prediction Results:")
        print(binding_predictions)

    print("\nScript finished.")
```

---

## 5. Next Steps: Further Analysis

Once you have downloaded the MHC sequences, you can proceed with a variety of bioinformatic analyses:

*   **Sequence Alignment:** Use tools like MAFFT, ClustalW, or MUSCLE to align sequences and identify conserved regions or variable sites.
*   **Phylogenetic Analysis:** Build evolutionary trees using PhyML, RAxML, or IQ-TREE to understand the relationships between different alleles.
*   **Variant Analysis:** Identify single nucleotide polymorphisms (SNPs) or other variations within your sequence set.
*   **Structural Modeling:** Use protein sequences to predict 3D structures with tools like AlphaFold to analyze the peptide-binding groove.
