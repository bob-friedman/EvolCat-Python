---
layout: default # Or whatever layout your Jekyll theme uses, if any. 'default' is common.
title: Accessing MHC Sequence Databases
---

# Accessing MHC (Major Histocompatibility Complex) Sequence Databases

Connecting to MHC (Major Histocompatibility Complex) sequence databases involves a few primary sources and methods. The most prominent databases are:

1.  **IMGT/HLA Database (IPD-IMGT/HLA):** The international ImMunoGeneTics information system® for Human Leukocyte Antigens. This is the gold standard for human MHC (HLA) sequences and nomenclature.
2.  **IPD-MHC Database:** For non-human MHC sequences from various species.
3.  **NCBI (GenBank/RefSeq):** A general nucleotide and protein sequence database that contains a vast amount of MHC sequences, though not as curated specifically for MHC as IMGT.
4.  **IEDB (Immune Epitope Database):** While primarily for epitopes, it contains MHC allele information and links to sequences, and offers tools for MHC binding prediction.

Below is a Python script demonstrating programmatic access, primarily using NCBI's Entrez API (via Biopython) for broad querying, and discussing how to obtain data from IMGT/IPD (which often involves downloading datasets).

## Python Script for Accessing MHC Data

```python
#!/usr/bin/env python3

import os
import ftplib
from Bio import Entrez, SeqIO
import requests # For IEDB or other REST APIs
import pandas as pd # For handling tabular data from IEDB

# --- Configuration ---
# ALWAYS tell NCBI who you are
Entrez.email = "your_email@example.com" # Replace with your actual email
API_KEY = None # Optional: Replace with your NCBI API key for higher request rates
if API_KEY:
    Entrez.api_key = API_KEY

# --- Helper Functions ---
def safe_filename(name):
    """Creates a safe filename by replacing non-alphanumeric characters."""
    return "".join(c if c.isalnum() or c in ('.', '_') else '_' for c in name)

# --- 1. NCBI (GenBank/RefSeq) ---
def query_ncbi_mhc(gene_name, organism, database="nucleotide", retmax=10):
    """
    Queries NCBI for MHC sequences.

    Args:
        gene_name (str): e.g., "HLA-A", "MHC class I", "BoLA-DRB3"
        organism (str): e.g., "Homo sapiens", "Bos taurus"
        database (str): "nucleotide" or "protein"
        retmax (int): Maximum number of records to retrieve initially

    Returns:
        list: A list of Entrez IDs.
    """
    print(f"\n[NCBI] Searching for '{gene_name}' in '{organism}' ({database})...")
    search_term = f'"{gene_name}"[Gene Name] OR "{gene_name}"[All Fields] AND "{organism}"[Organism] AND (MHC OR "major histocompatibility complex" OR HLA)[All Fields]'
    if "HLA" in gene_name.upper() and "Homo sapiens" not in organism:
         print(f"[NCBI] Warning: HLA is human-specific. Querying for '{gene_name}' in '{organism}' might yield limited results.")
    
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
    """
    Fetches sequences from NCBI given a list of IDs and saves them to a file.

    Args:
        ids (list): List of Entrez IDs.
        database (str): "nucleotide" or "protein".
        output_format (str): "fasta", "gb" (GenBank), etc.
        filename_prefix (str): Prefix for the output FASTA file.

    Returns:
        str: Path to the saved file, or None if failed.
    """
    if not ids:
        print("[NCBI] No IDs provided to fetch.")
        return None

    print(f"[NCBI] Fetching {len(ids)} sequences in {output_format} format...")
    try:
        # Fetch in batches if many IDs (NCBI might limit URL length for GET requests)
        batch_size = 100 
        all_records_content = ""
        for i in range(0, len(ids), batch_size):
            batch_ids = ids[i:i+batch_size]
            handle = Entrez.efetch(db=database, id=batch_ids, rettype=output_format, retmode="text")
            all_records_content += handle.read()
            handle.close()
        
        if not all_records_content.strip():
            print("[NCBI] No sequence data returned from efetch.")
            return None

        output_filename = f"{filename_prefix}_{database}.{output_format}"
        with open(output_filename, "w") as f:
            f.write(all_records_content)
        print(f"[NCBI] Sequences saved to {output_filename}")
        return output_filename
    except Exception as e:
        print(f"[NCBI] Error fetching sequences: {e}")
        return None

# --- 2. IMGT/HLA and IPD-MHC (Primarily via FTP Download) ---
# These databases are best accessed by downloading their curated datasets.
# Direct programmatic querying like NCBI is limited or not straightforward for general sequence retrieval.

def download_imgt_data(ftp_host, ftp_path, local_dir="imgt_data", filename_pattern=None):
    """
    Downloads data from an IMGT FTP site.

    Args:
        ftp_host (str): e.g., "ftp.ebi.ac.uk"
        ftp_path (str): Path on the FTP server, e.g., "/pub/databases/ipd/imgt/hla/fasta/"
        local_dir (str): Local directory to save files.
        filename_pattern (str, optional): If specified, only download files matching this pattern (e.g., "*.fasta").
    """
    print(f"\n[IMGT/IPD] Connecting to FTP: {ftp_host} path: {ftp_path}")
    os.makedirs(local_dir, exist_ok=True)
    
    try:
        with ftplib.FTP(ftp_host) as ftp:
            ftp.login() # Anonymous login
            ftp.cwd(ftp_path)
            
            filenames = ftp.nlst()
            downloaded_files = []

            for filename in filenames:
                if filename_pattern and not filename.endswith(filename_pattern.strip("*")): # Simple pattern matching
                    continue
                if filename in [".", ".."]:
                    continue

                local_filepath = os.path.join(local_dir, filename)
                
                # Check if it's a directory (IMGT FTP can have subdirs)
                try:
                    ftp.cwd(filename) # Try to change directory
                    print(f"[IMGT/IPD] Skipping directory: {filename}")
                    ftp.cwd("..") # Go back
                    continue
                except ftplib.error_perm: # It's a file
                    print(f"[IMGT/IPD] Downloading {filename} to {local_filepath}...")
                    with open(local_filepath, 'wb') as f_local:
                        ftp.retrbinary(f'RETR {filename}', f_local.write)
                    downloaded_files.append(local_filepath)
            
            if downloaded_files:
                print(f"[IMGT/IPD] Downloaded {len(downloaded_files)} files to {local_dir}")
            else:
                print(f"[IMGT/IPD] No files matching pattern '{filename_pattern}' found or downloaded from {ftp_path}")

    except ftplib.all_errors as e:
        print(f"[IMGT/IPD] FTP Error: {e}")
    except Exception as e:
        print(f"[IMGT/IPD] General Error: {e}")


# --- 3. IEDB (Immune Epitope Database) ---
# IEDB has a REST API, useful for querying epitopes by MHC allele, or MHC binding predictions.
# Here's an example for getting MHC allele details, which might link to sequences.

def query_iedb_mhc_allele_info(mhc_allele_name):
    """
    Queries IEDB for information about a specific MHC allele.
    The IEDB API is more focused on epitopes and binding than raw sequence retrieval.
    This is an example of how to interact with its API.

    Args:
        mhc_allele_name (str): e.g., "HLA-A*02:01" (URL encoding might be needed for some tools)

    Returns:
        dict: JSON response from IEDB or None.
    """
    print(f"\n[IEDB] Querying for MHC allele: {mhc_allele_name}")
    # Example: Using the mhc_allele_details endpoint (hypothetical or simplified)
    # Actual IEDB API usage for specific tasks (like binding prediction) can be more complex.
    # For raw sequences, you'd typically get them from IMGT/NCBI and then use IEDB tools.

    # This specific endpoint might not exist or work this way; it's illustrative.
    # A more common use case would be submitting a sequence for binding prediction.
    # Let's try getting epitopes for an allele:
    # https://www.iedb.org/tools_api_v2_swagger/
    # GET /allele_details
    # Note: The IEDB API can be complex and might require specific formatting.
    # This is a simplified example. For MHC sequences, IMGT or NCBI are primary.
    
    # Let's try fetching epitopes for a given MHC allele as an example of IEDB API interaction
    # From IEDB API docs: http://tools-cluster-interface.iedb.org/tools_api/latest/mhcbinding_tutorial
    # This is for binding prediction, not direct sequence retrieval.
    # For allele specific information, often you'd go to their website or download their data.
    
    # Let's use a different IEDB API endpoint that might be more relevant for allele listing or details
    # (This is more of a placeholder as direct MHC sequence download from IEDB API isn't its primary function)
    # IEDB's strength is epitope data and prediction tools.
    
    # Example: search for epitopes restricted by a certain MHC molecule
    # This is more about epitopes than the MHC sequence itself.
    url_safe_allele = requests.utils.quote(mhc_allele_name)
    # Example: Fetch epitopes for a given MHC allele
    # This is a search query, not a direct allele sequence download.
    # The IEDB API has changed, and direct allele sequence query is not straightforward.
    # Usually, you get sequences from IMGT/NCBI and use IEDB for analysis.
    
    # Let's try a different approach: finding T cell assays for an MHC allele
    # This demonstrates API interaction but not sequence retrieval.
    # IEDB API can be complex, refer to: https://www.iedb.org/tools_api_v2_swagger/
    # For actual MHC sequences, stick to IMGT and NCBI.
    
    # Let's try a search for MHC alleles on IEDB (web scrape simulation, not ideal)
    # Or better, use their downloadable files for MHC allele lists.
    # For this script, we'll focus on NCBI and IMGT for sequences.
    # IEDB is more for *using* MHC information.

    # A more practical IEDB interaction: list MHC alleles of a certain type
    # (using a hypothetical endpoint or one that might require specific parameters)
    # For now, this function is a placeholder as IEDB's API is more for tools.
    
    print("[IEDB] Note: IEDB is primarily for epitope data and analysis tools.")
    print("[IEDB] For MHC sequences, IMGT/HLA, IPD-MHC, and NCBI are the primary sources.")
    print("[IEDB] Example: To find epitopes for an allele, you'd use their web interface or specific API tools.")
    
    # A more realistic IEDB API example: MHC-I binding prediction
    # method:ann, sequence:SLYNTVATL, allele:HLA-A*01:01, length:9
    # This is for prediction, not sequence retrieval.
    iedb_api_url = "http://tools-cluster-interface.iedb.org/tools_api/mhci/" # Example endpoint
    payload = {
        'method': 'recommended', # Or specific methods like 'ann', 'smm'
        'sequence_text': 'SYFPEITHI\nNLVPMVATV', # Example peptide(s)
        'allele': mhc_allele_name, # e.g., "HLA-A*02:01"
        # 'length': '9' # Can specify length or let it infer
    }
    headers = {'Content-Type': 'application/x-www-form-urlencoded', 'Accept': 'application/json'}
    
    try:
        print(f"[IEDB] Submitting peptide binding prediction for allele {mhc_allele_name} (example).")
        response = requests.post(iedb_api_url, data=payload, headers=headers)
        response.raise_for_status() # Raise an exception for bad status codes
        
        # The response is typically text, needs parsing.
        # For JSON: results = response.json()
        # For TSV/CSV: pd.read_csv(io.StringIO(response.text), sep='\t')
        print(f"[IEDB] Prediction API Response (first 500 chars):\n{response.text[:500]}...")
        # You would parse this text (often TSV) into a pandas DataFrame
        # from io import StringIO
        # df = pd.read_csv(StringIO(response.text), sep='\t')
        # print(df.head())
        return response.text # Or parsed data
    except requests.exceptions.RequestException as e:
        print(f"[IEDB] API request error: {e}")
        if hasattr(e, 'response') and e.response is not None:
            print(f"[IEDB] Response content: {e.response.text}")
        return None
    except Exception as e:
        print(f"[IEDB] General error with IEDB API: {e}")
        return None


# --- 4. Basic Analysis (Example: Reading FASTA and printing info) ---
def analyze_sequences_from_file(fasta_file):
    """
    Performs very basic analysis on a FASTA file (e.g., count sequences, print IDs).
    """
    if not fasta_file or not os.path.exists(fasta_file):
        print(f"[Analysis] File not found: {fasta_file}")
        return

    print(f"\n[Analysis] Analyzing sequences from: {fasta_file}")
    count = 0
    total_length = 0
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            count += 1
            total_length += len(record.seq)
            # print(f"> {record.id} (Length: {len(record.seq)})")
            # print(f"{record.seq[:60]}...") # Print first 60 bases
        print(f"[Analysis] Found {count} sequences.")
        if count > 0:
            print(f"[Analysis] Average length: {total_length / count:.2f} bp/aa")
    except Exception as e:
        print(f"[Analysis] Error parsing FASTA file {fasta_file}: {e}")

# --- Main Execution ---
if __name__ == "__main__":
    # --- NCBI Example ---
    # 1. Search for Human HLA-A nucleotide sequences
    hla_a_ids = query_ncbi_mhc(gene_name="HLA-A", organism="Homo sapiens", database="nucleotide", retmax=5)
    if hla_a_ids:
        # 2. Fetch these sequences
        hla_fasta_file = fetch_sequences_from_ncbi(hla_a_ids, database="nucleotide", output_format="fasta", filename_prefix="human_hla_a")
        if hla_fasta_file:
            analyze_sequences_from_file(hla_fasta_file)

    # 3. Search for Bovine MHC Class I protein sequences
    bola_ids = query_ncbi_mhc(gene_name="MHC class I", organism="Bos taurus", database="protein", retmax=3)
    if bola_ids:
        bola_protein_fasta = fetch_sequences_from_ncbi(bola_ids, database="protein", output_format="fasta", filename_prefix="bovine_mhc_protein")
        if bola_protein_fasta:
            analyze_sequences_from_file(bola_protein_fasta)

    # --- IMGT/IPD Example (Downloading data) ---
    # For IMGT/HLA (Human) - Nucleotide sequences
    # Check current FTP paths on https://www.ebi.ac.uk/ipd/imgt/hla/download/
    # The exact path might change, verify on the IMGT/EBI website.
    # Common files: hla_nuc.fasta (all nucleotide), hla_prot.fasta (all protein)
    # download_imgt_data(
    #     ftp_host="ftp.ebi.ac.uk",
    #     ftp_path="/pub/databases/ipd/imgt/hla/fasta/", # This path might change!
    #     local_dir="imgt_hla_data",
    #     filename_pattern="hla_nuc.fasta" # Download specific file
    # )
    # if os.path.exists("imgt_hla_data/hla_nuc.fasta"):
    #     analyze_sequences_from_file("imgt_hla_data/hla_nuc.fasta")

    # For IPD-MHC (Non-human) - Example: Dog (Canis lupus familiaris)
    # Check current FTP paths on https://www.ebi.ac.uk/ipd/mhc/download/
    # download_imgt_data(
    #     ftp_host="ftp.ebi.ac.uk",
    #     ftp_path="/pub/databases/ipd/mhc/dla/fasta/", # DLA is dog MHC, path might change!
    #     local_dir="ipd_mhc_dla_data",
    #     filename_pattern="dla_nuc.fasta" # Or specific gene files if available
    # )
    # if os.path.exists("ipd_mhc_dla_data/dla_nuc.fasta"):
    #     analyze_sequences_from_file("ipd_mhc_dla_data/dla_nuc.fasta")

    print("\n--- IMGT/IPD FTP Download Note ---")
    print("The IMGT/IPD FTP download examples are commented out.")
    print("Please uncomment and VERIFY THE FTP PATHS on the official EBI/IPD websites before running.")
    print("FTP structure and filenames can change.")
    print("Example IMGT/HLA path: /pub/databases/ipd/imgt/hla/ (then explore subdirectories like 'fasta', 'Allele_lists')")
    print("Example IPD-MHC path: /pub/databases/ipd/mhc/ (then explore by species, e.g., 'macaque', 'dla')")


    # --- IEDB Example ---
    # This example shows interaction with IEDB's prediction API, not direct sequence retrieval.
    # For sequences, use NCBI/IMGT. Then use IEDB for functional analysis.
    iedb_results = query_iedb_mhc_allele_info(mhc_allele_name="HLA-A*02:01")
    # if iedb_results:
    #     print("\n[IEDB] Further processing of IEDB results would happen here.")
    #     # For example, if it was a TSV, parse with pandas:
    #     # from io import StringIO
    #     # iedb_df = pd.read_csv(StringIO(iedb_results), sep='\t')
    #     # print(iedb_df.head())


    print("\nScript finished.")
    print("Remember to replace 'your_email@example.com' with your actual email for NCBI Entrez.")
    print("For higher NCBI request rates, consider getting an API key.")
```

## Explanation and How to Use the Script

1.  **Prerequisites:**
    *   Install Python 3.
    *   Install necessary libraries:
        ```bash
        pip install biopython requests pandas
        ```

2.  **Configuration:**
    *   `Entrez.email`: **Crucial!** In the Python script, replace `"your_email@example.com"` with your real email address. NCBI requires this to contact you if there's an issue with your script.
    *   `API_KEY`: Optional. If you make many requests to NCBI, get an API key from your NCBI account settings for higher rate limits and assign it to the `API_KEY` variable in the script.

3.  **NCBI Interaction (`query_ncbi_mhc`, `fetch_sequences_from_ncbi`):**
    *   The `query_ncbi_mhc` function searches NCBI's nucleotide or protein databases.
        *   `gene_name`: Can be specific (e.g., `"HLA-A*02:01"`) or general (e.g., `"MHC class II beta chain"`).
        *   `organism`: Scientific name (e.g., `"Homo sapiens"`, `"Mus musculus"`).
        *   `database`: Accepts `"nucleotide"` or `"protein"`.
        *   `retmax`: Defines how many initial IDs to fetch (the search might find more).
    *   The `fetch_sequences_from_ncbi` function downloads sequences for the given IDs into a FASTA (or GenBank, etc.) file.

4.  **IMGT/IPD Interaction (`download_imgt_data`):**
    *   IMGT/HLA and IPD-MHC are best accessed by downloading their curated datasets via FTP.
    *   The `download_imgt_data` function provides a basic way to connect to their FTP server (hosted by EBI) and download files.
    *   **Important:**
        *   FTP paths (like `ftp_path`) and filenames can change! **Always verify the correct paths on the official IMGT/IPD websites** before running the download parts of the script.
        *   Official Download Websites:
            *   IMGT/HLA Downloads: [https://www.ebi.ac.uk/ipd/imgt/hla/download/](https://www.ebi.ac.uk/ipd/imgt/hla/download/)
            *   IPD-MHC Downloads: [https://www.ebi.ac.uk/ipd/mhc/download/](https://www.ebi.ac.uk/ipd/mhc/download/)
        *   The script examples for IMGT/IPD are commented out by default. You'll need to uncomment them and provide correct, up-to-date paths if you wish to use them.
    *   Commonly downloaded files are FASTA files containing all nucleotide (`*_nuc.fasta`) or protein (`*_prot.fasta`) sequences for a given locus or species.

5.  **IEDB Interaction (`query_iedb_mhc_allele_info`):**
    *   IEDB's strength is in immune epitope data and MHC binding prediction tools, not primarily raw MHC sequence hosting (though it links to them).
    *   The example function `query_iedb_mhc_allele_info` demonstrates how you *might* interact with an IEDB API endpoint. The example shown uses the MHC-I binding prediction tool.
    *   For actual MHC sequences, you'd typically get them from IMGT or NCBI and then use IEDB's tools for analysis (e.g., predicting which peptides bind to a specific MHC allele sequence you've obtained).
    *   IEDB API documentation can be found at [http://tools-cluster-interface.iedb.org/tools_api/latest/](http://tools-cluster-interface.iedb.org/tools_api/latest/) or the newer [https://www.iedb.org/tools_api_v2_swagger/](https://www.iedb.org/tools_api_v2_swagger/).

6.  **Basic Analysis (`analyze_sequences_from_file`):**
    *   This is a placeholder function demonstrating how to use `Bio.SeqIO` to parse a downloaded FASTA file, count sequences, and print basic information. This is where you would add more sophisticated analysis (alignments, motif finding, phylogenetic analysis, etc.).

7.  **Main Execution Block (`if __name__ == "__main__":`)**
    *   This block in the script shows example usage of the defined functions.
    *   NCBI examples are set to run by default if the script is executed directly.
    *   IMGT/IPD download examples are commented out – you'll need to uncomment them and update paths to use them.
    *   The IEDB example demonstrates calling one of its tool APIs.

## Further Analysis (Beyond this script)

Once you have the sequences (e.g., in a FASTA file from the script):

*   **Sequence Alignment:** Use Biopython's `Bio.Align` or external tools like ClustalW, MAFFT, MUSCLE.
*   **Phylogenetic Analysis:** Tools like PhyML, RAxML, IQ-TREE (often used after alignment).
*   **Motif Discovery:** MEME suite.
*   **Variant Calling:** If you have raw sequencing reads, align them to reference MHC sequences.
*   **MHC Typing:** Specialized tools for assigning HLA/MHC alleles from sequencing data (e.g., OptiType, HLA-LA, HISAT-genotype).
*   **Structural Analysis:** If you have protein sequences, predict structures (AlphaFold) or analyze existing PDB structures.

This script provides a solid foundation for accessing MHC sequence data programmatically. Remember to consult the documentation of each database for the most up-to-date access methods and data formats.

---
