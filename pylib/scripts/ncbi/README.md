# NCBI Python Scripts

## Overview

These Python scripts are modern implementations of original Perl scripts designed to interact with NCBI (National Center for Biotechnology Information) services. They provide functionality for performing BLAST searches and querying the NCBI Nucleotide database (via Entrez/E-utilities) to retrieve GenBank records.

## Dependencies

The primary external Python library required is `requests` for making HTTP requests.

### Installation

You can install the necessary dependencies using pip:

```bash
pip install requests
```

## Scripts

Below is a description of each script, its usage, and expected input/output.

---

### 1. `blast_ncbi_seq_def.py`

*   **Description:** This script takes a FASTA file containing protein sequences and an E-value threshold as input. It submits each sequence to the NCBI BLASTp service (searching the 'nr' database), polls for results, and then prints the full BLAST report for each sequence once ready.

*   **Usage:**
    ```bash
    python ncbi/blast_ncbi_seq_def.py <fasta_file_path> <e_value>
    ```

*   **Input:**
    *   `<fasta_file_path>`: Path to a valid FASTA file containing one or more protein sequences.
    *   `<e_value>`: The expect value threshold for the BLAST search (e.g., `0.01`, `1e-5`).

*   **Output:**
    *   Prints the input parameters.
    *   For each sequence:
        *   Prints a message indicating which sequence is being processed.
        *   Prints polling status (attempt number and status like "WAITING", "READY").
        *   Once the BLAST search is complete, prints the full text-based BLAST report.
    *   Prints error messages if issues occur (e.g., "bad input" for unknown RID, "timed out at 20 minutes", HTTP errors).

---

### 2. `query_ncbi_gi.py`

*   **Description:** This script queries the NCBI Nucleotide database for a given search term. It retrieves a list of UIDs (Unique Identifiers) matching the term and then, for each UID, fetches and prints its corresponding GenBank record.

*   **Usage:**
    ```bash
    python ncbi/query_ncbi_gi.py "<search_term>"
    ```
    *Note: It's recommended to enclose the search term in quotes if it contains spaces or special characters.*

*   **Input:**
    *   `<search_term>`: The term to search for in the NCBI Nucleotide database (e.g., "BRCA1 human", "COVID-19").

*   **Output:**
    *   Prints the search term being used.
    *   Prints the list of UIDs found for the search term, or a message if no UIDs are found.
    *   For each UID:
        *   Prints a message indicating which UID is being fetched.
        *   Prints the full GenBank record starting from the "LOCUS" line.
    *   Prints error messages if issues occur (e.g., HTTP errors, "LOCUS line not found").
    *   Prints "Script finished." at the end, after a final 5-second delay.

---

## Running Tests (Optional but Recommended)

Unit tests are provided for both scripts to ensure their functionality.

You can run the tests using the `unittest` module from the parent directory of `ncbi/`:

*   To run all tests within the `ncbi` directory:
    ```bash
    python -m unittest discover ncbi
    ```

*   Alternatively, to run tests for specific files:
    ```bash
    python -m unittest ncbi/test_blast_ncbi_seq_def.py ncbi/test_query_ncbi_gi.py
    ```

The tests mock external NCBI calls and verify various aspects of the scripts, including argument parsing, API interaction logic, output formatting, and error handling.
