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
    (When installed via setup.py, it can be run as `query_ncbi_gi_py "<search_term>"`)
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

### 3. `blast_ncbi_tabular.py`

*   **Description:** This script submits a sequence to NCBI BLAST, polls for results, and prints them in a tabular format. It allows customization of the BLAST database, program, and other parameters. This is a Python implementation of the original Perl script `script1_blast_ncbi_tab.pl`.

*   **Usage:**
    ```bash
    python blast_ncbi_tabular.py --query "<sequence>" [options]
    ```
    (When installed via setup.py, it can be run as `blast_ncbi_tabular_py --query "<sequence>" [options]`)

*   **Command-Line Arguments:**
    *   `--query`: (Required) Input sequence string.
    *   `--database`: BLAST database. Default: `nr`.
    *   `--program`: BLAST program. Default: `blastp`.
    *   `--filter`: Filter option (e.g., 'L' for low complexity, 'F' for no filter). Default: `L`.
    *   `--hitlist_size`: Number of hits to return. Default: `20`.
    *   `--evalue`: E-value threshold. Default: `0.01`.
    *   `--poll_interval`: Seconds between polls for results. Default: `10`.
    *   `--max_poll_attempts`: Maximum poll attempts before timeout. Default: `120` (e.g., 120 * 10s = 20 minutes).
    *   `--ncbi_url`: NCBI BLAST URL. Default: `https://www.ncbi.nlm.nih.gov/blast/Blast.cgi`.

*   **Output:**
    *   Prints informational messages to `stderr` (RID, polling status).
    *   Prints a formatted header: `#	Subject ID (Processed)	Identity	Length	E-value`.
    *   Prints tabular BLAST results to `stdout`, with subject IDs processed for readability (e.g., `gb|AAA123.1|` becomes `gb:AAA123.1`).
    *   If no hits are found, a message is printed.
    *   Error messages are printed to `stderr`.

---

### 4. `query_ncbi_entrez.py`

*   **Description:** This script performs a search on NCBI Entrez databases using ESearch and then fetches the corresponding records using EFetch. It allows specification of search and fetch databases, retrieval types, and modes. This is a Python implementation of the original Perl script `script2_query_ncbi.pl`.

*   **Usage:**
    ```bash
    python query_ncbi_entrez.py --term "<search_term>" [options]
    ```
    (When installed via setup.py, it can be run as `query_ncbi_entrez_py --term "<search_term>" [options]`)

*   **Command-Line Arguments:**
    *   `--term`: (Required) Search term (e.g., 'insulin human').
    *   `--search_db`: Database for ESearch. Default: `protein`.
    *   `--fetch_db`: Database for EFetch. Default: Matches `search_db` if not specified.
    *   `--rettype`: Retrieval type for EFetch (e.g., 'gb', 'fasta', 'xml'). Default: `gb`.
    *   `--retmode`: Retrieval mode for EFetch (e.g., 'text', 'xml'). Default: `text`.
    *   `--email`: User email for NCBI E-utils (recommended). Default: `None`.
    *   `--api_key`: NCBI API key (recommended for higher request rates). Default: `None`.
    *   `--max_retries`: Max retries for HTTP requests. Default: `3`.
    *   `--base_eutils_url`: Base URL for NCBI E-utilities. Default: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/`.

*   **Output:**
    *   Prints informational messages to `stderr` (UIDs found, fetching status).
    *   Prints fetched records (e.g., GenBank, FASTA) directly to `stdout`.
    *   Error messages are printed to `stderr`.

---

## Running Tests

Unit tests are provided for all NCBI scripts and are located within this directory (e.g., `test_blast_ncbi_seq_def.py`). These tests use Python's built-in `unittest` framework and mock external NCBI calls.

Tests for NCBI scripts can be run using the project's general testing command (see the main `README.md` "Testing" section), or more specifically by targeting this directory from the project root:

```bash
python -m unittest discover pylib/scripts/ncbi/
```

This command will automatically find and execute all files named `test_*.py` within the `pylib/scripts/ncbi/` directory. This includes tests for all scripts described herein, such as `test_blast_ncbi_seq_def.py`, `test_query_ncbi_gi.py`, `test_blast_ncbi_tabular.py`, and `test_query_ncbi_entrez.py`.

Running these tests is recommended to verify script functionality, especially after any modifications.
