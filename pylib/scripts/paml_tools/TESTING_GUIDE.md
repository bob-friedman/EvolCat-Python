# Testing Guide for calculate_dn_ds.py

This guide explains how to run a simulated test for the `calculate_dn_ds.py` script using mock data. This allows testing the script's parsing and output formatting logic without requiring a PAML installation.

## Prerequisites

1.  Ensure the following files are in the same directory:
    *   `calculate_dn_ds.py` (the script to be tested)
    *   `sample_aligned.fasta` (the sample input FASTA file)
    *   `mock_yn00_output.txt` (the mock `yn00` output data)
    *   `expected_script_output.txt` (the expected output from the script for this test case)

2.  Python 3 must be installed to run the script.
3.  The script relies on the BioPython library (`Bio.SeqIO`, `Bio.Phylo.PAML`). If not installed, then may need to install it:
    ```bash
    pip install biopython
    ```

## Running the Simulated Test

Follow these steps to execute the test:

1.  **Open terminal or command prompt.**

2.  **Navigate to the directory** containing the script and test files.

3.  **Set Environment Variables:**
    *   To activate the simulation mode in the script:
        ```bash
        export SIMULATE_PAML=true
        ```
        (On Windows, use `set SIMULATE_PAML=true`)
    *   To tell the script to use mock `yn00` output file:
        ```bash
        export MOCK_YN00_FILE_PATH=mock_yn00_output.txt
        ```
        (On Windows, use `set MOCK_YN00_FILE_PATH=mock_yn00_output.txt`)

4.  **Run the script:**
    Execute `calculate_dn_ds.py` with `sample_aligned.fasta` as its argument:
    ```bash
    python calculate_dn_ds.py sample_aligned.fasta
    ```

5.  **Verify the Output:**
    The script will print its results to standard output. Compare this output with the content of `expected_script_output.txt`.

    Also, one can redirect the script's output to a file and then use a diff tool:
    ```bash
    python calculate_dn_ds.py sample_aligned.fasta > actual_output.txt
    ```
    Then compare `actual_output.txt` with `expected_script_output.txt`. On Linux or macOS, then may use `diff`:
    ```bash
    diff actual_output.txt expected_script_output.txt
    ```
    If there are no differences, the `diff` command will produce no output, indicating the test passed. If there are differences, `diff` will show them.

    On Windows, then may try `FC`:
    ```bash
    FC actual_output.txt expected_script_output.txt
    ```

    The output lines from the script starting with "INFO:", "WARNING:", or "ERROR:" are for logging and should not be part of the `actual_output.txt` if only comparing the dN/dS results. The `expected_script_output.txt` only contains the tab-separated data lines. May need to manually filter the script's output if it includes these informational lines intermingled with the results, or adjust the comparison. However, the current script structure prints informational lines separately from the final results block.

    The `expected_script_output.txt` should look like this:
    ```
    Seq1_ID	Seq2_ID	dN	dS	dN_dS_ratio
    Seq1	Seq2	0.0123	0.0456	0.200
    Seq1	Seq3	0.0010	0.0100	0.100
    Seq1	Seq4	0.0015	0.0100	0.150
    Seq2	Seq3	0.0300	0.0600	0.300
    Seq2	Seq4	0.0350	0.0700	0.350
    Seq3	Seq4	0.0005	0.0100	0.050
    ```
    The script first prints some "INFO:" lines, then a header "Results (Seq1_ID...", then the data. For comparison, one should typically compare from the header line onwards.

## Important Notes
*   The `MOCK_YN00_FILE_PATH` environment variable makes the script use the content of `mock_yn00_output.txt` as if it were the output from PAML's `yn00`.
*   If `SIMULATE_PAML` is set to `false`, the script will attempt to run the actual `yn00` program from PAML, which must be installed and in the system's PATH. This testing guide focuses only on the simulated run.

## Note on PAML's `codeml` Program

The `calculate_dn_ds.py` script is designed to use the `yn00` program from the PAML package. `yn00` is specifically used for calculating pairwise dN/dS rates between sequences.

PAML also includes a more powerful and versatile program called `codeml`. `codeml` can perform more complex dN/dS analyses, including:
*   Detecting positive selection at specific codon sites within an alignment.
*   Testing different evolutionary models.
*   Analyzing dN/dS rates along different branches of a phylogenetic tree.

If research requires these more advanced types of dN/dS analysis, then may explore using PAML's `codeml` program directly or through other tools that might wrap it. The `calculate_dn_ds.py` script in its current form is focused on the straightforward pairwise estimation provided by `yn00`.
