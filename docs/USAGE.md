# Script Usage Guide

This document provides usage instructions for the converted Python scripts.

## General Notes

- All scripts can be run from the command line using `python3 path/to/script.py [arguments]`.
- Ensure you have Python 3 and the required libraries (Biopython, Matplotlib) installed.
- Use the `-h` or `--help` flag with any script to see its specific command-line options.

---

## `pylib/scripts/gb2fasta.py`

Converts a GenBank file to FASTA format.

**Usage:**
```bash
python3 pylib/scripts/gb2fasta.py <input_genbank_file>
```
- `<input_genbank_file>`: Path to the input file in GenBank format.
- FASTA output is printed to standard output.

**Example:**
```bash
python3 pylib/scripts/gb2fasta.py my_sequences.gb > my_sequences.fasta
```

---

## `pylib/scripts/translate_seq.py`

Translates nucleotide sequences from a FASTA file to protein sequences.

**Usage:**
```bash
python3 pylib/scripts/translate_seq.py --input_file <fasta_file> [--frame <1|2|3>] [--table <ncbi_table_id>]
```
- `--input_file <fasta_file>`: Path to the input FASTA file (required).
- `--frame <1|2|3>`: Translation frame (default: 1).
- `--table <ncbi_table_id>`: NCBI translation table ID (default: 1, Standard Code).
- Translated protein sequences are printed to standard output, with headers indicating the original sequence ID and translation parameters.

**Example:**
```bash
python3 pylib/scripts/translate_seq.py --input_file cds.fna --frame 1 --table 11 > proteins.faa
```

---

## `pylib/scripts/join_tables.py`

Joins two tab-delimited files based on the first column. This mimics a left join where the first file is the left table.

**Usage:**
```bash
python3 pylib/scripts/join_tables.py <file1_path> <file2_path>
```
- `<file1_path>`: Path to the first tab-delimited file (left table).
- `<file2_path>`: Path to the second tab-delimited file (right table).
- Joined output is printed to standard output.

**Example:**
```bash
python3 pylib/scripts/join_tables.py annotations.tsv metadata.tsv > combined_data.tsv
```

---

## `pylib/scripts/dot_plot.py`

Generates a dot plot from two FASTA sequences and saves it as a PNG image. It first calculates matches and stores them in an intermediate "dot file".

**Usage:**
```bash
python3 pylib/scripts/dot_plot.py -1 <seqfile1> -2 <seqfile2> -w <wordlen> -s <step> -o <outfile.png> -d <dotfile.txt> [-t <title>]
```
- `-1, --seqfile1 <seqfile1>`: Path to the first input FASTA file (required).
- `-2, --seqfile2 <seqfile2>`: Path to the second input FASTA file (required).
- `-w, --wordlen <wordlen>`: Word size for matching (integer, required).
- `-s, --step <step>`: Step length for moving the word (integer, required).
- `-o, --outfile <outfile.png>`: Output PNG file name (required).
- `-d, --dotfile <dotfile.txt>`: Path to the intermediate file to store match coordinates (required).
- `-t, --title <title>`: Title for the plot (string, default: empty string).

**Example:**
```bash
python3 pylib/scripts/dot_plot.py -1 seqA.fasta -2 seqB.fasta -w 10 -s 5 -o comparison.png -d matches.txt -t "SeqA vs SeqB"
```
---

## `pylib/scripts/fas2csv.py`

Converts a FASTA file to a CSV format, where each line is `name,sequence`.

**Usage:**
```bash
python3 pylib/scripts/fas2csv.py <input_fasta_file>
```
- `<input_fasta_file>`: Path to the input FASTA file.
- CSV output is printed to standard output.

**Example:**
```bash
python3 pylib/scripts/fas2csv.py my_sequences.fasta > my_sequences.csv
```
---

## `pylib/scripts/fas2meg.py`

Converts a FASTA file to MEGA (.meg) format.

**Usage:**
```bash
python3 pylib/scripts/fas2meg.py <input_fasta_file>
```
- `<input_fasta_file>`: Path to the input FASTA file.
- MEGA formatted output is printed to standard output.

**Example:**
```bash
python3 pylib/scripts/fas2meg.py my_sequences.fasta > my_sequences.meg
```
---

## `pylib/scripts/fas2phy.py`

Converts a FASTA file to PHYLIP (sequential) format.
**Important:** All sequences in the input FASTA file must be of the same length for PHYLIP format.

**Usage:**
```bash
python3 pylib/scripts/fas2phy.py <input_fasta_file>
```
- `<input_fasta_file>`: Path to the input FASTA file.
- PHYLIP formatted output is printed to standard output.

**Example:**
```bash
python3 pylib/scripts/fas2phy.py my_aligned_sequences.fasta > my_sequences.phy
```
---

## `pylib/scripts/phy2fas.py`

Converts a PHYLIP file (interleaved or sequential) to FASTA format.

**Usage:**
```bash
python3 pylib/scripts/phy2fas.py <input_phylip_file>
```
- `<input_phylip_file>`: Path to the input PHYLIP file.
- FASTA formatted output is printed to standard output.

**Example:**
```bash
python3 pylib/scripts/phy2fas.py my_sequences.phy > my_sequences.fasta
```
---

## `pylib/scripts/phy2meg.py`

Converts a PHYLIP file (interleaved or sequential) to MEGA (.meg) format.

**Usage:**
```bash
python3 pylib/scripts/phy2meg.py <input_phylip_file>
```
- `<input_phylip_file>`: Path to the input PHYLIP file.
- MEGA formatted output is printed to standard output.

**Example:**
```bash
python3 pylib/scripts/phy2meg.py my_sequences.phy > my_sequences.meg
```
---

## `pylib/scripts/gbCDS.py`

Extracts information (Locus ID, nucleotide sequence, translation) for the first CDS feature in each record of a GenBank file.

**Usage:**
```bash
python3 pylib/scripts/gbCDS.py <input_genbank_file>
```
- `<input_genbank_file>`: Path to the input GenBank file.
- Output is printed to standard output, with Locus ID, CDS nucleotide sequence, and protein translation each on a new line for every record. If a record has no CDS, it prints Locus ID, "No CDS feature found", and "N/A".

**Example:**
```bash
python3 pylib/scripts/gbCDS.py my_genome.gb
```
---

## `pylib/scripts/clean_fasta_name.py`

Cleans FASTA headers from an input file and prints the modified FASTA to standard output.
The cleaning process involves:
1.  Replacing all occurrences of one or more whitespace characters with a single underscore (`_`).
2.  Replacing all period characters (`.`) with underscores (`_`).
3.  Converting the entire header string to uppercase.

**Usage:**
```bash
python3 pylib/scripts/clean_fasta_name.py <input_fasta_file>
```
- `<input_fasta_file>`: Path to the input FASTA file.
- Modified FASTA content is printed to standard output.

**Example:**
```bash
python3 pylib/scripts/clean_fasta_name.py original.fasta > cleaned.fasta
```
---

## `pylib/scripts/extract_region.py`

Extracts a specified region from sequences in a FASTA file.

**Usage:**
```bash
python3 pylib/scripts/extract_region.py <input_fasta_file> <start_bp> <end_bp> [--output_file <output.fasta>]
```
- `<input_fasta_file>`: Path to the input FASTA file (required).
- `<start_bp>`: Start base pair of the region to extract (1-based, integer, required).
- `<end_bp>`: End base pair of the region to extract (1-based, integer, required).
- `--output_file <output.fasta>`: Optional path to the output FASTA file. If not provided, output is to standard output.

**Details:**
- Coordinates are 1-based and inclusive.
- If the specified region is out of bounds for a sequence, a warning is printed to `stderr` for that sequence, and it's skipped.
- The output FASTA headers will be modified to include the region extracted (e.g., `original_id_region_start-end`).

**Example:**
```bash
python3 pylib/scripts/extract_region.py my_sequences.fasta 100 200 --output_file extracted_regions.fasta
# To print to stdout:
python3 pylib/scripts/extract_region.py my_sequences.fasta 50 150
```
---

## `pylib/scripts/merge_fastas.py`

Merges multiple FASTA files. Sequences with identical headers (full line after '>') are concatenated.
The order of sequences in the output is based on the first appearance of each header across all input files.

**Usage:**
```bash
python3 pylib/scripts/merge_fastas.py <input_fasta_file1> [<input_fasta_file2> ...]
```
- `<input_fasta_file1>`: Path to the first input FASTA file (required).
- `[<input_fasta_file2> ...]`: Paths to additional input FASTA files (optional).
- Merged FASTA content is printed to standard output.

**Example:**
```bash
python3 pylib/scripts/merge_fastas.py file1.fasta file2.fasta file3.fasta > merged_output.fasta
```
---

## `pylib/scripts/rev_comp.py`

Generates the reverse complement of sequences from a FASTA file and prints the result to standard output.

**Usage:**
```bash
python3 pylib/scripts/rev_comp.py <input_fasta_file>
```
- `<input_fasta_file>`: Path to the input FASTA file.
- Reverse complemented sequences are printed to standard output in FASTA format.
- Headers are modified: `original_id` becomes `original_id_rc`, and `(reverse complement)` is appended to the description.

**Example:**
```bash
python3 pylib/scripts/rev_comp.py my_sequences.fasta > my_sequences_rc.fasta
```
---

## `pylib/scripts/nogaps.py`

Removes columns containing non-alphabetic characters (e.g., gaps `-`, `?`) from an aligned FASTA file.
All input sequences must be of the same length. Output is to standard output.

**Usage:**
```bash
python3 pylib/scripts/nogaps.py <input_fasta_file>
```
- `<input_fasta_file>`: Path to the input aligned FASTA file.
- FASTA output with gap columns removed is printed to standard output.

**Example:**
```bash
python3 pylib/scripts/nogaps.py my_aligned_sequences.fasta > my_aligned_sequences_nogaps.fasta
```
---

## `pylib/scripts/print_duplicate_column_values.py`

Reads tab-delimited data from input file(s) (or standard input) and prints values from the first two columns that appear more than once.
The counts are cumulative across both columns. If a value appears once in column 1 and then again in column 2, its second appearance will cause it to be printed.

**Usage:**
```bash
python3 pylib/scripts/print_duplicate_column_values.py [<input_file1> <input_file2> ...]
```
- `[<input_file1> <input_file2> ...]`: Optional paths to input tab-delimited files. If not provided, reads from standard input.
- Duplicate values (as defined above) are printed to standard output, each on a new line.

**Example:**
Given a file `data.tsv`:
```
A	B
C	A
D	E
B	F
```
Command: `python3 pylib/scripts/print_duplicate_column_values.py data.tsv`
Output:
```
A
B
```
(Explanation: 'A' is printed when encountered the second time (in col2). 'B' is printed when encountered the second time (in col1 of the last row)).

---

## `pylib/scripts/sort_numbers.py`

Reads numbers (one per line) from input file(s) (or standard input), sorts them numerically, and prints them to standard output.
Non-numeric lines are skipped with a warning message printed to `stderr`.

**Usage:**
```bash
python3 pylib/scripts/sort_numbers.py [<input_file1> <input_file2> ...]
```
- `[<input_file1> <input_file2> ...]`: Optional paths to input files containing numbers, one per line. If not provided, reads from standard input.
- Sorted numbers are printed to standard output, each on a new line.

**Example:**
Given a file `numbers.txt`:
```
10
2.5
-3
1
```
Command: `python3 pylib/scripts/sort_numbers.py numbers.txt`
Output:
```
-3
1
2.5
10
```
---

## `pylib/scripts/transpose_text_matrix.py`

Transposes a character matrix where each input line is a row and each character in the line is an element.
Reads from input file(s) or standard input if no files are specified.
All input lines must have the same length.

**Usage:**
```bash
python3 pylib/scripts/transpose_text_matrix.py [<input_file1> <input_file2> ...]
```
- `[<input_file1> <input_file2> ...]`: Optional paths to input files. If not provided, reads from standard input.
- The transposed character matrix is printed to standard output. Characters in each output row are joined by a single space.

**Example:**
Given a file `text_matrix.txt`:
```
ABC
DEF
```
Command: `python3 pylib/scripts/transpose_text_matrix.py text_matrix.txt`
Output:
```
A D
B E
C F
```
---

## `pylib/scripts/iupac_to_regexp.py`

Converts IUPAC ambiguity strings (one per line from an input file) to their equivalent regular expressions.
Removes `^` characters from input strings before conversion. Unknown IUPAC codes are treated as literals and a warning is printed to `stderr`.

**Usage:**
```bash
python3 pylib/scripts/iupac_to_regexp.py <input_file>
```
- `<input_file>`: Path to the input file containing IUPAC strings, one per line.
- Regular expressions are printed to standard output, one for each processed input line.

**Example:**
Given an input file `iupac_codes.txt`:
```
RYM
N^K
XG
```
Command: `python3 pylib/scripts/iupac_to_regexp.py iupac_codes.txt`
Output:
```
[GA][CT][AC]
[ACGT][GT]
XG
```
(A warning about 'X' would also be printed to stderr).

---

## `pylib/scripts/transpose_tsv.py`

Transposes a tab-delimited matrix. Reads from input file(s) or standard input if no files are specified.
All rows in the input must have a consistent number of columns.

**Usage:**
```bash
python3 pylib/scripts/transpose_tsv.py [<input_file1> <input_file2> ...]
```
- `[<input_file1> <input_file2> ...]`: Optional paths to input tab-delimited files. If not provided, reads from standard input.
- The transposed matrix is printed to standard output, with elements in each row joined by tabs.

**Example:**
Given a file `matrix.tsv`:
```
A1	B1	C1
A2	B2	C2
```
Command: `python3 pylib/scripts/transpose_tsv.py matrix.tsv`
Output:
```
A1	A2
B1	B2
C1	C2
```
---

## `pylib/scripts/table_to_binary.py`

Converts tab-delimited table data to a binary format. Each element in the table is converted to "1" if its numeric value is greater than or equal to 1.0, and "0" otherwise. Non-numeric elements are treated as 0.
The script reads from input file(s) or standard input if no files are specified.

**Usage:**
```bash
python3 pylib/scripts/table_to_binary.py [<input_file1> <input_file2> ...]
```
- `[<input_file1> <input_file2> ...]`: Optional paths to input tab-delimited files. If not provided, reads from standard input.
- Binary-converted rows are printed to standard output, with elements joined by tabs.

**Example:**
Given a file `data.tsv`:
```
Name	Value1	Value2	Value3
ItemA	0.5	1.0	3
ItemB	text	0	-2
ItemC	2	5	0.9
```
Command: `python3 pylib/scripts/table_to_binary.py data.tsv`
Output:
```
0	1	1	1
0	0	0	0
1	1	0	0
```
---

## `pylib/scripts/count_kmers.py`

Counts k-mer frequencies in sequences from a FASTA file.
Only sequences containing valid DNA characters (A, T, C, G, case-insensitive) are processed.
Other sequences are skipped with a warning.

**Usage:**
```bash
python3 pylib/scripts/count_kmers.py <input_fasta_file> [--kmer_len <length>]
```
- `<input_fasta_file>`: Path to the input FASTA file (required).
- `--kmer_len <length>`: Length of the k-mers to count (integer, default: 5).
- Output is printed to standard output, with each line in the format `kmer:count`, sorted by count (descending) and then k-mer (alphabetically).

**Example:**
Given a file `sequences.fasta`:
```
>seq1
ATGCATGC
>seq2
GATTACA
>seq3
ATGN AT
```
Command: `python3 pylib/scripts/count_kmers.py sequences.fasta --kmer_len 3`
Output (order might vary for ties in count, but ATG should be before TGC if counts are same):
```
ATG:2
TGC:2
GCA:1
CAT:1
GAT:1
ATT:1
TTA:1
TAC:1
ACA:1
```
(A warning for `seq3` would also be printed to stderr).

---

---

## `pylib/scripts/approximate_string_match.py`

Performs approximate string matching using edit distance to find the best match(es) of a pattern within a text.

**Usage:**
```bash
python3 pylib/scripts/approximate_string_match.py --pattern <pattern_string> --text <text_string> [--print_matrix]
```
- `--pattern <pattern_string>`: The pattern string (required).
- `--text <text_string>`: The text string in which to search (required).
- `--print_matrix`: Optional flag to print the full edit distance matrix.

**Example:**
```bash
python3 pylib/scripts/approximate_string_match.py --pattern "ACTG" --text "ACGTXACTG" 
# To also print the distance matrix:
# python3 pylib/scripts/approximate_string_match.py --pattern "ACTG" --text "ACGTXACTG" --print_matrix
```
---
---

## `pylib/scripts/parse_blast_text.py`

Parses a standard text BLAST output file and prints structured information for each hit and HSP.

**Usage:**
```bash
python3 pylib/scripts/parse_blast_text.py <input_blast_file>
```
- `<input_blast_file>`: Path to the input BLAST output file (text format).

**Output Format:**
For each query, then for each subject hit:
```
Subject_ID:
     Subject Description line
     Length = Subject Length
              [E-value, Identities_Num, Alignment_Span, Gaps_Num, Orientation, Q_Start, Q_End, S_Start, S_End]
              ... (more HSPs if present)
```
(Coordinates are 1-based, inclusive.)

**Example:**
```bash
python3 pylib/scripts/parse_blast_text.py my_blast_run.txt
```
---
---

## `pylib/scripts/blast_to_table.py`

Parses a standard text BLAST output file and outputs selected fields as a tab-delimited table, excluding self-hits.

**Usage:**
```bash
python3 pylib/scripts/blast_to_table.py <input_blast_file>
```
- `<input_blast_file>`: Path to the input BLAST output file (text format).

**Output Columns (Tab-delimited):**
`Query	Subject_Locus	Subject_Accession	Subject_Length	Subject_Description	P_value_Mantissa	P_value_Exponent	Percent_Identities	Percent_Positives	Q_Start	Q_End	S_Start	S_End`
(Coordinates are 1-based, inclusive. Subject_Accession is 'N/A' in current version.)

**Example:**
```bash
python3 pylib/scripts/blast_to_table.py my_blast_run.txt > blast_summary.tsv
```
---
---

## `pylib/scripts/find_blast_top_pairs.py`

Processes a tab-delimited BLAST-like input file to find and print top-scoring unique pairs.
A pair is considered unique if neither its query ID nor its subject ID has appeared in a previously reported top pair. Self-hits (query ID equals subject ID) are excluded.

**Usage:**
```bash
python3 pylib/scripts/find_blast_top_pairs.py [input_file] [options]
```
- `[input_file]`: Path to the input tab-delimited file. If omitted, reads from standard input.

**Options:**
- `--score_col_idx INT`: Index of the score column (0-based, default: 2). Higher score is better.
- `--filter_col_idx INT`: Index of the column used for filtering (0-based, default: 3).
- `--filter_threshold FLOAT`: Minimum value for the filter column to consider a hit (default: 50.0).
- `--query_id_col_idx INT`: Index of the Query ID column (0-based, default: 0).
- `--subject_id_col_idx INT`: Index of the Subject ID column (0-based, default: 1).

**Output Format:**
Prints the first four fields of the original input lines that meet the criteria, tab-delimited.

**Example:**
Assuming `blast_results.tsv` with QueryID, SubjectID, BitScore, AlignmentLength in columns 0, 1, 2, 3:
```bash
python3 pylib/scripts/find_blast_top_pairs.py blast_results.tsv --filter_threshold 100
```
---
---

## `pylib/scripts/calculate_k2p.py`

Calculates Kimura 2-Parameter (K2P) distance, its standard error (SE), and the transition/transversion (Ts/Tv) ratio between all pairs of DNA sequences in a FASTA file. Assumes sequences are aligned and of equal length.

**Usage:**
```bash
python3 pylib/scripts/calculate_k2p.py <input_fasta_file>
```
- `<input_fasta_file>`: Path to the input FASTA file containing aligned DNA sequences.

**Output Format (Tab-delimited):**
`Gene1   Gene2   K2P_Distance    SE_K2P  Ts_Tv_Ratio`
(Calculated values are "N/A" if conditions for formulas are not met.)

**Example:**
```bash
python3 pylib/scripts/calculate_k2p.py aligned_sequences.fasta > k2p_results.tsv
```
---
---

## `pylib/scripts/calculate_dna_distances.py`

Calculates various DNA distance measures between all pairs of sequences in a FASTA file. Assumes sequences are aligned and of equal length. Sites with non-DNA characters (not A, T, C, G) are skipped.

**Usage:**
```bash
python3 pylib/scripts/calculate_dna_distances.py <input_fasta_file>
```
- `<input_fasta_file>`: Path to the input FASTA file containing aligned DNA sequences.

**Output Columns (Tab-delimited):**
`Gene1   Gene2   Ts/Tv_Ratio JC_Distance SE_JC K2P_Distance    SE_K2P  TamuraNei_Distance  SE_TamuraNei    TajimaNei_Distance  SE_TajimaNei`
(Calculated values are "N/A" if conditions for formulas are not met. SE for Tamura-Nei is "N/A" in current version.)

**Example:**
```bash
python3 pylib/scripts/calculate_dna_distances.py aligned_sequences.fasta > dna_distances_summary.tsv
```
---
---

## `pylib/scripts/find_tandem_repeats.py`

Searches for tandem repeats of a specified pattern within a given DNA sequence, allowing for a "fudge factor" (maximum allowed separation between consecutive repeats).

**Usage:**
```bash
python3 pylib/scripts/find_tandem_repeats.py --pattern <pattern_string> [--sequence <dna_string> | --sequence_file <file_path>] [--fudge <int>] [--verbose]
```
- `--pattern <pattern_string>`: The repeat pattern string (required).
- `--sequence <dna_string>`: The DNA sequence string. Required if `--sequence_file` is not used.
- `--sequence_file <file_path>`: Path to a file containing the DNA sequence. If FASTA, the first record's sequence is used. If plain text, the first non-empty line is used. Overrides `--sequence` if both provided.
- `--fudge <int>`: Maximum separation allowed between repeats to be considered tandem (default: 1).
- `--verbose`: If set, prints detailed messages during repeat discovery.

**Output Format:**
Prints a summary of found repeat blocks:
```
<N> repeat<s>
    <M> instance<s> found
        Offset <start> - <end> (<length>)
        ...
```
(Offsets are 0-based.)

**Example:**
```bash
python3 pylib/scripts/find_tandem_repeats.py --pattern "CAT" --sequence "CATCATCATNNNCATCAT" --fudge 3
```
---
---

## `pylib/scripts/calculate_k2p.py`

Calculates Kimura 2-Parameter (K2P) distance, its standard error (SE), and the transition/transversion (Ts/Tv) ratio between all pairs of DNA sequences in a FASTA file. Assumes sequences are aligned and of equal length.

**Usage:**
```bash
python3 pylib/scripts/calculate_k2p.py <input_fasta_file>
```
- `<input_fasta_file>`: Path to the input FASTA file containing aligned DNA sequences.

**Output Format (Tab-delimited):**
`Gene1   Gene2   K2P_Distance    SE_K2P  Ts_Tv_Ratio`
(Calculated values are "N/A" if conditions for formulas are not met, e.g., log of non-positive.)

**Example:**
```bash
python3 pylib/scripts/calculate_k2p.py aligned_sequences.fasta > k2p_results.tsv
```
---
---

## `pylib/scripts/calculate_dna_distances.py`

Calculates various DNA distance measures between all pairs of sequences in a FASTA file. Assumes sequences are aligned and of equal length. Sites with non-DNA characters (not A, T, C, G) are skipped.

**Usage:**
```bash
python3 pylib/scripts/calculate_dna_distances.py <input_fasta_file>
```
- `<input_fasta_file>`: Path to the input FASTA file containing aligned DNA sequences.

**Output Columns (Tab-delimited):**
`Gene1   Gene2   Ts/Tv_Ratio JC_Distance SE_JC K2P_Distance    SE_K2P  TamuraNei_Distance  SE_TamuraNei    TajimaNei_Distance  SE_TajimaNei`
(Calculated values are "N/A" if conditions for formulas are not met. SE for Tamura-Nei is currently "N/A".)

**Example:**
```bash
python3 pylib/scripts/calculate_dna_distances.py aligned_sequences.fasta > dna_distances_summary.tsv
```
---

## `pylib/scripts/calculate_nucleotide_diversity.py`

Calculates nucleotide diversity (pi) from a FASTA alignment file.

Nucleotide diversity (pi) is defined as "the average number of nucleotide differences per site between any two DNA sequences chosen randomly from the sample population".

**Dependencies:** BioPython

**Usage:**
```bash
python3 pylib/scripts/calculate_nucleotide_diversity.py <fasta_file>
```
- `<fasta_file>`: Path to the input FASTA alignment file. The sequences within this file must be aligned (i.e., all sequences must have the same length).

**Output:**
The script prints the following information to standard output:
- Number of sequences
- Alignment length
- Total pairwise differences
- Number of pairs
- Nucleotide diversity (pi) (formatted to 6 decimal places)

**Example Output:**
```
Number of sequences: 3
Alignment length: 100
Total pairwise differences: 15
Number of pairs: 3
Nucleotide diversity (pi): 0.050000
```

## `pylib/scripts/paml_tools/calculate_dn_ds.py`

Calculates pairwise dN/dS (non-synonymous to synonymous substitution rates) ratios for sequences from an aligned FASTA file. It uses the `yn00` program from the PAML package.

**PAML Dependency:**
- Requires PAML (Phylogenetic Analysis by Maximum Likelihood) to be installed, and the `yn00` executable from PAML must be in the system's PATH for actual calculations (i.e., not in simulation mode).
- PAML can be obtained from: http://abacus.gene.ucl.ac.uk/software/paml.html

**Usage:**
```bash
python pylib/scripts/paml_tools/calculate_dn_ds.py <aligned_fasta_file>
```

**Input Argument:**
- `<aligned_fasta_file>`: Path to the input FASTA file. This file must contain two or more coding sequences that have already been aligned.

**Output Format:**
Prints tab-separated values to standard output with the following columns:
`Seq1_ID	Seq2_ID	dN	dS	dN_dS_ratio`
(where `dN_dS_ratio` is the 'omega' value from `yn00`)

**Simulation Mode:**
The script can run in a simulation mode (default) for testing without requiring a PAML installation. Set the environment variable `SIMULATE_PAML=false` to run actual PAML `yn00`.
If `SIMULATE_PAML=true` (or is not set), you can also set the `MOCK_YN00_FILE_PATH` environment variable to a file mimicking `yn00` output for the simulation.

**Guidance Note:**
Note: This script is designed for pairwise dN/dS estimation using `yn00`. For more complex analyses, such as detecting positive selection at specific codon sites or using different evolutionary models along a phylogeny, consider using the `pylib/scripts/paml_tools/calculate_site_specific_ds_dn.py` script (which uses PAML's `codeml` program).

---
=======
**Example:**
See TESTING_GUIDE in the directory containing this script.

---

## `pylib/scripts/calculate_nucleotide_diversity.py`

Calculates nucleotide diversity (pi) from a FASTA alignment file.

Nucleotide diversity (pi) is defined as "the average number of nucleotide differences per site between any two DNA sequences chosen randomly from the sample population".

**Dependencies:** BioPython

**Usage:**
```bash
python3 pylib/scripts/calculate_nucleotide_diversity.py <fasta_file>
```
- `<fasta_file>`: Path to the input FASTA alignment file. The sequences within this file must be aligned (i.e., all sequences must have the same length).

**Output:**
The script prints the following information to standard output:
- Number of sequences
- Alignment length
- Total pairwise differences
- Number of pairs
- Nucleotide diversity (pi) (formatted to 6 decimal places)

**Example Output:**
```
Number of sequences: 3
Alignment length: 100
Total pairwise differences: 15
Number of pairs: 3
Nucleotide diversity (pi): 0.050000
```

## `pylib/scripts/paml_tools/calculate_site_specific_ds_dn.py`

This script wraps PAML's `codeml` program to perform site-specific dN/dS (ratio of non-synonymous to synonymous substitution rates) analysis on coding sequence alignments. It helps in identifying amino acid sites that may be under positive selection.

**Prerequisites:**
- PAML (specifically the `codeml` executable) must be installed and accessible in the system's PATH, or its location must be provided using the `--paml_path` argument.
- Python dependency: BioPython.

**Command-Line Usage:**
```bash
python pylib/scripts/paml_tools/calculate_site_specific_ds_dn.py \
       --alignment <path_to_fasta_alignment> \
       --tree <path_to_newick_tree> \
       --model <paml_model> \
       --outfile_prefix <prefix_for_outputs> \
       [options]
```

**Input Arguments:**
- `--alignment <path_to_fasta_alignment>` (required): Path to the input coding sequence alignment file in FASTA format.
- `--tree <path_to_newick_tree>` (required): Path to the phylogenetic tree file in Newick format. (Must be compatible with the alignment).
- `--model <paml_model>` (required): The PAML model. Supported models:
    - `M0` (one-ratio): Assumes a single dN/dS ratio across all sites.
    - `M1a` (neutral): Assumes two classes of sites: conserved (0 < ω0 < 1) and strictly neutral (ω1 = 1).
    - `M2a` (selection): Extends M1a by adding a third class of sites allowing for positive selection (ω2 > 1).
    - `M3` (discrete): Allows dN/dS to vary among sites, with k discrete categories of ratios.
    - `M7` (beta): Assumes a beta distribution for dN/dS ratios between 0 and 1.
    - `M8` (beta&w>1): Extends M7 by adding an extra site class that allows for positive selection (ωs > 1).
- `--outfile_prefix <prefix_for_outputs>` (required): Prefix for naming output files (e.g., `my_analysis` results in `my_analysis.mlc`, `my_analysis_site_analysis.tsv`).

**Optional Arguments:**
- `--paml_path <path_to_codeml>`: Specify the path to the `codeml` executable.
- `--verbose`: Enable detailed PAML `codeml` output.
- `--cleandata <0_or_1>`: PAML `cleandata` option (default: 1, removes alignment sites with gaps/ambiguities before analysis).
- `--num_site_categories <integer>`: Number of site categories for the M3 (discrete) model (default: 3).

**Output:**
- Main output files:
    - `<outfile_prefix>.mlc`: The primary, detailed output file from PAML `codeml`.
    - `<outfile_prefix>_site_analysis.tsv`: A tab-separated summary of site-specific dN/dS ratios. For models like M2a and M8, this includes Bayes Empirical Bayes (BEB) posterior probabilities for sites identified as potentially under positive selection. Example columns: "Site", "AminoAcid", "dN_dS", "PosteriorProbability_PositiveSelection", "Note".
- Standard Output: The script prints a summary of model parameters (e.g., log-likelihood (lnL), kappa, omega values) to the console.
- Other Files: PAML may generate additional files (e.g., `rst`, `rub`, `2NG.dN`, `2NG.dS`) in the working directory.

**Important Considerations:**
- Interpreting dN/dS ratios: dN/dS > 1 suggests positive (Darwinian) selection; dN/dS = 1 suggests neutral evolution; dN/dS < 1 suggests purifying (negative) selection.
- The choice of PAML model is critical and should be guided by the specific biological questions and hypotheses. Model comparison using Likelihood Ratio Tests (LRTs) is common practice (e.g., comparing M1a vs. M2a, or M7 vs. M8) but is not performed by this script directly.
- Note: This script is designed for site-specific evolutionary analyses using PAML's `codeml` program. For straightforward pairwise dN/dS calculations between sequences, the `pylib/scripts/paml_tools/calculate_dn_ds.py` script (which uses PAML's `yn00`) is more direct.
---
---

## `pylib/scripts/build_tree_from_distances.py`

Builds a phylogenetic tree from a distance matrix file using Neighbor-Joining (NJ) or UPGMA methods.

**Input:**
- A distance matrix file. The script expects a PHYLIP-style distance matrix (square or lower-triangular). This can be a standalone file containing just the matrix, or a matrix embedded within a simple Nexus `DISTANCES` block.
  Example of a standalone PHYLIP lower-triangular matrix (`my_distances.phy`):
  ```
     4
  Alpha
  Beta      0.20
  Gamma     0.50  0.40
  Delta     0.80  0.70  0.60
  ```
  Or a PHYLIP square matrix:
  ```
     4
  Alpha     0.00  0.20  0.50  0.80
  Beta      0.20  0.00  0.40  0.70
  Gamma     0.50  0.40  0.00  0.60
  Delta     0.80  0.70  0.60  0.00
  ```
  The script uses Biopython's Nexus parser, which is flexible for these formats. For best results with standalone PHYLIP files, ensure they strictly follow the format (number of taxa on the first line, then matrix data).

**Output:**
- A phylogenetic tree file in Newick format.

**Usage:**
```bash
python3 pylib/scripts/build_tree_from_distances.py <distance_matrix_file> [--method <method>] [--outfile <output_tree_file>] [--informat <format>]
```

**Arguments:**
- `<distance_matrix_file>`: (Required) Path to the input distance matrix file.
- `--method <method>`: Tree construction method. Choices:
    - `nj` (Neighbor-Joining, default)
    - `upgma` (UPGMA)
- `--outfile <output_tree_file>`: Output file name for the tree in Newick format. (Default: `phylogenetic_tree.nwk`)
- `--informat <format>`: Input distance matrix file format.
    - `nexus` (Default. Used for PHYLIP-style matrices, which can be standalone or within a Nexus `DISTANCES` block).

**Example:**
Assuming `my_distances.phy` contains a PHYLIP-formatted distance matrix:
```bash
# Build a Neighbor-Joining tree
python3 pylib/scripts/build_tree_from_distances.py my_distances.phy --method nj --outfile nj_tree.nwk

# Build a UPGMA tree
python3 pylib/scripts/build_tree_from_distances.py my_distances.phy --method upgma --outfile upgma_tree.nwk
```

**Notes:**
- This script relies on Biopython's `Bio.Phylo` and `Bio.Nexus` modules.
- Ensure your input distance matrix is correctly formatted. The Nexus parser can be particular. If you have issues with a standalone PHYLIP file, wrapping it in a minimal Nexus block might help:
  ```nexus
  #NEXUS
  BEGIN DISTANCES;
    FORMAT TRIANGLE=LOWER NODIAGONAL;
    MATRIX
      Alpha
      Beta      0.20
      Gamma     0.50  0.40
      Delta     0.80  0.70  0.60
    ;
  END;
  ```
---

## `pylib/scripts/analyze_msa.py`

Analyzes a Multiple Sequence Alignment (MSA) file to calculate a consensus sequence, derive basic statistics, or convert the alignment to a different file format.

**Input:**
- A Multiple Sequence Alignment (MSA) file in a format supported by Biopython's `Bio.Align.parse` (e.g., FASTA, Clustal, PHYLIP, Nexus, Stockholm).

**Output:**
- Text-based report (consensus sequence, statistics) printed to standard output or a specified file (`--output_report`).
- Optionally, a new file containing the MSA in a different format (if conversion is requested).

**Usage:**
```bash
python3 pylib/scripts/analyze_msa.py <msa_file> --informat <format> [task_options] [--output_report <report_file>]
```

**Arguments:**

*   `<msa_file>`: (Required) Path to the input Multiple Sequence Alignment file.
*   `--informat <format>`: (Required) Format of the input MSA file (e.g., `fasta`, `clustal`, `phylip`, `nexus`, `stockholm`).

**Task Options (At least one task must be specified):**

*   `--get_consensus`: Calculate and display the consensus sequence.
    *   `--consensus_threshold <float>`: (Optional, default: 0.5) Minimum frequency for a character to be considered in the consensus for a column.
    *   `--consensus_ambiguous_char <char>`: (Optional, default: 'X') Character to use for positions in the consensus where no single character meets the threshold, or where diversity meets `consensus_require_multiple`. The script attempts to change this to 'N' if the input is inferred as DNA/RNA and the default 'X' is active.
    *   `--consensus_require_multiple <int>`: (Optional, default: 1) If the number of different characters at a position (excluding gaps) is greater than or equal to this value, the `consensus_ambiguous_char` will be used even if one character meets the `consensus_threshold`. This helps represent highly variable but not entirely ambiguous positions. Set to a high value (e.g., 100) to effectively disable this and rely solely on the threshold.

*   `--get_stats`: Calculate and display basic MSA statistics. Output includes:
    *   Number of sequences.
    *   Alignment length (columns).
    *   Overall GC content (if input is inferred as DNA/RNA, otherwise "N/A").
    *   Percentage of columns consisting entirely of gaps.
    *   Percentage of columns containing at least one gap.

*   `--convert_to <output_format>`: Convert the input MSA to a different format.
    *   `--outfile <converted_msa_filepath>`: (Required if `--convert_to` is used) Path to save the converted MSA file. The `<output_format>` should be a format supported by Biopython's `Bio.Align.write`.

**General Optional Arguments:**

*   `--output_report <report_filepath>`: Path to a file where the text-based output (consensus, statistics) will be written. If not provided, output is printed to the standard output.

**Examples:**

1.  **Get consensus sequence from a Clustal alignment:**
    ```bash
    python3 pylib/scripts/analyze_msa.py my_alignment.aln --informat clustal --get_consensus
    ```

2.  **Get consensus with a 70% threshold and 'N' as ambiguous character:**
    ```bash
    python3 pylib/scripts/analyze_msa.py my_dna_alignment.fasta --informat fasta --get_consensus --consensus_threshold 0.7 --consensus_ambiguous_char N
    ```

3.  **Get MSA statistics:**
    ```bash
    python3 pylib/scripts/analyze_msa.py my_alignment.phy --informat phylip --get_stats
    ```

4.  **Convert a FASTA alignment to PHYLIP format:**
    ```bash
    python3 pylib/scripts/analyze_msa.py my_alignment.fasta --informat fasta --convert_to phylip --outfile my_alignment.phy
    ```

5.  **Get both consensus and stats, and save the report to a file:**
    ```bash
    python3 pylib/scripts/analyze_msa.py my_alignment.stk --informat stockholm --get_consensus --get_stats --output_report analysis_report.txt
    ```

**Notes:**
- The script processes the first alignment found in the input file.
- Alphabet inference (DNA/RNA vs. Protein) is basic and used for GC content calculation and defaulting the ambiguous consensus character.
- For consensus generation, gap characters (`-`, `.`) in a column are ignored when calculating frequencies but a column of only gaps will result in a gap in the consensus.

---

## `pylib/scripts/scan_sequences_for_motif.py`

Scans sequences in a FASTA file for occurrences of a known motif (simple string or IUPAC degenerate string) on both strands.

**Input:**
- A FASTA file containing nucleotide sequences.
- A motif string (e.g., `GATTACA` or `RRYYTGCAAW`).

**Output:**
- A tab-delimited text report listing each found motif occurrence, including sequence ID, motif, start and end positions (1-based), strand, and the matched sequence segment.
- Output is to standard output by default, or to a specified file.

**Usage:**
```bash
python3 pylib/scripts/scan_sequences_for_motif.py <sequence_file> --motif <motif_string> [--informat <format>] [--output_report <report_filepath>]
```

**Arguments:**

*   `<sequence_file>`: (Required) Path to the input FASTA file containing sequences to scan.
*   `--motif <motif_string>`: (Required) The motif to search for. This can be a simple DNA sequence (e.g., `AGCT`) or a degenerate IUPAC string (e.g., `RRYYNNNNNGATT`).
*   `--informat <format>`: (Optional) Format of the input sequence file. Default: `fasta`. (This argument is passed to Biopython's `SeqIO.parse`).
*   `--output_report <report_filepath>`: (Optional) Path to a file where the output report will be written. If not provided, the report is printed to standard output.

**Output Report Format (Tab-delimited):**

`Sequence_ID	Motif_Found	Start	End	Strand	Matched_Sequence`

-   `Sequence_ID`: Identifier of the sequence from the input FASTA file.
-   `Motif_Found`: The motif string provided by the user.
-   `Start`: 1-based start position of the motif match on the original forward strand.
-   `End`: 1-based end position of the motif match on the original forward strand.
-   `Strand`: `+` if found on the forward strand, `-` if found on the reverse complement strand.
-   `Matched_Sequence`: The actual DNA sequence segment that matched the motif.

**Example:**

Given a FASTA file `my_sequences.fasta`:
```fasta
>SeqA
GATTACAGATTACANNNRRYYTGCAAWGATTACA
>SeqB
CCCTTTGCAAYYKMNNNGTAATCG
```

Command:
```bash
python3 pylib/scripts/scan_sequences_for_motif.py my_sequences.fasta --motif RRYYTGCAAW --output_report motif_scan_results.tsv
```

Expected content of `motif_scan_results.tsv` (example):
```tsv
Sequence_ID	Motif_Found	Start	End	Strand	Matched_Sequence
SeqA	RRYYTGCAAW	19	28	+	AGCTTGCAAT
SeqB	RRYYTGCAAW	7	16	-	ATTGCAAGGG
```
*(Note: The exact matched sequence for the reverse strand will be the sequence as it appeared on the reverse complement, matching the forward version of the motif string).*

A summary message indicating the total number of motifs found and sequences processed will be printed to standard error.

**Notes:**
- The script converts the IUPAC motif string to a regular expression for searching.
- Searches are case-insensitive.
- The script will report all non-overlapping matches.