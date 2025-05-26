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

## `pylib/scripts/paml_tools/calculate_ds_dn.py`

Calculates synonymous (dS) and non-synonymous (dN) substitution rates (more accurately, counts of synonymous/non-synonymous substitutions and potential sites) between pairs of DNA sequences from a FASTA file.

Calculates pairwise dN/dS ratios for sequences from an aligned FASTA file.

This script utilizes BioPython to interface with PAML's `yn00` program to perform
the dN/dS calculations. The script then parses the output generated by `yn00`
to extract dN (non-synonymous substitution rate), dS (synonymous substitution
rate), and omega (dN/dS ratio) for each pair of sequences.

**Usage:**
```bash
python pylib/scripts/calculate_ds_dn.py <input_fasta_file> [ratio] [--genetic_code <standard|mito>]
```

Arguments:
  fasta_file: Path to the input FASTA file. This file must contain two or more
              coding sequences that have already been aligned.

PAML Dependency:
  For actual dN/dS calculations (i.e., not simulating), this script requires
  PAML (Phylogenetic Analysis by Maximum Likelihood) to be installed, and the
  `yn00` executable from PAML must be in the system's PATH. If `yn00` is not
  found, the script will output an error. PAML can be obtained from:
  http://abacus.gene.ucl.ac.uk/software/paml.html

**Arguments:**
- `<input_fasta_file>`: Path to the input FASTA file containing aligned DNA sequences. Sequences must be of equal length after removing gaps/Ns and be a multiple of 3.
- `[ratio]`: Optional. The transition/transversion ratio (R). Can be a floating-point number (e.g., 0.5, 2.0) or the string 'R' to estimate the ratio from the data. Defaults to 0.5.
- `--genetic_code <standard|mito>`: Optional. The genetic code to use. Choices are 'standard' or 'mito'. Defaults to 'standard'.

**Output Format (Tab-delimited):**
`gene1_name	gene2_name	formatted_R_ratio	syn_subs_count	nonsyn_subs_count	potential_syn_sites	potential_nonsyn_sites`
- `formatted_R_ratio`: R value used, formatted to two decimal places, or "undef" if estimation failed.
- `syn_subs_count`: Count of synonymous substitutions.
- `nonsyn_subs_count`: Count of non-synonymous substitutions.
- `potential_syn_sites`: Number of potential synonymous sites.
- `potential_nonsyn_sites`: Number of potential non-synonymous sites.
- If denominators are zero during calculation for a pair, specific messages like "denom_0_syn=..." may be printed for that pair instead of the counts.

**Example:**
See TESTING_GUIDE in the directory containing this script.
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
