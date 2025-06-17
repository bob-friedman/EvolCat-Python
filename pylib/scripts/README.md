# extract_cds_region.py

## Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Dependencies](#dependencies)
- [Input Format](#input-format)
- [Usage](#usage)
  - [Mode 1: Range Extraction](#mode-1-range-extraction)
    - [Arguments](#arguments-range-extraction)
    - [Example](#example-range-extraction)
    - [Output Format (Range Extraction)](#output-format-range-extraction)
  - [Mode 2: Single Position Analysis](#mode-2-single-position-analysis)
    - [Arguments](#arguments-single-position-analysis)
    - [Example](#example-single-position-analysis)
    - [Output Format (Single Position Analysis)](#output-format-single-position-analysis)
- [Error Handling](#error-handling)
- [Notes](#notes)

## Overview

`extract_cds_region.py` is a Python script designed to process GenBank files to extract and analyze coding sequence (CDS) information. It offers two primary modes of operation:

1.  **Range Extraction Mode**: Extracts a nucleotide sequence from a specified region (start and end positions) within a GenBank record. If this region overlaps with any CDS features, the script will also extract the corresponding nucleotide sequence from the CDS and translate it into its protein sequence.
2.  **Single Position Mode**: Identifies the codon and corresponding amino acid at a specific nucleotide position within a GenBank record, provided that position falls within a CDS feature.

The script leverages the BioPython library for parsing GenBank files and performing sequence manipulations.

## Features

-   Parses GenBank files (can handle multi-record files).
-   Extracts nucleotide sequences from user-defined genomic ranges.
-   Identifies CDS features overlapping with the queried range.
-   Extracts corresponding CDS nucleotide sequences, handling strand information (forward/reverse) and compound locations (e.g., exons).
-   Translates overlapping CDS sequences to protein sequences using the appropriate NCBI translation table specified in the GenBank record (defaults to table 1).
-   Pinpoints the codon and translated amino acid for a single genomic nucleotide position.
-   Handles `codon_start` qualifiers in CDS features for correct frame determination.
-   Provides informative output, including record ID, CDS ID, nucleotide sequences, and protein sequences.
-   Includes error handling for invalid inputs, file issues, and out-of-bounds queries.

## Dependencies

-   **Python 3**: The script is written for Python 3.
-   **BioPython**: Used for GenBank parsing, sequence objects, and translation. Ensure BioPython is installed in your Python environment. You can typically install it using pip:
    ```bash
    pip install biopython
    ```

## Input Format

The script requires a GenBank file (`.gb`, `.gbk`, `.genbank`) as input. This file should conform to the standard GenBank format. Key information used by the script from the GenBank file includes:
-   The main nucleotide sequence (`ORIGIN`).
-   CDS features in the `FEATURES` table, including:
    -   `location` (e.g., `51..62`, `complement(50..149)`, `join(51..80,101..130)`)
    -   `protein_id` or `locus_tag` (used as CDS identifier)
    -   `codon_start` (e.g., `1`, `2`, or `3`; defaults to `1` if not present)
    -   `transl_table` (NCBI translation table ID; defaults to `1` if not present)

## Usage

The script is run from the command line.

### Mode 1: Range Extraction

This mode is activated when `--start_pos` and `--end_pos` arguments are provided.

#### Arguments (Range Extraction)

```
usage: extract_cds_region.py genbank_file --start_pos START_POS --end_pos END_POS

Required arguments:
  genbank_file          Path to the input GenBank file.
  --start_pos START_POS Start base pair of the region to extract (1-based).
  --end_pos END_POS     End base pair of the region to extract (1-based).
```

#### Example (Range Extraction)

Suppose you have a GenBank file named `my_sequence.gb` and you want to extract the region from base 50 to 150:

```bash
python extract_cds_region.py my_sequence.gb --start_pos 50 --end_pos 150
```

#### Output Format (Range Extraction)

The script prints output to standard output. For each record in the GenBank file:
1.  It indicates which record is being processed (e.g., `Processing record: TestEntry (length 202)`).
2.  It prints the nucleotide sequence from the user-specified query range (e.g., `Nucleotide sequence in query range (51-62): ATGAAACTTAAT`).
3.  For each CDS feature found:
    -   It prints the CDS ID and its location (e.g., `Found CDS: XYZ_001 located at [50:62]`).
    -   If the CDS overlaps with the query range, it prints the nucleotide sequence derived from the CDS that corresponds to the overlapping part (e.g., `Nucleotide sequence from CDS (XYZ_001) overlapping query: ATGAAACTTAAT`).
    -   It then prints the translated protein sequence from this overlapping CDS portion, indicating the frame adjustment and translation table used (e.g., `Protein sequence from CDS XYZ_001 (frame adj: 0, table: 1): MKLN`).
4.  If no CDS features overlap with the query range in a record, a message like `No CDS features found whose exonic parts overlap with the query range...` is printed.

Warnings or errors (e.g., sequence containing 'N' characters, translation errors) are printed to standard output or standard error.

### Mode 2: Single Position Analysis

This mode is activated when the `--single_position` argument is provided.

#### Arguments (Single Position Analysis)

```
usage: extract_cds_region.py genbank_file --single_position SINGLE_POSITION

Required arguments:
  genbank_file            Path to the input GenBank file.
  --single_position SINGLE_POSITION
                          Specify a single 1-based nucleotide position to find its codon.
```
*Note: If `--start_pos` or `--end_pos` are provided along with `--single_position`, they will be ignored.*

#### Example (Single Position Analysis)

To find information about nucleotide position 57 in `my_sequence.gb`:

```bash
python extract_cds_region.py my_sequence.gb --single_position 57
```

#### Output Format (Single Position Analysis)

For each record in the GenBank file where the target position falls within a CDS:
1.  Processing information for the record and target position.
2.  If the position is found within a CDS feature, detailed information is printed:
    ```
    ---------------------------------------------------
    Record ID:                    TestSingleMidFwd.1
    CDS ID:                       XYZ_S1
    CDS Location:                 [50:150]
    Target Nucleotide Position:   57 (genomic, 1-based)
    CDS Nucleotide Index:         6 (0-based, within full CDS sequence)
    Target Base in Codon:       1
    Codon Sequence:               CTT
    Translation Table ID:         1
    Translated Codon:             L
    ---------------------------------------------------
    ```
    -   **Record ID**: Identifier of the GenBank record.
    -   **CDS ID**: Identifier of the CDS feature (protein_id or locus_tag).
    -   **CDS Location**: Genomic location of the CDS feature.
    -   **Target Nucleotide Position**: The user-specified 1-based genomic position.
    -   **CDS Nucleotide Index**: The 0-based index of the target nucleotide within the conceptual full coding sequence of the CDS (after accounting for strand and joins).
    -   **Target Base in Codon**: The position (1, 2, or 3) of the target nucleotide within its codon.
    -   **Codon Sequence**: The three-nucleotide sequence of the codon containing the target base.
    -   **Translation Table ID**: The NCBI translation table used.
    -   **Translated Codon**: The single-letter amino acid code for the translated codon (or 'X' for unknown, '*' for stop, or an error message).
3.  If the target position is within a CDS but occurs before the effective translation start (due to `codon_start` qualifier), a message indicating this is printed, and codon/translation details are omitted.
4.  If the target position is part of an incomplete codon at the end of a CDS, this is noted.
5.  If the target position is not found within any CDS feature in any processed record, a message like `Target nucleotide position X was not found within any CDS feature...` is printed.

## Error Handling

The script includes checks for:
-   File not found.
-   Invalid start/end positions (e.g., start > end, non-positive values).
-   Query range out of bounds for a record.
-   Invalid single position (e.g., non-positive).
-   Use of invalid NCBI translation table IDs in GenBank features.
-   Other unexpected errors during processing.

Error messages are typically printed to standard error (`stderr`). For critical errors like file not found or invalid arguments, the script may exit with a non-zero status code.

## Notes

-   Positions are 1-based for user input (`--start_pos`, `--end_pos`, `--single_position`) and in most output messages referring to genomic coordinates, consistent with common biological conventions. Internally, 0-based indexing is often used for sequence manipulation.
-   The script assumes the GenBank file is correctly formatted. Malformed files might lead to parsing errors or unexpected behavior.
-   For CDS features with compound locations (e.g., `join(exon1, exon2)`), the script correctly reconstructs the full coding sequence before processing.
-   For reverse strand CDS features (e.g., `complement(...)`), the script correctly uses the reverse complement for extraction and translation.
-   The "frame adjustment" mentioned in the range extraction output refers to how the extracted CDS portion aligns with codon boundaries, considering the `codon_start` qualifier and the offset of the query within the CDS.
```
