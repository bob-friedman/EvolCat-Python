# Test Suite for `extract_cds_region.py`

## Overview

This directory contains tests for the `extract_cds_region.py` script located in `pylib/scripts/`. The primary test script is `test_extract_cds_region.py`.

## Purpose

The `test_extract_cds_region.py` script aims to ensure the correctness and robustness of the `extract_cds_region.py` script across various scenarios. It validates both the "Range Extraction" mode and the "Single Position Analysis" mode.

## Running Tests

The tests are built using Python's `unittest` framework. To run the tests:

1.  Navigate to the root directory of this repository.
2.  Execute the test script directly or use `unittest` discovery:

    **Directly:**
    ```bash
    python -m pylib.tests.test_extract_cds_region
    ```

    **Using unittest discovery (from the project root directory):**
    ```bash
    python -m unittest discover -s pylib/tests
    ```
    Or, if `pylib` is in your `PYTHONPATH`:
    ```bash
    python -m unittest discover -s pylib.tests
    ```

    If you are in the `pylib/tests/` directory, you might be able to run:
    ```bash
    python -m unittest test_extract_cds_region.py
    ```

The tests will run, and `unittest` will report the status of each test (e.g., pass, fail, error), along with a summary.

## Test Case Coverage

The test suite `TestExtractCDSRegion` in `test_extract_cds_region.py` covers a range of scenarios, including but not limited to:

### Range Extraction Mode (`--start_pos`, `--end_pos`):
-   **CDS on Forward Strand**:
    -   Query fully overlaps the CDS.
-   **CDS on Reverse Strand**:
    -   Query fully overlaps the CDS (handles `complement()` locations).
-   **Partial Overlaps**:
    -   Query overlaps only the beginning of a CDS.
    -   Query overlaps only the end of a CDS (implicit, covered by other partials).
    -   Query is fully contained within a CDS but smaller than the CDS.
    -   CDS is fully contained within the query range.
-   **No CDS Overlap**:
    -   Query region does not intersect with any CDS features.
-   **Multiple CDS Features**:
    -   Handling of files with more than one CDS (implicit in processing logic).
-   **Compound Locations (Exons)**:
    -   Query overlapping parts of a CDS defined by `join()` (e.g., spanning exon boundaries).
-   **Frame Adjustments**:
    -   CDS with `codon_start` not equal to 1.
    -   Extraction resulting in partial codons at the beginning/end of the translated segment.
-   **Invalid Inputs**:
    -   Start position greater than end position.
    -   Query range out of the bounds of the sequence record.
    -   Non-positive start/end positions.
-   **Edge Cases**:
    -   Querying very short regions.
    -   CDS at the very beginning or end of a sequence.

### Single Position Mode (`--single_position`):
-   **Position within CDS (Forward Strand)**:
    -   Target nucleotide in the middle of a codon.
    -   Target nucleotide at the start/end of a codon.
    -   CDS with `codon_start=2` or `codon_start=3`.
-   **Position within CDS (Reverse Strand)**:
    -   Target nucleotide in the middle of a codon on a reverse strand CDS.
-   **Compound Locations (Exons)**:
    -   Target nucleotide within an exon of a joined CDS (e.g., in the second exon).
-   **Boundary Conditions**:
    -   Target nucleotide is part of an incomplete codon at the very end of a CDS.
    -   Target nucleotide is within a CDS feature but *before* the actual translation start indicated by `codon_start`.
-   **Position Not in CDS**:
    -   Target nucleotide falls outside any annotated CDS feature.
-   **Invalid Inputs**:
    -   Non-positive position.
    -   Position out of sequence bounds.

Each test method typically creates a temporary GenBank file with specific content tailored to the scenario being tested, runs `extract_cds_region.py` as a subprocess, and then asserts expectations on the script's standard output and standard error.
```
