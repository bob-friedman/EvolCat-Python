# Standalone Tools

This directory contains standalone, supported Python scripts for various sequence analysis tasks. These tools are organized into functional categories.

## Directory Structure

- **`metrics/`**: Scripts for calculating evolutionary metrics.
- **`fasta_utils/`**: General utilities for handling FASTA files.
- **`msa_utils/`**: Utilities for analyzing and manipulating Multiple Sequence Alignments (MSA).
- **`test_data/`**: Sample datasets and configuration files for testing the scripts.

---

## 1. Metrics (`metrics/`)

These scripts estimate substitution rates and other metrics. They have been tested with theoretical data but require further validation for empirical use.

### `dsdn_dist.py`
Estimates dS/dN (synonymous to non-synonymous substitution rate) metrics using the Nei & Gojobori (1986) method. See the [Evolutionary Metrics Guide](../../../guides/evolutionary_metrics_guide.md) for more details.

**Usage:**
```bash
python metrics/dsdn_dist.py <input_fasta> [ratio|R]
```

### `pnc_pnr_dist.py`
Estimates pNC/pNR (conservative to radical non-synonymous substitution rate) metrics using the Hughes et al. (1990) method. See the [Evolutionary Metrics Guide](../../../guides/evolutionary_metrics_guide.md) for more details on the method and its application.

**Usage:**
```bash
python metrics/pnc_pnr_dist.py <input_fasta> <property_file> [ratio|R]
```
*Sample property file `test_data/property_charge` is provided.*

---

## 2. FASTA Utilities (`fasta_utils/`)

General-purpose scripts for processing FASTA files.

### `clean_fasta_name.py`
Cleans FASTA headers by replacing whitespace and periods with underscores and converting to uppercase.
**Usage:**
```bash
python fasta_utils/clean_fasta_name.py <input_fasta_file> > <output_fasta_file>
```

### `merge_fastas.py`
Merges multiple FASTA files. Sequences with identical headers are concatenated.
**Usage:**
```bash
python fasta_utils/merge_fastas.py <file1.fasta> <file2.fasta> ... > <merged.fasta>
```

### `translate_seq.py`
Translates nucleotide sequences from a FASTA file to protein sequences.
**Usage:**
```bash
python fasta_utils/translate_seq.py --input_file <input_fasta_file> [options]
```

---

## 3. MSA Utilities (`msa_utils/`)

Tools for analyzing or converting sequence alignments.

### `analyze_msa.py`
Analyzes an MSA file to calculate a consensus sequence, compute basic statistics, or convert to a different format.
**Usage:**
```bash
python msa_utils/analyze_msa.py <msa_file> --informat <format> [options]
```

### `nogaps.py`
Removes columns containing non-alphabetic characters (gaps, etc.) from an aligned FASTA file.
**Usage:**
```bash
python msa_utils/nogaps.py <input_fasta_file> > <output_fasta_file>
```

### `pal2nal_enforce.py`
Enforces an amino acid sequence alignment onto corresponding unaligned nucleotide sequences, producing a codon-aware nucleotide alignment.
**Usage:**
```bash
python msa_utils/pal2nal_enforce.py -a <aligned_aa.fasta> -n <unaligned_nt.fasta> -o <output_nt_aligned.fasta>
```
*Test data: `test_data/test_aa_aligned.fasta`, `test_data/test_nt_unaligned.fasta`.*

### `phy2fas.py`
Converts a PHYLIP file (interleaved or sequential) to FASTA format.
**Usage:**
```bash
python msa_utils/phy2fas.py <input_phylip_file> > <output_fasta_file>
```

---

## Dependencies

Most scripts depend on **BioPython**. You can install it via:
```bash
pip install biopython
```
Some scripts are standalone and do not require external libraries.
