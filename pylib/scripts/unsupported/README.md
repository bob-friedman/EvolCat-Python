<a name="top"></a>
# Unsupported Tools

Included here are standalone unsupported Python scripts for various sequence analysis tasks.

## Metrics
- `dsdn_dist.py`: Estimates dS/dN metrics.
- `pnc_pnr_dist.py`: Estimates pNC/pNR metrics.

There is a sample test dataset (`test.fasta`) and a file for the amino acid property of interest (`property_charge`, used for pNC/pNR). These metrics and their calculations were tested by use of theoretical sequence data and reference code as a verification guide, but they require further testing for use in production since the testing was not yet expanded to empirical sequence data. There is also some error checking that the sequence data is in a format that is expected from the code.

## General Utilities
These scripts are standalone versions of tools found in the main `pylib/scripts` directory.

### `analyze_msa.py`
Analyzes a Multiple Sequence Alignment (MSA) file to calculate a consensus sequence, compute basic statistics, or convert to a different MSA format.
**Usage:**
```bash
python analyze_msa.py <msa_file> --informat <format> [options]
```

### `clean_fasta_name.py`
Cleans FASTA headers by replacing whitespace and periods with underscores and converting to uppercase.
**Usage:**
```bash
python clean_fasta_name.py <input_fasta_file> > <output_fasta_file>
```

### `merge_fastas.py`
Merges multiple FASTA files. Sequences with identical headers are concatenated.
**Usage:**
```bash
python merge_fastas.py <file1.fasta> <file2.fasta> ... > <merged.fasta>
```

### `pal2nal_enforce.py`
Enforces an amino acid sequence alignment onto corresponding unaligned nucleotide sequences, producing a codon-aware nucleotide alignment. Gaps in the AA alignment are converted to '---' in the nucleotide output. Sample data files `test_aa_aligned.fasta`, `test_nt_unaligned.fasta`, and `test_nt_aligned_output.fasta` (reference output) are provided for testing.
**Usage:**
```bash
python pal2nal_enforce.py -a <aligned_aa.fasta> -n <unaligned_nt.fasta> -o <output_nt_aligned.fasta>
```

### `nogaps.py`
Removes columns containing non-alphabetic characters (gaps, etc.) from an aligned FASTA file.
**Usage:**
```bash
python nogaps.py <input_fasta_file> > <output_fasta_file>
```

### `phy2fas.py`
Converts a PHYLIP file (interleaved or sequential) to FASTA format.
**Usage:**
```bash
python phy2fas.py <input_phylip_file> > <output_fasta_file>
```

### `translate_seq.py`
Translates nucleotide sequences from a FASTA file to protein sequences.
**Usage:**
```bash
python translate_seq.py --input_file <input_fasta_file> [options]
```
