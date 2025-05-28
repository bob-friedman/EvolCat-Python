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

THEN
scan_sequences_for_motif.py

#!/usr/bin/env python3

"""
Scans sequences in a FASTA file for a given motif (IUPAC or simple string)
and reports all occurrences.
"""

import argparse
import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq # For reverse_complement

# --- Helper Function: IUPAC to Regex ---
def iupac_to_regex(iupac_string):
    """Converts an IUPAC string to a Python regex pattern."""
    iupac_map = {
        'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'U', # U is often treated as T in DNA contexts
        'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
        'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
        'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]',
    }
    regex_pattern = []
    for char in iupac_string.upper(): # Convert to uppercase to match keys
        if char in iupac_map:
            regex_pattern.append(iupac_map[char])
        else:
            # Treat unknown characters as literals and warn
            sys.stderr.write(f"Warning: Unknown IUPAC code '{char}' in motif treated as literal.\n")
            regex_pattern.append(re.escape(char)) # Escape to be safe if it's a regex special char
    return "".join(regex_pattern)

# --- Main Logic ---
def main():
    parser = argparse.ArgumentParser(
        description="Scan sequences for a given motif (IUPAC or simple string).",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "sequence_file",
        help="Path to the input sequence file (e.g., FASTA)."
    )
    parser.add_argument(
        "--motif",
        required=True,
        help="Motif string to search for (IUPAC or simple string)."
    )
    parser.add_argument(
        "--informat",
        default="fasta",
        help="Format of the input sequence file (default: fasta)."
    )
    parser.add_argument(
        "--output_report",
        default=None, # Default to None, will be handled as sys.stdout
        help="Path to the output report file (default: print to standard output)."
    )

    args = parser.parse_args()

    # Convert motif to regex
    try:
        regex_motif = iupac_to_regex(args.motif)
        compiled_regex = re.compile(regex_motif, re.IGNORECASE)
    except Exception as e:
        sys.stderr.write(f"Error: Could not compile regex from motif '{args.motif}'. Details: {e}\n")
        sys.exit(1)

    # Setup output stream
    output_stream = sys.stdout
    if args.output_report:
        try:
            output_stream = open(args.output_report, 'w')
        except IOError as e:
            sys.stderr.write(f"Error: Could not open output report file '{args.output_report}': {e}\n")
            sys.exit(1)

    # Write header
    header = "Sequence_ID\tMotif_Found\tStart\tEnd\tStrand\tMatched_Sequence\n"
    output_stream.write(header)

    total_motifs_found = 0
    sequences_processed = 0

    try:
        for record in SeqIO.parse(args.sequence_file, args.informat):
            sequences_processed += 1
            seq_str = str(record.seq)
            seq_len = len(seq_str)

            # Forward Strand Search
            for match in compiled_regex.finditer(seq_str):
                total_motifs_found += 1
                start_1based = match.start() + 1
                end_1based = match.end() # match.end() is already exclusive, so it's correct for 1-based end
                output_stream.write(
                    f"{record.id}\t{args.motif}\t{start_1based}\t{end_1based}\t+\t{match.group(0)}\n"
                )

            # Reverse Complement Strand Search
            rc_seq_str = str(record.seq.reverse_complement())
            for match in compiled_regex.finditer(rc_seq_str):
                total_motifs_found += 1
                # Convert coordinates back to forward strand
                # match.start() and match.end() are 0-based on the rc_seq_str
                rc_match_start_0based = match.start()
                rc_match_end_0based = match.end()

                original_start_1based = seq_len - rc_match_end_0based + 1
                original_end_1based = seq_len - rc_match_start_0based
                
                output_stream.write(
                    f"{record.id}\t{args.motif}\t{original_start_1based}\t{original_end_1based}\t-\t{match.group(0)}\n"
                )
                
    except FileNotFoundError:
        sys.stderr.write(f"Error: Sequence file '{args.sequence_file}' not found.\n")
        if args.output_report: output_stream.close()
        sys.exit(1)
    except ValueError as e: # Biopython often raises ValueError for format issues
        sys.stderr.write(f"Error parsing sequence file '{args.sequence_file}' with format '{args.informat}'. Details: {e}\n")
        if args.output_report: output_stream.close()
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"An unexpected error occurred: {e}\n")
        if args.output_report: output_stream.close()
        sys.exit(1)
    finally:
        if args.output_report and output_stream is not sys.stdout:
            output_stream.close()

    sys.stderr.write(f"Found {total_motifs_found} motif(s) in {sequences_processed} sequence(s).\n")

if __name__ == '__main__':
    main()
```
