import unittest
from unittest.mock import patch, MagicMock, mock_open, call
import os
import sys
import io
from contextlib import redirect_stdout, redirect_stderr

# Add paths for script import
current_dir = os.path.dirname(os.path.abspath(__file__))
pylib_root = os.path.dirname(os.path.dirname(current_dir))
if pylib_root not in sys.path:
    sys.path.insert(0, pylib_root)

# Import the script to be tested
try:
    from pylib.scripts import nogaps
except ModuleNotFoundError:
    scripts_dir = os.path.join(pylib_root, "pylib", "scripts")
    if scripts_dir not in sys.path:
        sys.path.insert(0, scripts_dir)
    import nogaps


# Mock Biopython's SeqRecord for creating test sequences
class MockSeqRecord:
    def __init__(self, seq_str, id="test_id", name="", description="test_desc"):
        self.seq = seq_str  # The script accesses record.seq directly as a string
        self.id = id
        self.name = name if name else id
        self.description = description

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        # This simplified __str__ is for easier assertion if we check string representation of the mock.
        # Bio.SeqIO.write will primarily use .id, .description, .seq attributes.
        return f"ID: {self.id}, Seq: {self.seq}"


class TestNoGaps(unittest.TestCase):

    def setUp(self):
        pass

    def _run_main_with_mocked_parse(
        self,
        fasta_file_arg,
        parse_fasta_return_value=None,
        parse_fasta_side_effect=None,
    ):
        """Helper to run main with mocked parse_fasta_file and captured output."""
        full_argv = ["nogaps.py", fasta_file_arg]

        with patch.object(sys, "argv", full_argv), patch(
            "pylib.utils.seq_parser.parse_fasta_file",
            return_value=parse_fasta_return_value,
            side_effect=parse_fasta_side_effect,
        ) as mock_parse, patch("Bio.SeqIO.write") as mock_seqio_write, redirect_stdout(
            io.StringIO()
        ) as captured_stdout, redirect_stderr(
            io.StringIO()
        ) as captured_stderr, patch(
            "sys.exit"
        ) as mock_sys_exit:

            try:
                nogaps.main()
            except SystemExit:
                pass

            return {
                "stdout": captured_stdout.getvalue(),
                "stderr": captured_stderr.getvalue(),
                "mock_exit": mock_sys_exit,
                "mock_parse": mock_parse,
                "mock_seqio_write": mock_seqio_write,
            }

    # --- 1. Argument Parsing Tests (Covered in previous turn, re-verified) ---
    def test_arg_parser_input_file_required(self):
        with patch.object(sys, "argv", ["nogaps.py"]):
            with self.assertRaises(SystemExit) as cm, redirect_stderr(io.StringIO()):
                nogaps.main()
            self.assertEqual(cm.exception.code, 2)

    # --- 2. FASTA Input and Validation ---
    def test_input_file_not_found(self):
        """Test FileNotFoundError for the input FASTA file."""
        results = self._run_main_with_mocked_parse(
            "non_existent.fasta", parse_fasta_side_effect=FileNotFoundError
        )
        results["mock_parse"].assert_called_once_with("non_existent.fasta")
        self.assertIn("Error: File 'non_existent.fasta' not found.", results["stderr"])
        results["mock_exit"].assert_called_once_with(1)

    def test_empty_input_file_no_sequences(self):
        """Test with an empty input file or a file with no sequences."""
        results = self._run_main_with_mocked_parse(
            "empty.fasta", parse_fasta_return_value=[]
        )
        self.assertIn("Error: No sequences found in 'empty.fasta'.", results["stderr"])
        results["mock_exit"].assert_called_once_with(1)

    def test_sequences_of_varying_lengths_error(self):
        """Test with sequences of varying lengths."""
        mock_records = [
            MockSeqRecord("ACGT", "seq1"),
            MockSeqRecord("ACG", "seq2"),  # Shorter
        ]
        results = self._run_main_with_mocked_parse(
            "varied.fasta", parse_fasta_return_value=mock_records
        )
        self.assertIn(
            "Error: Sequences in 'varied.fasta' are not of equal length.",
            results["stderr"],
        )
        results["mock_exit"].assert_called_once_with(1)

    def test_valid_aligned_sequences_equal_length(self):
        """Test with valid aligned sequences (equal length), should proceed without length error."""
        mock_records = [MockSeqRecord("ACGT", "seq1"), MockSeqRecord("GATT", "seq2")]
        results = self._run_main_with_mocked_parse(
            "aligned.fasta", parse_fasta_return_value=mock_records
        )
        # Expect no error related to length, and SeqIO.write to be called.
        self.assertEqual(results["stderr"], "")
        results["mock_seqio_write"].assert_called()  # If processing completes

    def test_sequences_of_zero_length(self):
        """Test with sequences of zero length (all are zero length)."""
        mock_records = [MockSeqRecord("", "seq1"), MockSeqRecord("", "seq2")]
        results = self._run_main_with_mocked_parse(
            "zero_len.fasta", parse_fasta_return_value=mock_records
        )
        # Expect no error related to length, as they are equal.
        # Output should be empty sequences.
        self.assertEqual(results["stderr"], "")
        results["mock_seqio_write"].assert_called()
        # Check what was written: should be two records with empty sequences.
        written_records = results["mock_seqio_write"].call_args[0][0]
        self.assertEqual(len(list(written_records)), 2)
        # The generator needs to be consumed to check individual records
        output_recs_list = list(results["mock_seqio_write"].call_args[0][0])
        self.assertEqual(output_recs_list[0].seq, "")
        self.assertEqual(output_recs_list[1].seq, "")

    # --- 3. Gap/Invalid Column Identification Logic & 4. Sequence Reconstruction ---
    def test_no_invalid_columns(self):
        mock_records = [MockSeqRecord("ACGT", "s1"), MockSeqRecord("GATT", "s2")]
        results = self._run_main_with_mocked_parse(
            "no_gaps.fasta", parse_fasta_return_value=mock_records
        )
        written_records = list(results["mock_seqio_write"].call_args[0][0])
        self.assertEqual(written_records[0].seq, "ACGT")
        self.assertEqual(written_records[1].seq, "GATT")

    def test_all_columns_are_invalid(self):
        mock_records = [MockSeqRecord("----", "s1"), MockSeqRecord("????", "s2")]
        results = self._run_main_with_mocked_parse(
            "all_gaps.fasta", parse_fasta_return_value=mock_records
        )
        written_records = list(results["mock_seqio_write"].call_args[0][0])
        self.assertEqual(written_records[0].seq, "")
        self.assertEqual(written_records[1].seq, "")

    def test_invalid_columns_at_start_middle_end(self):
        # Col 0: - (invalid)
        # Col 1: A,G (valid)
        # Col 2: ?,N (invalid, N is treated as valid by isalpha) -> Script considers N valid.
        # Col 3: -,G (invalid)
        # Col 4: T,T (valid)
        # Col 5: *, (invalid)
        # Expected: A T for s1, G T for s2 (if N is invalid)
        # If N is valid by isalpha(): A N T for s1, G N T for s2.
        # Script uses str.isalpha() - this means N is VALID.
        # Invalid cols: 0, 3, 5
        # Valid cols: 1, 2, 4
        # s1: - A ? - T *  -> A ? T (if ? isalpha()=false, N isalpha()=true)
        # s2: - G N G T ,  -> G N T
        # The script's rule is `not char.isalpha()`. `?`, `-`, `*`, `,` are not alpha.
        # So column 2 (with N) is valid.
        # s1: -A?-T* -> AT (col 1='A', col 4='T')
        # s2: -GN-T, -> GT (col 1='G', col 4='T')
        # Let's re-evaluate based on `char.isalpha()`:
        # Col 0: '-' (s1), '-' (s2) -> invalid
        # Col 1: 'A' (s1), 'G' (s2) -> valid
        # Col 2: 'N' (s1), 'N' (s2) -> valid (N is alpha)
        # Col 3: '?' (s1), '-' (s2) -> invalid
        # Col 4: 'T' (s1), 'T' (s2) -> valid
        # Col 5: '*' (s1), ',' (s2) -> invalid
        # Valid columns are 1, 2, 4.
        # s1 sequence becomes: s1[1]s1[2]s1[4] = ANT
        # s2 sequence becomes: s2[1]s2[2]s2[4] = GNT
        mock_records = [MockSeqRecord("-AN?T*", "s1"), MockSeqRecord("-GN-T,", "s2")]
        results = self._run_main_with_mocked_parse(
            "mixed_gaps.fasta", parse_fasta_return_value=mock_records
        )
        written_records = list(results["mock_seqio_write"].call_args[0][0])
        self.assertEqual(written_records[0].seq, "ANT")
        self.assertEqual(written_records[1].seq, "GNT")

    def test_mix_of_different_invalid_chars(self):
        # Col 0: valid (A,G)
        # Col 1: invalid (-,?)
        # Col 2: valid (C,N) (N is alpha)
        # Col 3: invalid (*, whitespace)
        # Col 4: valid (T,T)
        # Expected: s1 -> ACT, s2 -> GNT
        mock_records = [MockSeqRecord("A-C*T", "s1"), MockSeqRecord("G?N T", "s2")]
        results = self._run_main_with_mocked_parse(
            "various_invalids.fasta", parse_fasta_return_value=mock_records
        )
        written_records = list(results["mock_seqio_write"].call_args[0][0])
        self.assertEqual(written_records[0].seq, "ACT")
        self.assertEqual(written_records[1].seq, "GNT")

    def test_numeric_chars_are_invalid(self):
        # Col 1: 1,2 (invalid)
        # Expected: s1 -> AC, s2 -> GT
        mock_records = [MockSeqRecord("A1C", "s1"), MockSeqRecord("G2T", "s2")]
        results = self._run_main_with_mocked_parse(
            "numeric_gaps.fasta", parse_fasta_return_value=mock_records
        )
        written_records = list(results["mock_seqio_write"].call_args[0][0])
        self.assertEqual(written_records[0].seq, "AC")
        self.assertEqual(written_records[1].seq, "GT")

    # --- 5. Output Verification ---
    def test_output_is_fasta_format_and_preserves_ids_desc(self):
        mock_records = [
            MockSeqRecord("A-C", "s1", description="desc1"),
            MockSeqRecord("G-T", "s2", description="desc2"),
        ]
        results = self._run_main_with_mocked_parse(
            "output_test.fasta", parse_fasta_return_value=mock_records
        )

        # Check that Bio.SeqIO.write was called
        results["mock_seqio_write"].assert_called_once()

        # The first argument to SeqIO.write is a generator of SeqRecord-like objects.
        # The second argument is the output handle (sys.stdout in this case).
        # The third is 'fasta'.
        gen_records, out_handle, fmt = results["mock_seqio_write"].call_args[0]
        self.assertEqual(fmt, "fasta")
        self.assertEqual(out_handle, sys.stdout)  # Script writes to sys.stdout

        # Consume the generator to check its contents
        output_records_list = list(gen_records)
        self.assertEqual(len(output_records_list), 2)

        # Check first record
        self.assertEqual(output_records_list[0].id, "s1")
        self.assertEqual(output_records_list[0].description, "desc1")
        self.assertEqual(output_records_list[0].seq, "AC")  # Reconstructed seq

        # Check second record
        self.assertEqual(output_records_list[1].id, "s2")
        self.assertEqual(output_records_list[1].description, "desc2")
        self.assertEqual(output_records_list[1].seq, "GT")

    def test_output_empty_if_all_columns_gapped_or_no_valid_input_seqs(self):
        # Case 1: All columns gapped
        mock_records_all_gaps = [MockSeqRecord("--", "s1"), MockSeqRecord("??", "s2")]
        results_all_gaps = self._run_main_with_mocked_parse(
            "all_gaps_output.fasta", parse_fasta_return_value=mock_records_all_gaps
        )
        written_records_all_gaps = list(
            results_all_gaps["mock_seqio_write"].call_args[0][0]
        )
        self.assertEqual(written_records_all_gaps[0].seq, "")
        self.assertEqual(written_records_all_gaps[1].seq, "")

        # Case 2: No input sequences (already tested, but confirm SeqIO.write not called if exit)
        results_no_seqs = self._run_main_with_mocked_parse(
            "no_seqs_output.fasta", parse_fasta_return_value=[]
        )
        results_no_seqs["mock_seqio_write"].assert_not_called()

    # --- 6. Error Handling (already covered FileNotFoundError, unequal lengths, no sequences) ---
    # No other specific exceptions seem to be explicitly handled by nogaps.py beyond sys.exit calls.


if __name__ == "__main__":
    unittest.main(argv=["first-arg-is-ignored"], exit=False, verbosity=2) # Ensure single newline at EOF
