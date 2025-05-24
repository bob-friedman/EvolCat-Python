import unittest
from unittest.mock import patch, MagicMock, mock_open
import os
import sys
import io
from contextlib import redirect_stdout, redirect_stderr

# Add paths for script and utility_functions import
current_dir = os.path.dirname(os.path.abspath(__file__))
pylib_root = os.path.dirname(os.path.dirname(current_dir)) 
if pylib_root not in sys.path:
    sys.path.insert(0, pylib_root)

# Import the script to be tested
try:
    from pylib.scripts import find_tandem_repeats
except ModuleNotFoundError:
    scripts_dir = os.path.join(pylib_root, "pylib", "scripts")
    if scripts_dir not in sys.path:
        sys.path.insert(0, scripts_dir)
    import find_tandem_repeats

class MockSeqRecord:
    def __init__(self, seq_str, id="test_id", name='', description=''):
        self.seq = seq_str 
        self.id = id
        self.name = name
        self.description = description
    def __str__(self): 
        return str(self.seq)

class TestFindTandemRepeats(unittest.TestCase):

    def setUp(self):
        pass

    # --- 1. Argument Parsing Tests (Covered in previous turn, re-verified) ---
    def test_arg_parser_pattern_required(self):
        with patch.object(sys, 'argv', ['find_tandem_repeats.py', '--sequence', 'ACGT']):
            with self.assertRaises(SystemExit) as cm, redirect_stderr(io.StringIO()):
                find_tandem_repeats.main()
            self.assertEqual(cm.exception.code, 2)

    def test_arg_parser_sequence_source_required(self):
        with patch.object(sys, 'argv', ['find_tandem_repeats.py', '--pattern', 'AA']):
            with self.assertRaises(SystemExit) as cm, redirect_stderr(io.StringIO()) as r_err:
                find_tandem_repeats.main()
            self.assertEqual(cm.exception.code, 1)
            self.assertIn("Must provide either --sequence or --sequence_file", r_err.getvalue())


    def test_arg_parser_fudge_default(self):
        with patch.object(sys, 'argv', ['find_tandem_repeats.py', '--pattern', 'AA', '--sequence', 'AAAA']), \
             patch('pylib.scripts.find_tandem_repeats.process_sequence', return_value=None) as mock_process:
            find_tandem_repeats.main()
            args_used = mock_process.call_args[0][1]
            self.assertEqual(args_used.fudge, 0)

    def test_arg_parser_verbose_default_false(self):
        with patch.object(sys, 'argv', ['find_tandem_repeats.py', '--pattern', 'AA', '--sequence', 'AAAA']), \
             patch('pylib.scripts.find_tandem_repeats.process_sequence', return_value=None) as mock_process:
            find_tandem_repeats.main()
            args_used = mock_process.call_args[0][1]
            self.assertFalse(args_used.verbose)

    def test_arg_parser_verbose_true(self):
        with patch.object(sys, 'argv', ['find_tandem_repeats.py', '--pattern', 'AA', '--sequence', 'AAAA', '--verbose']), \
             patch('pylib.scripts.find_tandem_repeats.process_sequence', return_value=None) as mock_process:
            find_tandem_repeats.main()
            args_used = mock_process.call_args[0][1]
            self.assertTrue(args_used.verbose)

    # --- 2. Sequence Input Tests (Covered in previous turn, re-verified) ---
    @patch('pylib.scripts.find_tandem_repeats.process_sequence', return_value=None)
    def test_input_from_sequence_argument(self, mock_process_sequence):
        test_seq = "GATTACA"
        with patch.object(sys, 'argv', ['find_tandem_repeats.py', '--pattern', 'GA', '--sequence', test_seq]):
            find_tandem_repeats.main()
        call_args = mock_process_sequence.call_args[0]
        self.assertEqual(call_args[0], "command_line_sequence") 
        self.assertEqual(call_args[2], test_seq.upper())

    @patch('pylib.utils.seq_parser.parse_fasta_file')
    @patch('pylib.scripts.find_tandem_repeats.process_sequence', return_value=None)
    def test_input_from_fasta_file(self, mock_process_sequence, mock_parse_fasta):
        fasta_content_id = "test_seq_id"
        fasta_content_seq = "CGATCGAT"
        mock_parse_fasta.return_value = [MockSeqRecord(fasta_content_seq, id=fasta_content_id)]
        with patch.object(sys, 'argv', ['find_tandem_repeats.py', '--pattern', 'CG', '--sequence_file', 'dummy.fasta']):
            find_tandem_repeats.main()
        mock_parse_fasta.assert_called_once_with('dummy.fasta')
        call_args = mock_process_sequence.call_args[0]
        self.assertEqual(call_args[0], fasta_content_id)
        self.assertEqual(call_args[2], fasta_content_seq.upper())

    @patch('builtins.open', new_callable=mock_open, read_data="#comment\nactual_sequence_line\n  another line")
    @patch('pylib.scripts.find_tandem_repeats.process_sequence', return_value=None)
    def test_input_from_plain_text_file(self, mock_process_sequence, mock_file_open):
        file_path = "dummy.txt"
        with patch.object(sys, 'argv', ['find_tandem_repeats.py', '--pattern', 'SEQ', '--sequence_file', file_path]):
            find_tandem_repeats.main()
        mock_file_open.assert_called_once_with(file_path, 'r')
        call_args = mock_process_sequence.call_args[0]
        self.assertEqual(call_args[0], os.path.basename(file_path))
        self.assertEqual(call_args[2], "actual_sequence_line".upper().strip())

    @patch('pylib.utils.seq_parser.parse_fasta_file', side_effect=FileNotFoundError("File not found mock"))
    @patch('sys.exit')
    def test_input_sequence_file_not_found(self, mock_sys_exit, mock_parse_fasta):
        file_path = "non_existent.fasta"
        with patch.object(sys, 'argv', ['find_tandem_repeats.py', '--pattern', 'AA', '--sequence_file', file_path]), \
             redirect_stderr(io.StringIO()) as stderr_capture:
            find_tandem_repeats.main()
        self.assertIn(f"Error: Sequence file '{file_path}' not found.", stderr_capture.getvalue())
        mock_sys_exit.assert_called_once_with(1)

    # --- 3. Pattern Handling (Covered in previous turn, re-verified) ---
    def test_pattern_empty_error(self):
        with patch.object(sys, 'argv', ['find_tandem_repeats.py', '--pattern', '', '--sequence', 'ACGT']), \
             self.assertRaises(SystemExit) as cm, \
             redirect_stderr(io.StringIO()) as stderr_capture:
            find_tandem_repeats.main()
        self.assertEqual(cm.exception.code, 1)
        self.assertIn("Error: Pattern cannot be empty.", stderr_capture.getvalue())

    @patch('pylib.scripts.find_tandem_repeats.process_sequence', return_value=None)
    def test_pattern_and_sequence_case_insensitivity(self, mock_process_sequence):
        with patch.object(sys, 'argv', ['find_tandem_repeats.py', '--pattern', 'at', '--sequence', 'aTaTcGat']):
            find_tandem_repeats.main()
        args_passed_to_process = mock_process_sequence.call_args[0][1]
        sequence_content_passed = mock_process_sequence.call_args[0][2]
        self.assertEqual(args_passed_to_process.pattern, "AT")
        self.assertEqual(sequence_content_passed, "ATATCGAT")

    # --- 4. Repeat Finding Logic (Core Algorithm) & 5. Output Verification ---
    def _run_and_capture_output(self, args_list):
        with patch.object(sys, 'argv', args_list), \
             redirect_stdout(io.StringIO()) as stdout_capture, \
             redirect_stderr(io.StringIO()): # Suppress verbose if not testing for it specifically
            try:
                find_tandem_repeats.main()
            except SystemExit: # Catch sys.exit if it occurs, but we primarily check stdout
                pass 
            return stdout_capture.getvalue()

    def test_no_repeats_found(self):
        output = self._run_and_capture_output(['find_tandem_repeats.py', '--pattern', 'XXX', '--sequence', 'ACGTACGT'])
        self.assertIn("Sequence: command_line_sequence\tPattern: XXX\tFudge: 0", output)
        self.assertIn("Found 0 instances of pattern XXX", output)
        self.assertIn("Found 0 blocks of repeats", output)

    def test_simple_exact_repeats_aaa(self):
        # Pattern AAA in CGAAACG
        # Expected: Match at offset 2 (AAA). 1 instance. 1 block.
        # Block 1: 1 repeat of AAA, Start: 2, End: 4
        output = self._run_and_capture_output(['find_tandem_repeats.py', '--pattern', 'AAA', '--sequence', 'CGAAACG'])
        self.assertIn("Found 1 instances of pattern AAA", output)
        self.assertIn("Found 1 blocks of repeats", output)
        self.assertIn("Block 1: 1 repeats of AAA. Start: 2 End: 4", output)

    def test_simple_exact_repeats_at_at_at(self):
        # Pattern AT in CGATATATATGCG
        # AT at 2, AT at 4, AT at 6, AT at 8. Fudge 0.
        # Match 1: 'AT' at 2. Block starts. repeats_in_block=1. last_match_start=2.
        # Match 2: 'AT' at 4. Is 4 <= 2 + len(AT) + fudge(0)? 4 <= 2+2+0? Yes. repeats_in_block=2. last_match_start=4.
        # Match 3: 'AT' at 6. Is 6 <= 4 + 2 + 0? Yes. repeats_in_block=3. last_match_start=6.
        # Match 4: 'AT' at 8. Is 8 <= 6 + 2 + 0? Yes. repeats_in_block=4. last_match_start=8.
        # End of seq. Block 1: 4 repeats of AT. Start: 2, End: 8 + 2 - 1 = 9.
        output = self._run_and_capture_output(['find_tandem_repeats.py', '--pattern', 'AT', '--sequence', 'CGATATATATGCG'])
        self.assertIn("Found 4 instances of pattern AT", output)
        self.assertIn("Found 1 blocks of repeats", output)
        self.assertIn("Block 1: 4 repeats of AT. Start: 2 End: 9", output)

    def test_fudge_factor_joins_blocks(self):
        # Pattern XX, fudge 1: XX_XX (X X space X X) -> XXNXX
        # XX at 0. XX at 3.
        # Match 1: 'XX' at 0. Block starts. repeats=1. last_match_start=0.
        # Match 2: 'XX' at 3. Is 3 <= 0 + len(XX) + fudge(1)? 3 <= 0+2+1? Yes, 3 <= 3. repeats=2. last_match_start=3.
        # End. Block 1: 2 repeats. Start: 0, End: 3 + 2 - 1 = 4.
        output = self._run_and_capture_output(['find_tandem_repeats.py', '--pattern', 'XX', '--sequence', 'XXNXX', '--fudge', '1'])
        self.assertIn("Found 2 instances of pattern XX", output) # Instances are still distinct matches
        self.assertIn("Found 1 blocks of repeats", output)
        self.assertIn("Block 1: 2 repeats of XX. Start: 0 End: 4", output)

    def test_fudge_factor_separates_blocks(self):
        # Pattern XX, fudge 0: XX_XX
        # Match 1: 'XX' at 0. Block starts. repeats=1. last_match_start=0.
        # Match 2: 'XX' at 3. Is 3 <= 0 + 2 + 0? No, 3 > 2. New block.
        # Block 1: 1 repeat. Start 0, End 1.
        # Block 2: 1 repeat. Start 3, End 4.
        output = self._run_and_capture_output(['find_tandem_repeats.py', '--pattern', 'XX', '--sequence', 'XXNXX', '--fudge', '0'])
        self.assertIn("Found 2 instances of pattern XX", output)
        self.assertIn("Found 2 blocks of repeats", output)
        self.assertIn("Block 1: 1 repeats of XX. Start: 0 End: 1", output)
        self.assertIn("Block 2: 1 repeats of XX. Start: 3 End: 4", output)

    def test_fudge_factor_perfect_tandem(self):
        # Pattern XX, fudge 1: XXXX
        # XX at 0. XX at 1 (if pattern can overlap like this, script searches from match_start+1), XX at 2.
        # Script's search_start_idx = match_start_offset + 1
        # Pat XX, Seq XXXX
        # Match 1: XX at 0. last_match_start=0. repeats=1. search_idx=1.
        # Match 2: XX at 1. 1 <= 0+2+1=3. Yes. last_match_start=1. repeats=2. search_idx=2.
        # Match 3: XX at 2. 2 <= 1+2+1=4. Yes. last_match_start=2. repeats=3. search_idx=3.
        # End. Block 1: 3 repeats. Start 0, End 2+2-1 = 3.
        output = self._run_and_capture_output(['find_tandem_repeats.py', '--pattern', 'XX', '--sequence', 'XXXX', '--fudge', '1'])
        self.assertIn("Found 3 instances of pattern XX", output)
        self.assertIn("Found 1 blocks of repeats", output)
        self.assertIn("Block 1: 3 repeats of XX. Start: 0 End: 3", output)
        
    def test_fudge_factor_ta_various(self):
        # Pattern TA, fudge 1: TA N TA
        output1 = self._run_and_capture_output(['find_tandem_repeats.py', '--pattern', 'TA', '--sequence', 'TANTANNTA', '--fudge', '1'])
        # TA at 0. TA at 3. TA at 6.
        # M1: TA at 0. last=0, reps=1.
        # M2: TA at 3. 3 <= 0+2+1=3. Yes. last=3, reps=2.
        # M3: TA at 6. 6 <= 3+2+1=6. Yes. last=6, reps=3.
        # Block 1: 3 repeats. Start 0, End 6+2-1=7.
        self.assertIn("Block 1: 3 repeats of TA. Start: 0 End: 7", output1)

        # Pattern TA, fudge 2: TA N N TA
        output2 = self._run_and_capture_output(['find_tandem_repeats.py', '--pattern', 'TA', '--sequence', 'TANNTA', '--fudge', '2'])
        # TA at 0. TA at 4.
        # M1: TA at 0. last=0, reps=1.
        # M2: TA at 4. 4 <= 0+2+2=4. Yes. last=4, reps=2.
        # Block 1: 2 repeats. Start 0, End 4+2-1=5.
        self.assertIn("Block 1: 2 repeats of TA. Start: 0 End: 5", output2)

        # Pattern TA, fudge 1: TA N N TA (should be 2 blocks)
        output3 = self._run_and_capture_output(['find_tandem_repeats.py', '--pattern', 'TA', '--sequence', 'TANNTA', '--fudge', '1'])
        # M1: TA at 0. last=0, reps=1.
        # M2: TA at 4. 4 <= 0+2+1=3. No. New Block.
        # Block 1: 1 repeat. Start 0, End 1.
        # Block 2: 1 repeat. Start 4, End 5.
        self.assertIn("Found 2 blocks of repeats", output3)
        self.assertIn("Block 1: 1 repeats of TA. Start: 0 End: 1", output3)
        self.assertIn("Block 2: 1 repeats of TA. Start: 4 End: 5", output3)

    def test_overlapping_pattern_instances_aaaaa_aa(self):
        # Sequence AAAAA, pattern AA, fudge 0
        # M1: AA at 0. last=0, reps=1. search_idx=1.
        # M2: AA at 1. 1 <= 0+2+0=2. Yes. last=1, reps=2. search_idx=2.
        # M3: AA at 2. 2 <= 1+2+0=3. Yes. last=2, reps=3. search_idx=3.
        # M4: AA at 3. 3 <= 2+2+0=4. Yes. last=3, reps=4. search_idx=4.
        # Block 1: 4 repeats. Start 0, End 3+2-1=4.
        output = self._run_and_capture_output(['find_tandem_repeats.py', '--pattern', 'AA', '--sequence', 'AAAAA', '--fudge', '0'])
        self.assertIn("Found 4 instances of pattern AA", output)
        self.assertIn("Found 1 blocks of repeats", output)
        self.assertIn("Block 1: 4 repeats of AA. Start: 0 End: 4", output)

    def test_multiple_distinct_blocks(self):
        # Sequence: ATATATxxxxxAGAGAG, pattern AT (fudge 0), then pattern AG (fudge 0)
        seq = "ATATATXXXXXAGAGAG"
        output_at = self._run_and_capture_output(['find_tandem_repeats.py', '--pattern', 'AT', '--sequence', seq, '--fudge', '0'])
        # AT: M1 at 0, M2 at 2, M3 at 4. Block of 3. Start 0, End 4+2-1=5.
        self.assertIn("Block 1: 3 repeats of AT. Start: 0 End: 5", output_at)
        
        output_ag = self._run_and_capture_output(['find_tandem_repeats.py', '--pattern', 'AG', '--sequence', seq, '--fudge', '0'])
        # AG: M1 at 11, M2 at 13, M3 at 15. Block of 3. Start 11, End 15+2-1=16.
        self.assertIn("Block 1: 3 repeats of AG. Start: 11 End: 16", output_ag)

    def test_edge_cases_pattern_location(self):
        # Pattern at start
        output_start = self._run_and_capture_output(['find_tandem_repeats.py', '--pattern', 'XXX', '--sequence', 'XXXACGT', '--fudge', '0'])
        self.assertIn("Block 1: 1 repeats of XXX. Start: 0 End: 2", output_start)
        
        # Pattern at end
        output_end = self._run_and_capture_output(['find_tandem_repeats.py', '--pattern', 'XXX', '--sequence', 'ACGTXXX', '--fudge', '0'])
        self.assertIn("Block 1: 1 repeats of XXX. Start: 4 End: 6", output_end)

        # Sequence length equal to pattern length
        output_equal = self._run_and_capture_output(['find_tandem_repeats.py', '--pattern', 'ACGT', '--sequence', 'ACGT', '--fudge', '0'])
        self.assertIn("Block 1: 1 repeats of ACGT. Start: 0 End: 3", output_equal)

    def test_verbose_output(self):
        # Check for some key verbose messages.
        # Requires capturing stderr or patching print if verbose prints to stdout.
        # Script prints verbose messages to stdout.
        args = ['find_tandem_repeats.py', '--pattern', 'TA', '--sequence', 'TANTATA', '--fudge', '1', '--verbose']
        with patch.object(sys, 'argv', args), \
             redirect_stdout(io.StringIO()) as stdout_capture:
            find_tandem_repeats.main()
        
        output = stdout_capture.getvalue()
        self.assertIn("VERBOSE: Sequence: TANTATA Pattern: TA Fudge: 1", output)
        self.assertIn("VERBOSE: Match for TA found at offset 0", output)
        # TA at 0. TA at 3. TA at 5.
        # M1: TA at 0. last=0, reps=1.
        # M2: TA at 3. 3 <= 0+2+1=3. Yes. reps=2, last=3.
        # M3: TA at 5. 5 <= 3+2+1=6. Yes. reps=3, last=5.
        # Block 1: 3 repeats. Start 0, End 5+2-1=6.
        self.assertIn("VERBOSE: Continuing current block. Repeats in block: 2", output) # For match at offset 3
        self.assertIn("VERBOSE: Continuing current block. Repeats in block: 3", output) # For match at offset 5
        self.assertIn("VERBOSE: Ending current block. Repeats: 3, Start: 0, End: 6", output)
        self.assertIn("Block 1: 3 repeats of TA. Start: 0 End: 6", output)


if __name__ == '__main__':
    unittest.main(argv=['first-arg-is-ignored'], exit=False, verbosity=2)
