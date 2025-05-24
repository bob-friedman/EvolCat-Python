import unittest
from unittest.mock import patch, MagicMock, mock_open
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
    from pylib.scripts import find_blast_top_pairs
except ModuleNotFoundError:
    scripts_dir = os.path.join(pylib_root, "pylib", "scripts")
    if scripts_dir not in sys.path:
        sys.path.insert(0, scripts_dir)
    import find_blast_top_pairs


class TestFindBlastTopPairs(unittest.TestCase):

    def setUp(self):
        # Common setup for tests if needed
        pass

    def _run_main_with_args(self, arg_list, stdin_data=None, mock_file_data=None):
        """Helper function to run the script's main() with specified args and stdin/file."""
        full_argv = ['find_blast_top_pairs.py'] + arg_list
        
        mock_stdin = io.StringIO(stdin_data if stdin_data is not None else "")
        
        # Determine if input_file is in args to decide whether to mock open
        input_file_path = None
        if '--input_file' in arg_list:
            try:
                input_file_path = arg_list[arg_list.index('--input_file') + 1]
            except IndexError: # Handle case where --input_file is last arg
                pass
        
        open_mock = mock_open(read_data=mock_file_data if mock_file_data is not None else "")
        
        # Patch open only if input_file is specified, otherwise patch stdin
        # The script logic is: if args.input_file, use open(), else use sys.stdin.
        # So, we always patch sys.stdin, and conditionally patch open based on args.
        
        with patch.object(sys, 'argv', full_argv), \
             patch('sys.stdin', mock_stdin), \
             patch('builtins.open', open_mock) if input_file_path else MagicMock(), \
             redirect_stdout(io.StringIO()) as captured_stdout, \
             redirect_stderr(io.StringIO()) as captured_stderr, \
             patch('sys.exit') as mock_sys_exit:
            
            try:
                find_blast_top_pairs.main()
            except SystemExit:
                pass 
            
            # If open was patched, assert it was called with the correct path
            if input_file_path and open_mock.called:
                 open_mock.assert_called_once_with(input_file_path, 'r')

            return captured_stdout.getvalue(), captured_stderr.getvalue(), mock_sys_exit


    # --- 1. Argument Parsing Tests (Covered in previous turn, re-verified) ---
    def test_arg_parser_defaults(self):
        with patch('pylib.scripts.find_blast_top_pairs.process_input_stream') as mock_process_stream:
            self._run_main_with_args([]) 
            self.assertTrue(mock_process_stream.called)
            passed_args = mock_process_stream.call_args[0][1]
            self.assertEqual(passed_args.query_id_col_idx, 0)
            self.assertEqual(passed_args.subject_id_col_idx, 1)
            self.assertEqual(passed_args.score_col_idx, 2)
            self.assertEqual(passed_args.filter_col_idx, 3)
            self.assertEqual(passed_args.filter_threshold, 0.0)
            self.assertIsNone(passed_args.input_file)

    def test_arg_parser_custom_values(self):
        custom_args_list = [
            '--input_file', 'test.blast', '--query_id_col_idx', '1',
            '--subject_id_col_idx', '2', '--score_col_idx', '10',
            '--filter_col_idx', '4', '--filter_threshold', '50.5'
        ]
        with patch('pylib.scripts.find_blast_top_pairs.process_input_stream') as mock_process_stream:
            self._run_main_with_args(custom_args_list, mock_file_data="") # Provide data for mocked open
            self.assertTrue(mock_process_stream.called)
            passed_args = mock_process_stream.call_args[0][1]
            self.assertEqual(passed_args.input_file, 'test.blast')
            self.assertEqual(passed_args.query_id_col_idx, 1); self.assertEqual(passed_args.subject_id_col_idx, 2)
            self.assertEqual(passed_args.score_col_idx, 10); self.assertEqual(passed_args.filter_col_idx, 4)
            self.assertEqual(passed_args.filter_threshold, 50.5)

    # --- 2. Input Handling ---
    def test_input_from_stdin(self):
        stdin_data = "q1\ts1\t100\t10\nextra\nq2\ts2\t90\t20\n"
        stdout, _, _ = self._run_main_with_args([], stdin_data=stdin_data)
        self.assertIn("q1\ts1\t100\t10", stdout)
        self.assertIn("q2\ts2\t90\t20", stdout)

    def test_input_from_file(self):
        file_data = "q_file\ts_file\t150\t25\n"
        stdout, _, _ = self._run_main_with_args(['--input_file', 'data.txt'], mock_file_data=file_data)
        self.assertIn("q_file\ts_file\t150\t25", stdout)

    def test_input_file_not_found_error(self):
        # Patch open to raise FileNotFoundError
        with patch('builtins.open', side_effect=FileNotFoundError("No such file")) as mock_bad_open:
            _, stderr, mock_exit = self._run_main_with_args(['--input_file', 'non_existent.txt'])
        
        mock_bad_open.assert_called_once_with('non_existent.txt', 'r')
        self.assertIn("Error: Input file 'non_existent.txt' not found.", stderr)
        mock_exit.assert_called_once_with(1)

    # --- 3. Line Parsing and Data Conversion ---
    def test_line_parsing_correct_indices(self):
        # q_idx=1, s_idx=0, score_idx=3, filter_idx=2. Threshold 5.
        # Input: s1 q1 10 100
        # Output: q1 s1 100 10 (first 4 fields of original, but logic uses score from score_idx)
        data = "s1\tq1\t10\t100.0\n" # subj, query, filter_val, score_val
        args = ['--query_id_col_idx', '1', '--subject_id_col_idx', '0', 
                '--score_col_idx', '3', '--filter_col_idx', '2', '--filter_threshold', '5.0']
        stdout, _, _ = self._run_main_with_args(args, stdin_data=data)
        self.assertIn("s1\tq1\t10\t100.0", stdout) # Original line printed

    def test_line_parsing_skip_fewer_columns(self):
        data = "q1\ts1\t100\n" # Missing filter_val if default indices (0,1,2,3)
        _, stderr, _ = self._run_main_with_args([], stdin_data=data)
        self.assertIn("Skipping line (not enough columns or conversion error)", stderr)
        self.assertIn("Original line: q1\ts1\t100", stderr)

    def test_line_parsing_skip_non_numeric_score(self):
        data = "q1\ts1\tNON_NUMERIC\t10\n"
        _, stderr, _ = self._run_main_with_args([], stdin_data=data)
        self.assertIn("Skipping line (not enough columns or conversion error)", stderr)
        self.assertIn("Original line: q1\ts1\tNON_NUMERIC\t10", stderr)

    def test_line_parsing_skip_non_numeric_filter_value(self):
        data = "q1\ts1\t100\tNON_NUMERIC_FILTER\n"
        _, stderr, _ = self._run_main_with_args([], stdin_data=data)
        self.assertIn("Skipping line (not enough columns or conversion error)", stderr)
        self.assertIn("Original line: q1\ts1\t100\tNON_NUMERIC_FILTER", stderr)

    def test_line_parsing_self_hits_excluded(self):
        data = "q1\tq1\t100\t10\n" # Self-hit
        stdout, _, _ = self._run_main_with_args([], stdin_data=data)
        self.assertEqual(stdout, "") # No output as self-hit is excluded

    # --- 4. Sorting Logic (implicitly tested by Top Pair Selection) ---
    # The core logic sorts by score. If top pair selection works, sorting is implicit.

    # --- 5. Top Pair Selection Logic ---
    def test_top_pair_simple_case(self):
        data = "q1\ts1\t100\t0\nq1\ts2\t90\t0\n" # q1 prefers s1
        stdout, _, _ = self._run_main_with_args([], stdin_data=data)
        self.assertIn("q1\ts1\t100\t0", stdout)
        self.assertNotIn("q1\ts2\t90\t0", stdout)

    def test_top_pair_subject_claimed(self):
        # s1 is hit by q1 (score 100) and q2 (score 90). q1 should claim s1.
        data = "q1\ts1\t100\t0\nq2\ts1\t90\t0\n"
        stdout, _, _ = self._run_main_with_args([], stdin_data=data)
        self.assertIn("q1\ts1\t100\t0", stdout)
        self.assertNotIn("q2\ts1\t90\t0", stdout)

    def test_top_pair_bidirectional_best(self):
        # q1-s1 is best for both q1 and s1.
        data = "q1\ts1\t100\t0\nq1\ts2\t50\t0\nq2\ts1\t60\t0\n"
        stdout, _, _ = self._run_main_with_args([], stdin_data=data)
        self.assertIn("q1\ts1\t100\t0", stdout)
        self.assertNotIn("q1\ts2\t50\t0", stdout)
        self.assertNotIn("q2\ts1\t60\t0", stdout)
    
    def test_top_pair_second_best_selected_after_conflict(self):
        # q1-s1 (100) -> q1 claims s1
        # q2-s1 (90)  -> q2 cannot claim s1 (claimed by q1)
        # q2-s2 (80)  -> q2 can claim s2
        # q3-s2 (70)  -> q3 cannot claim s2 (claimed by q2)
        data = "q1\ts1\t100\t0\nq2\ts1\t90\t0\nq2\ts2\t80\t0\nq3\ts2\t70\t0\n"
        stdout, _, _ = self._run_main_with_args([], stdin_data=data)
        self.assertIn("q1\ts1\t100\t0", stdout)
        self.assertIn("q2\ts2\t80\t0", stdout)
        self.assertNotIn("q2\ts1\t90\t0", stdout)
        self.assertNotIn("q3\ts2\t70\t0", stdout)

    def test_top_pair_order_of_input_does_not_matter_for_score(self):
        data1 = "q1\ts1\t100\t0\nq2\ts1\t90\t0\n" # q1-s1 is best for s1
        data2 = "q2\ts1\t90\t0\nq1\ts1\t100\t0\n" # Same, q1-s1 still best for s1
        
        stdout1, _, _ = self._run_main_with_args([], stdin_data=data1)
        stdout2, _, _ = self._run_main_with_args([], stdin_data=data2)
        
        self.assertIn("q1\ts1\t100\t0", stdout1)
        self.assertNotIn("q2\ts1\t90\t0", stdout1)
        self.assertIn("q1\ts1\t100\t0", stdout2)
        self.assertNotIn("q2\ts1\t90\t0", stdout2)
        self.assertEqual(stdout1, stdout2)

    # --- 6. Filtering Logic ---
    def test_filter_threshold_applied_correctly(self):
        # Default filter_col_idx=3, threshold=0.0
        # Line 1: filter_val=10 > 5.0 -> Should pass
        # Line 2: filter_val=5.0 == 5.0 -> Should pass
        # Line 3: filter_val=4.9 < 5.0 -> Should NOT pass
        data = "q1\ts1\t100\t10.0\nq2\ts2\t90\t5.0\nq3\ts3\t80\t4.9\n"
        args = ['--filter_threshold', '5.0']
        stdout, _, _ = self._run_main_with_args(args, stdin_data=data)
        self.assertIn("q1\ts1\t100\t10.0", stdout)
        self.assertIn("q2\ts2\t90\t5.0", stdout)
        self.assertNotIn("q3\ts3\t80\t4.9", stdout)

    def test_filter_threshold_default_zero(self):
        # Default threshold is 0.0.
        # Line 1: filter_val=0.1 > 0.0 -> Pass
        # Line 2: filter_val=0.0 == 0.0 -> Pass
        # Line 3: filter_val=-1.0 < 0.0 -> Fail
        data = "q1\ts1\t100\t0.1\nq2\ts2\t90\t0.0\nq3\ts3\t80\t-1.0\n"
        stdout, _, _ = self._run_main_with_args([], stdin_data=data)
        self.assertIn("q1\ts1\t100\t0.1", stdout)
        self.assertIn("q2\ts2\t90\t0.0", stdout)
        self.assertNotIn("q3\ts3\t80\t-1.0", stdout)

    # --- 7. Output Verification ---
    def test_output_tab_delimited_first_four_fields(self):
        data = "q1\ts1\t100.0\t60.0\textra_field1\textra_field2\n"
        stdout, _, _ = self._run_main_with_args([], stdin_data=data)
        # Script prints the first 4 fields of the original line
        expected_output = "q1\ts1\t100.0\t60.0\n" 
        self.assertEqual(stdout, expected_output)

    def test_output_fewer_than_four_fields_input(self):
        # If input has < 4 fields, but enough for configured indices and passes parsing.
        # Script tries to print fields[0:4]. If fewer, it prints what's available.
        data_2_fields = "q1\ts1\n" # Score/Filter will fail parsing, so this won't be output.
                                 # Let's make it parseable but short for output.
        data_3_fields_parseable = "q1\ts1\t50\n" # Score=50, Filter will fail parsing with default filter_idx=3
        
        # Test with custom indices to make data_3_fields_parseable actually parse and output
        # query=0, subject=1, score=2. No filter used (filter_threshold=0, filter_idx could be out of bounds but not used).
        # If filter_idx is out of bounds, it should warn.
        # Let's make filter_idx also 2, and threshold very low.
        args_custom_idx = ['--score_col_idx', '2', '--filter_col_idx', '2', '--filter_threshold', '-100']
        
        stdout_3_fields, stderr_3_fields, _ = self._run_main_with_args(args_custom_idx, stdin_data=data_3_fields_parseable)
        # The line "q1\ts1\t50" has 3 fields. fields[0:4] will give "q1","s1","50".
        # Output should be "q1\ts1\t50"
        self.assertIn("q1\ts1\t50", stdout_3_fields) # Original line gets printed.
        self.assertEqual(stdout_3_fields.strip().count('\t'), 2) # 2 tabs for 3 fields

        # Test if line parsing fails due to insufficient columns for specified indices
        data_short_for_indices = "q1\ts1\n" # query_idx=0, subject_idx=1. score_idx=2 (default) is missing.
        stdout_short, stderr_short, _ = self._run_main_with_args([], stdin_data=data_short_for_indices)
        self.assertEqual(stdout_short, "") # No output
        self.assertIn("Skipping line (not enough columns or conversion error)", stderr_short)


    # --- 8. Error Handling & Warnings (already covered some) ---
    def test_warning_for_malformed_line_non_numeric(self):
        data = "q1\ts1\tBAD_SCORE\t10\n"
        _, stderr, _ = self._run_main_with_args([], stdin_data=data)
        self.assertIn("Skipping line (not enough columns or conversion error): ValueError", stderr.upper()) # Check for ValueError part
        self.assertIn("Original line: q1\ts1\tBAD_SCORE\t10", stderr)

    def test_warning_for_malformed_line_too_few_cols_for_indices(self):
        # Default indices: Q=0, S=1, Score=2, Filter=3
        data = "q1\ts1\t100\n" # Missing column for default filter_idx=3
        _, stderr, _ = self._run_main_with_args([], stdin_data=data)
        self.assertIn("Skipping line (not enough columns or conversion error): IndexError", stderr.upper())
        self.assertIn("Original line: q1\ts1\t100", stderr)

if __name__ == '__main__':
    unittest.main(argv=['first-arg-is-ignored'], exit=False, verbosity=2)
