import unittest
from unittest.mock import patch, MagicMock, mock_open, call
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
    from pylib.scripts import dot_plot
except ModuleNotFoundError:
    scripts_dir = os.path.join(pylib_root, "pylib", "scripts")
    if scripts_dir not in sys.path:
        sys.path.insert(0, scripts_dir)
    import dot_plot

# Mock Biopython's SeqRecord for creating test sequences
class MockSeqRecord:
    def __init__(self, seq_str, id="test_id", name='', description=''):
        self.seq_str = seq_str 
        self.id = id
        self.name = name if name else id # Ensure name is populated like SeqRecord
        self.description = description
        # For simplicity in tests, make self.seq point to an object that behaves like Bio.Seq.Seq
        # The script uses str(record.seq).upper() or record.seq.upper()
        # and record.seq.reverse_complement()
        # The patched Bio.Seq.Seq will handle this.
        # self.seq = dot_plot.Seq(seq_str) # This would use the actual Bio.Seq if not for setUp patch
        self.seq = MockBioSeq(seq_str) # Use our MockBioSeq directly for clarity in record

    def __str__(self): 
        return self.seq_str
    
    def __len__(self): # Script uses len(seq_obj)
        return len(self.seq_str)


# Mock for Bio.Seq.Seq object
class MockBioSeq:
    def __init__(self, data):
        self._data = str(data).upper()

    def __str__(self):
        return self._data

    def __len__(self):
        return len(self._data)
    
    def upper(self): # Method called by script
        return self # Return self as it's already upper and behaves like a Seq obj
    
    def toString(self): # if script uses this (older biopython)
        return self._data

    def reverse_complement(self): # Method called by script
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', '-':'-'}
        rev_comp_str = "".join(complement.get(base, base) for base in reversed(self._data))
        return MockBioSeq(rev_comp_str)


class TestDotPlot(unittest.TestCase):

    def setUp(self):
        # Patch Bio.Seq.Seq to return our MockBioSeq
        # This means whenever `Seq(...)` is called in the script, it returns a MockBioSeq instance.
        self.seq_patcher = patch('Bio.Seq.Seq', side_effect=lambda data, *args, **kwargs: MockBioSeq(data))
        self.mock_bio_seq_class = self.seq_patcher.start()
        self.addCleanup(self.seq_patcher.stop)
        
        # Common args for many tests
        self.common_args_dict = {
            'seqfile1': 's1.fasta', 'seqfile2': 's2.fasta',
            'wordlen': 3, 'step': 1,
            'outfile': 'plot.png', 'dotfile': 'dots.txt',
            'title': 'Test Plot', 'reverse_complement2': False,
            'min_run_length': 1 # Default in DotPlotter class if not from args
        }

    def _get_default_args(self, **kwargs):
        """Helper to get a Namespace object with default args, overridden by kwargs."""
        args_dict = self.common_args_dict.copy()
        args_dict.update(kwargs)
        return unittest.mock.Namespace(**args_dict)


    # --- 1. Argument Parsing Tests (Continued) ---
    def test_arg_parser_all_required_args_present(self):
        args_list = [
            'dot_plot.py',
            '--seqfile1', 's1.fasta', '--seqfile2', 's2.fasta',
            '--wordlen', '3', '--step', '1',
            '--outfile', 'plot.png', '--dotfile', 'dots.txt'
        ]
        # If this doesn't raise SystemExit, argparse is happy
        try:
            with patch.object(sys, 'argv', args_list):
                 # Need to mock downstream to prevent full run
                with patch('pylib.scripts.dot_plot.DotPlotter'):
                    dot_plot.main()
        except SystemExit as e:
            self.fail(f"Argparse failed with required args present: {e}")


    def test_arg_parser_default_title(self):
        args_list_no_title = [
            'dot_plot.py', '--seqfile1', 's1.fasta', '--seqfile2', 's2.fasta',
            '--wordlen', '3', '--step', '1', '--outfile', 'p.png', '--dotfile', 'd.txt'
        ]
        parsed_args = dot_plot.parse_arguments(args_list_no_title[1:])
        self.assertEqual(parsed_args.title, "Dot Plot")

    def test_arg_parser_custom_title(self):
        custom_title = "My Custom Dot Plot"
        args_list_custom_title = [
            'dot_plot.py', '--seqfile1', 's1.fasta', '--seqfile2', 's2.fasta',
            '--wordlen', '3', '--step', '1', '--outfile', 'p.png', '--dotfile', 'd.txt',
            '--title', custom_title
        ]
        parsed_args = dot_plot.parse_arguments(args_list_custom_title[1:])
        self.assertEqual(parsed_args.title, custom_title)

    def test_arg_parser_reverse_complement2_default_false(self):
        args_list = ['dot_plot.py', '--seqfile1', 's1.fasta', '--seqfile2', 's2.fasta',
                     '--wordlen', '3', '--step', '1', '--outfile', 'p.png', '--dotfile', 'd.txt']
        parsed_args = dot_plot.parse_arguments(args_list[1:])
        self.assertFalse(parsed_args.reverse_complement2)

    def test_arg_parser_reverse_complement2_true(self):
        args_list = ['dot_plot.py', '--seqfile1', 's1.fasta', '--seqfile2', 's2.fasta',
                     '--wordlen', '3', '--step', '1', '--outfile', 'p.png', '--dotfile', 'd.txt',
                     '--reverse_complement2']
        parsed_args = dot_plot.parse_arguments(args_list[1:])
        self.assertTrue(parsed_args.reverse_complement2)

    # --- 2. Sequence Input Tests ---
    @patch('pylib.utils.seq_parser.parse_fasta_file')
    @patch('pylib.scripts.dot_plot.DotPlotter') # Mock the class
    @patch('sys.exit')
    def test_sequence_input_file_not_found(self, mock_exit, MockDotPlotter, mock_parse_fasta):
        # Test FileNotFoundError for seqfile1
        mock_parse_fasta.side_effect = FileNotFoundError("s1.fasta not found")
        args = self._get_default_args()
        with patch.object(sys, 'argv', ['dot_plot.py'] + [f'--{k}' if isinstance(v, bool) and v else f'--{k}={v}' for k,v in vars(args).items() if not (isinstance(v, bool) and not v)]):
             with redirect_stderr(io.StringIO()) as err:
                dot_plot.main()
        self.assertIn("Error: Sequence file 's1.fasta' not found.", err.getvalue())
        mock_exit.assert_called_with(1)

        # Test FileNotFoundError for seqfile2
        mock_exit.reset_mock()
        mock_parse_fasta.side_effect = [MockSeqRecord("ACGT", "s1"), FileNotFoundError("s2.fasta not found")]
        with patch.object(sys, 'argv', ['dot_plot.py'] + [f'--{k}' if isinstance(v, bool) and v else f'--{k}={v}' for k,v in vars(args).items() if not (isinstance(v, bool) and not v)]):
            with redirect_stderr(io.StringIO()) as err:
                dot_plot.main()
        self.assertIn("Error: Sequence file 's2.fasta' not found.", err.getvalue())
        mock_exit.assert_called_with(1)


    @patch('pylib.utils.seq_parser.parse_fasta_file')
    @patch('pylib.scripts.dot_plot.DotPlotter')
    @patch('sys.exit')
    def test_sequence_input_empty_fasta(self, mock_exit, MockDotPlotter, mock_parse_fasta):
        # Test empty FASTA for seqfile1
        mock_parse_fasta.return_value = [] # No records
        args = self._get_default_args()
        with patch.object(sys, 'argv', ['dot_plot.py'] + [f'--{k}' if isinstance(v, bool) and v else f'--{k}={v}' for k,v in vars(args).items() if not (isinstance(v, bool) and not v)]):
            with redirect_stderr(io.StringIO()) as err:
                dot_plot.main()
        self.assertIn("Error: No sequences found in 's1.fasta'.", err.getvalue())
        mock_exit.assert_called_with(1)

        # Test empty FASTA for seqfile2
        mock_exit.reset_mock()
        mock_parse_fasta.side_effect = [[MockSeqRecord("ACGT", "s1")], []]
        with patch.object(sys, 'argv', ['dot_plot.py'] + [f'--{k}' if isinstance(v, bool) and v else f'--{k}={v}' for k,v in vars(args).items() if not (isinstance(v, bool) and not v)]):
            with redirect_stderr(io.StringIO()) as err:
                dot_plot.main()
        self.assertIn("Error: No sequences found in 's2.fasta'.", err.getvalue())
        mock_exit.assert_called_with(1)

    @patch('pylib.utils.seq_parser.parse_fasta_file')
    @patch('pylib.scripts.dot_plot.DotPlotter')
    def test_sequence_input_case_conversion(self, MockDotPlotter, mock_parse_fasta):
        # DotPlotter instance will be created with these args
        mock_plotter_instance = MockDotPlotter.return_value
        mock_plotter_instance.generate_dot_file = MagicMock()
        mock_plotter_instance.read_dot_file_and_plot = MagicMock()

        seq1_lower = "acgt"
        seq2_mixed = "aAcCgGtT"
        mock_parse_fasta.side_effect = [
            [MockSeqRecord(seq1_lower, "s1")],
            [MockSeqRecord(seq2_mixed, "s2")]
        ]
        args = self._get_default_args()
        with patch.object(sys, 'argv', ['dot_plot.py'] + [f'--{k}' if isinstance(v, bool) and v else f'--{k}={v}' for k,v in vars(args).items() if not (isinstance(v, bool) and not v)]):
            dot_plot.main()
        
        # Check the sequences passed to DotPlotter constructor
        # DotPlotter(args, seq1_name, seq2_name, seq1_obj, seq2_obj)
        # The script does: seq1_obj = Seq(str(records1[0].seq).upper())
        # Our MockBioSeq gets uppercased data.
        constructor_args = MockDotPlotter.call_args[0]
        passed_seq1_obj = constructor_args[3]
        passed_seq2_obj = constructor_args[4]

        self.assertEqual(str(passed_seq1_obj), "ACGT")
        self.assertEqual(str(passed_seq2_obj), "AACCGGTT")


    # --- 3. Match Calculation and Dotfile Generation ---
    @patch('builtins.open', new_callable=mock_open)
    def test_generate_dot_file_no_matches(self, mock_file_open):
        args = self._get_default_args(wordlen=3)
        plotter = dot_plot.DotPlotter(args, "s1", "s2", MockBioSeq("AAA"), MockBioSeq("TTT"))
        plotter.generate_dot_file()
        mock_file_open.assert_called_once_with(args.dotfile, 'w')
        handle = mock_file_open()
        handle.write.assert_not_called() # No matches, no writes

    @patch('builtins.open', new_callable=mock_open)
    def test_generate_dot_file_simple_direct_repeats(self, mock_file_open):
        args = self._get_default_args(wordlen=3, step=1)
        # seq1 = "ABCDEFG", seq2 = "XYZABCXYZ" -> match "ABC" at (0,3)
        plotter = dot_plot.DotPlotter(args, "s1", "s2", MockBioSeq("ABCDEFG"), MockBioSeq("XYZABCXYZ"))
        plotter.generate_dot_file()
        handle = mock_file_open()
        # Expected: word "ABC" (len 3) in s1 at index 0 matches "ABC" in s2 at index 3
        handle.write.assert_any_call("0\t3\n") 

    @patch('builtins.open', new_callable=mock_open)
    def test_generate_dot_file_simple_inverted_repeat(self, mock_file_open):
        args = self._get_default_args(wordlen=3, step=1, reverse_complement2=True)
        # seq1 = "ABCDEFG", seq2_rc = "CBAKLM" (seq2 = "MLKCBA" -> "GHIJKLMCBA")
        # We want "ABC" in s1 to match "CBA" (rc of "GHI") in s2.
        # Let seq1 = "ABCXYZ", seq2 = "MLKABC" -> seq2_rc = "GCBTML" (if ABC maps to GCB)
        # No, seq2_rc needs to contain ABC. So seq2 contains GCB.
        # seq1 = "ABCXYZ", seq2 contains "GCB" -> seq2 = "MNOGCBDEF"
        # seq2_rc = "FEDAGCBTML" (ABC is found)
        plotter = dot_plot.DotPlotter(args, "s1", "s2", MockBioSeq("ABCXYZ"), MockBioSeq("MNOGCBDEF"))
        plotter.generate_dot_file()
        handle = mock_file_open()
        # "ABC" in s1 at 0. "ABC" in s2_rc at 3 (from "GCB" in original s2 at index 3, assuming GCB is reverse of ABC)
        # This depends on the MockBioSeq reverse_complement logic.
        # If seq2 = "XXXTCBXXX" (TCB is rc of GAT, not ABC)
        # If seq1 = "GATTACA", seq2 contains reverse of "TTA" (AAT) -> "TTA"
        # seq1 = "GATTACA", seq2 = "MYSEQTTAMAIN"
        # s2_rc = "NIAMTTKSEGYM" (if TTA in s2, AAT in s2_rc)
        # word "TTA" in s1 at 1.
        # If seq2_rc contains "TTA", then seq2 contains "AAT".
        # Let seq1 = "AAATTT", seq2 = "CCC AAAGGG" (AAAGGG rc is CCCTTT)
        # word "AAA" in s1 at 0. word "AAA" in s2_rc from "TTT" in s2.
        # If seq2 = "CCCTTTGGG", then seq2_rc = "CCCAAAGGG"
        # word "TTT" in s1 at 3. word "TTT" in s2_rc from "AAA" in s2.
        
        # Simpler: seq1 = "GAT", seq2 = "ATC" (rc of GAT)
        # word "GAT" (len 3) in s1 at 0.
        # s2_rc = "GAT". word "GAT" in s2_rc at 0.
        # Expected: 0 \t 0
        plotter_rc = dot_plot.DotPlotter(args, "s1rc", "s2rc", MockBioSeq("GAT"), MockBioSeq("ATC"))
        plotter_rc.generate_dot_file()
        handle.write.assert_any_call("0\t0\n")

    @patch('builtins.open', new_callable=mock_open)
    def test_generate_dot_file_wordlen_step(self, mock_file_open):
        args = self._get_default_args(wordlen=2, step=2)
        # seq1 = "ABABABA", seq2 = "CDABABEF"
        # words s1: AB(0), AB(2), AB(4)
        # words s2: CD(0), AB(2), AB(4)
        # Matches: s1_AB(0) vs s2_AB(2) -> 0,2
        #          s1_AB(0) vs s2_AB(4) -> 0,4
        #          s1_AB(2) vs s2_AB(2) -> 2,2
        #          s1_AB(2) vs s2_AB(4) -> 2,4
        #          s1_AB(4) vs s2_AB(2) -> 4,2
        #          s1_AB(4) vs s2_AB(4) -> 4,4
        plotter = dot_plot.DotPlotter(args, "s1", "s2", MockBioSeq("ABABABA"), MockBioSeq("CDABABEF"))
        plotter.generate_dot_file()
        handle = mock_file_open()
        expected_calls = [call("0\t2\n"), call("0\t4\n"), 
                          call("2\t2\n"), call("2\t4\n"),
                          call("4\t2\n"), call("4\t4\n")]
        handle.write.assert_has_calls(expected_calls, any_order=True)

    @patch('builtins.open', new_callable=mock_open)
    def test_generate_dot_file_skip_n_chars(self, mock_file_open):
        args = self._get_default_args(wordlen=3, step=1)
        # seq1 = "ABCNNNXYZ", seq2 = "DEFABCXYZ"
        # Word "ABC" in s1 (idx 0) matches s2 (idx 3) -> 0,3
        # Word "XYZ" in s1 (idx 6) matches s2 (idx 6) -> 6,6
        # Word "BCN" contains N, skipped. "CN_N", "N_NN" also skipped.
        plotter = dot_plot.DotPlotter(args, "s1", "s2", MockBioSeq("ABCNNNXYZ"), MockBioSeq("DEFABCXYZ"))
        plotter.generate_dot_file()
        handle = mock_file_open()
        handle.write.assert_any_call("0\t3\n")
        handle.write.assert_any_call("6\t6\n")
        # Check that words with N were not written (implicitly by not being called)
        # e.g. if BCN matched something, it shouldn't be written.
        # A more direct test would be to inspect the `words1` dict inside generate_dot_file,
        # but that's white-box. Black-box: ensure no unexpected calls.
        # Total calls should be 2 for the two valid matches.
        self.assertEqual(handle.write.call_count, 2)


    # --- 4. Plot Generation (Mocking Matplotlib) ---
    @patch('matplotlib.pyplot.savefig')
    @patch('matplotlib.pyplot.show') # if called
    @patch('matplotlib.pyplot.scatter')
    @patch('matplotlib.pyplot.subplots')
    @patch('builtins.open', new_callable=mock_open)
    def test_plot_generation_from_dotfile(self, mock_file_open_dot, mock_subplots, mock_scatter, mock_show, mock_savefig):
        mock_ax = MagicMock()
        mock_fig = MagicMock()
        mock_subplots.return_value = (mock_fig, mock_ax)
        
        # Mock dotfile content: "0\t0\n10\t20\n"
        mock_file_open_dot.return_value.read.return_value = "0\t0\n10\t20\n"

        args = self._get_default_args(title="Plot Test", dotfile="fake_dots.txt")
        seq1_len, seq2_len = 100, 120 # For axis limits
        seq1_name, seq2_name = "SeqA", "SeqB"

        plotter = dot_plot.DotPlotter(args, seq1_name, seq2_name, MockBioSeq("A"*seq1_len), MockBioSeq("C"*seq2_len))
        plotter.read_dot_file_and_plot() # This is called by main after generate_dot_file

        mock_file_open_dot.assert_called_once_with(args.dotfile, 'r')
        
        # Check scatter call
        # x_coords = [0, 10], y_coords = [0, 20]
        scatter_call_args = mock_scatter.call_args[0]
        self.assertListEqual(list(scatter_call_args[0]), [0, 10]) # x_coords
        self.assertListEqual(list(scatter_call_args[1]), [0, 20]) # y_coords
        
        # Check plot attributes
        mock_ax.set_title.assert_called_with("Plot Test: SeqA vs SeqB")
        mock_ax.set_xlim.assert_called_with(0, seq1_len)
        mock_ax.set_ylim.assert_called_with(0, seq2_len)
        # Check for text labels (seq names)
        # ax.text(x, y, text, ...)
        # Example: ax.text(0.5, 1.01, seq1_name, ...)
        # Example: ax.text(-0.01, 0.5, seq2_name, ..., rotation='vertical')
        # This requires checking the `call_args_list` for multiple calls to `ax.text`
        found_s1_label = any(c[0][2] == seq1_name for c in mock_ax.text.call_args_list)
        found_s2_label = any(c[0][2] == seq2_name and 'vertical' in c[1].get('rotation','') for c in mock_ax.text.call_args_list)
        self.assertTrue(found_s1_label, "Seq1 name label not found on plot")
        self.assertTrue(found_s2_label, "Seq2 name label (vertical) not found on plot")

        mock_savefig.assert_called_once_with(args.outfile)

    @patch('matplotlib.pyplot.scatter')
    @patch('matplotlib.pyplot.subplots')
    @patch('builtins.open', new_callable=mock_open)
    def test_plot_empty_or_invalid_dotfile(self, mock_file_open_dot, mock_subplots, mock_scatter):
        mock_ax = MagicMock()
        mock_subplots.return_value = (MagicMock(), mock_ax)

        # Case 1: Empty dotfile
        mock_file_open_dot.return_value.read.return_value = ""
        args = self._get_default_args(dotfile="empty_dots.txt")
        plotter_empty = dot_plot.DotPlotter(args, "s1", "s2", MockBioSeq("A"), MockBioSeq("C"))
        plotter_empty.read_dot_file_and_plot()
        mock_scatter.assert_not_called() # No data to scatter

        # Case 2: Invalid data in dotfile
        mock_scatter.reset_mock()
        mock_file_open_dot.return_value.read.return_value = "invalid\ndata\n1\tA\n" # Line 3 bad
        with redirect_stderr(io.StringIO()) as err:
            plotter_invalid = dot_plot.DotPlotter(args, "s1", "s2", MockBioSeq("A"*10), MockBioSeq("C"*10))
            plotter_invalid.read_dot_file_and_plot()
        
        # Scatter might be called for valid lines before error, or not at all if error handling is strict.
        # Script tries to convert to int, so "1\tA" will raise ValueError.
        # Check if scatter was called with only valid data, or if an error was logged.
        # The script has a try-except ValueError and prints a warning.
        self.assertIn("Skipping invalid line in dotfile", err.getvalue())
        # If there were any valid lines before "1\tA", scatter might have been called.
        # If "invalid\ndata" causes scatter to not be called with anything, or if it's empty:
        # self.assertFalse(mock_scatter.called or not mock_scatter.call_args[0][0])
        # This depends on how many valid lines it processes before hitting bad one.
        # If "invalid" and "data" are skipped, and "1\tA" is skipped, then no points.
        self.assertFalse(mock_scatter.called)


    # --- 5. Error Handling ---
    @patch('pylib.utils.seq_parser.parse_fasta_file')
    @patch('sys.exit')
    def test_error_sequences_too_short_for_wordlen(self, mock_exit, mock_parse_fasta):
        mock_parse_fasta.side_effect = [
            [MockSeqRecord("A", "s1")], # Length 1
            [MockSeqRecord("CG", "s2")] # Length 2
        ]
        args = self._get_default_args(wordlen=3) # Wordlen 3
        
        # Test s1 too short
        with patch.object(sys, 'argv', ['dot_plot.py'] + [f'--{k}' if isinstance(v, bool) and v else f'--{k}={v}' for k,v in vars(args).items() if not (isinstance(v, bool) and not v)]):
            with redirect_stderr(io.StringIO()) as err:
                dot_plot.main()
        self.assertIn("Error: Sequence s1 (length 1) is shorter than word length (3).", err.getvalue())
        mock_exit.assert_called_with(1)

        # Test s2 too short (after s1 is okay)
        mock_exit.reset_mock()
        mock_parse_fasta.side_effect = [
            [MockSeqRecord("ACGTT", "s1")], # Length 5
            [MockSeqRecord("CG", "s2")]    # Length 2
        ]
        with patch.object(sys, 'argv', ['dot_plot.py'] + [f'--{k}' if isinstance(v, bool) and v else f'--{k}={v}' for k,v in vars(args).items() if not (isinstance(v, bool) and not v)]):
            with redirect_stderr(io.StringIO()) as err:
                dot_plot.main()
        self.assertIn("Error: Sequence s2 (length 2) is shorter than word length (3).", err.getvalue())
        mock_exit.assert_called_with(1)


    # --- 6. End-to-End Flow (Successful Run Mocked) ---
    @patch('pylib.utils.seq_parser.parse_fasta_file')
    @patch('builtins.open', new_callable=mock_open)
    @patch('matplotlib.pyplot.subplots')
    @patch('matplotlib.pyplot.scatter')
    @patch('matplotlib.pyplot.savefig')
    @patch('matplotlib.pyplot.show')
    @patch('sys.exit') # To ensure it doesn't exit test suite
    def test_end_to_end_successful_run(self, mock_sys_exit, mock_show, mock_savefig, 
                                     mock_scatter, mock_subplots, mock_file_open, mock_parse_fasta):
        # Setup mocks
        mock_ax = MagicMock()
        mock_fig = MagicMock()
        mock_subplots.return_value = (mock_fig, mock_ax)
        
        mock_parse_fasta.side_effect = [
            [MockSeqRecord("GATTACA", "seqA")],
            [MockSeqRecord("CATGATT", "seqB")]
        ]
        # Mock dotfile writing and then reading
        # For writing: mock_file_open used by generate_dot_file
        # For reading: mock_file_open used by read_dot_file_and_plot
        # We need to make sure the "read" part gets some data.
        # Let's assume some matches are written.
        # GAT at 0 in sA. GAT at 3 in sB. -> "0\t3\n"
        # TTA at 3 in sA. TTA not in sB.
        # ACA at 4 in sA. ACA not in sB.
        # mock_file_open.return_value.write # This is already a mock from mock_open
        # For reading, we need to mock read()
        mock_dot_file_content = "0\t3\n"
        
        # Side effect for open: first call is write (dotfile), second is read (dotfile)
        mock_dot_write_handle = MagicMock()
        mock_dot_read_handle = MagicMock()
        mock_dot_read_handle.read.return_value = mock_dot_file_content
        mock_file_open.side_effect = [mock_dot_write_handle, mock_dot_read_handle]


        args = self._get_default_args(wordlen=3, step=1, title="E2E Test")
        args_list = ['dot_plot.py'] + [f'--{k}' if isinstance(v, bool) and v else f'--{k}={v}' for k,v in vars(args).items() if not (isinstance(v, bool) and not v)]
        
        with patch.object(sys, 'argv', args_list):
            with redirect_stdout(io.StringIO()) as out, redirect_stderr(io.StringIO()) as err:
                dot_plot.main()

        # Verify sequence parsing calls
        self.assertEqual(mock_parse_fasta.call_count, 2)
        mock_parse_fasta.assert_any_call(args.seqfile1)
        mock_parse_fasta.assert_any_call(args.seqfile2)

        # Verify dotfile write
        mock_file_open.assert_any_call(args.dotfile, 'w')
        mock_dot_write_handle.write.assert_any_call("0\t3\n") # Check if "GAT" match was written

        # Verify dotfile read
        mock_file_open.assert_any_call(args.dotfile, 'r')
        mock_dot_read_handle.read.assert_called_once()

        # Verify plotting calls
        mock_subplots.assert_called_once()
        mock_scatter.assert_called_once_with([0], [3], s=1, c='black', marker='.')
        mock_ax.set_title.assert_called_with("E2E Test: seqA vs seqB")
        mock_savefig.assert_called_once_with(args.outfile)
        
        # Check for success message (if any)
        self.assertIn(f"Dot file saved to: {args.dotfile}", out.getvalue())
        self.assertIn(f"Plot image saved to: {args.outfile}", out.getvalue())
        
        mock_sys_exit.assert_not_called() # Should complete without exiting due to error


if __name__ == '__main__':
    unittest.main(argv=['first-arg-is-ignored'], exit=False, verbosity=2)
