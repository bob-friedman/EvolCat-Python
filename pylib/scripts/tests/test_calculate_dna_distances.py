import unittest
from unittest.mock import patch, MagicMock, mock_open, call
import os
import sys
import io
from contextlib import redirect_stdout, redirect_stderr
import math

# Add paths for script and utility_functions import
current_dir = os.path.dirname(os.path.abspath(__file__))
# Ensure 'pylib' is in sys.path to allow 'from pylib.utils import seq_parser'
# and 'from pylib.scripts import calculate_dna_distances'
pylib_root = os.path.dirname(os.path.dirname(current_dir)) 
if pylib_root not in sys.path:
    sys.path.insert(0, pylib_root)

# Import the script to be tested
try:
    from pylib.scripts import calculate_dna_distances
except ModuleNotFoundError:
    scripts_dir = os.path.join(pylib_root, "pylib", "scripts")
    if scripts_dir not in sys.path:
        sys.path.insert(0, scripts_dir)
    import calculate_dna_distances

# Mock Biopython's SeqRecord for creating test sequences
class MockSeqRecord:
    def __init__(self, seq_str, id, name='', description=''):
        # Store seq as a string, as the script will call .upper() on it
        self.seq = seq_str 
        self.id = id
        self.name = name
        self.description = description

    def __len__(self):
        return len(self.seq)

    def __str__(self): 
        return str(self.seq)

    # Add upper method to mimic Biopython's Seq object behavior if script uses it
    def upper(self):
        return self.seq.upper()


class TestCalculateDnaDistances(unittest.TestCase):

    def setUp(self):
        # Reset relevant parts of the module if they maintain state.
        # For calculate_dna_distances, this seems less of an issue as functions are mostly pure.
        # However, if there were global configurations or caches, reset them here.
        pass

    # --- 1. Argument Parsing Tests ---
    def test_arg_parser_requires_input_file(self):
        with patch.object(sys, 'argv', ['calculate_dna_distances.py']):
            with self.assertRaises(SystemExit) as cm:
                with redirect_stderr(io.StringIO()): 
                    calculate_dna_distances.main()
            self.assertEqual(cm.exception.code, 2)

    @patch('pylib.utils.seq_parser.parse_fasta_file', side_effect=FileNotFoundError("Mocked FileNotFoundError"))
    @patch('sys.exit')
    def test_arg_parser_accepts_input_file_and_handles_not_found(self, mock_sys_exit, mock_parse_fasta):
        dummy_fasta_path = "non_existent.fasta"
        with patch.object(sys, 'argv', ['calculate_dna_distances.py', dummy_fasta_path]):
            with redirect_stderr(io.StringIO()) as stderr_capture:
                calculate_dna_distances.main()
        
        mock_parse_fasta.assert_called_once_with(dummy_fasta_path)
        self.assertIn(f"Error: File not found: {dummy_fasta_path}", stderr_capture.getvalue())
        mock_sys_exit.assert_called_once_with(1)


    # --- 2. FASTA Parsing and Validation ---
    @patch('pylib.utils.seq_parser.parse_fasta_file')
    @patch('sys.exit') 
    def test_sequences_different_lengths_error(self, mock_sys_exit, mock_parse_fasta):
        mock_parse_fasta.return_value = [
            MockSeqRecord("ACGT", id="seq1"),
            MockSeqRecord("ACGTT", id="seq2") 
        ]
        with patch.object(sys, 'argv', ['calculate_dna_distances.py', 'dummy.fasta']):
            with redirect_stderr(io.StringIO()) as stderr_capture:
                calculate_dna_distances.main()
        
        mock_sys_exit.assert_called_once_with(1)
        self.assertIn("Error: Sequences must be of the same length for distance calculation.", stderr_capture.getvalue())

    @patch('pylib.utils.seq_parser.parse_fasta_file')
    @patch('sys.exit')
    def test_fewer_than_two_sequences_error(self, mock_sys_exit, mock_parse_fasta):
        mock_parse_fasta.return_value = [
            MockSeqRecord("ACGT", id="seq1") 
        ]
        with patch.object(sys, 'argv', ['calculate_dna_distances.py', 'dummy.fasta']):
            with redirect_stderr(io.StringIO()) as stderr_capture:
                calculate_dna_distances.main()
        
        mock_sys_exit.assert_called_once_with(1)
        self.assertIn("Error: Need at least two sequences to calculate distances.", stderr_capture.getvalue())

    @patch('pylib.utils.seq_parser.parse_fasta_file')
    @patch('sys.exit') # Mock exit to inspect output
    def test_empty_sequences_prints_na_distances(self, mock_sys_exit, mock_parse_fasta):
        mock_parse_fasta.return_value = [
            MockSeqRecord("", id="seq1"), 
            MockSeqRecord("", id="seq2")
        ]
        with patch.object(sys, 'argv', ['calculate_dna_distances.py', 'dummy.fasta']):
            with redirect_stdout(io.StringIO()) as stdout_capture:
                calculate_dna_distances.main()
        
        output = stdout_capture.getvalue()
        self.assertIn("Pair\tRaw_Diffs\tValid_Sites\tProp_Diff\tJC69_d\tJC69_SE", output) # Check header
        self.assertIn("seq1 vs seq2\t0\t0\tN/A\tN/A\tN/A", output) # Check data line for N/A
        # Script should exit cleanly, not with error code 1 if it prints N/A
        # Depending on how sys.exit is called (or not called for success)
        if mock_sys_exit.called:
            self.assertNotEqual(mock_sys_exit.call_args[0][0], 1)


    @patch('pylib.utils.seq_parser.parse_fasta_file')
    def test_sequences_with_invalid_chars_ignored(self, mock_parse_fasta):
        # Invalid chars like N, -, ? should be ignored for pair counts, affecting valid_sites
        mock_parse_fasta.return_value = [
            MockSeqRecord("ACGTN-?", id="seq1"),
            MockSeqRecord("ACGTN-?", id="seq2"), # Identical, 0 diffs, 4 valid sites
            MockSeqRecord("NNNNAAAA", id="seq3"),
            MockSeqRecord("----CCCC", id="seq4")  # NNNN vs ---- (0 valid), AAAA vs CCCC (4 diffs, 4 valid)
        ]
        with patch.object(sys, 'argv', ['calculate_dna_distances.py', 'dummy.fasta']):
            with redirect_stdout(io.StringIO()) as stdout_capture:
                calculate_dna_distances.main()
        
        output = stdout_capture.getvalue()
        self.assertIn("seq1 vs seq2\t0\t4\t0.0000", output) 
        # For seq3 vs seq4:
        # NNNN vs ---- : 0 valid sites
        # AAAA vs CCCC : 4 diffs, 4 valid sites
        # Total: 4 diffs, 4 valid sites. Prop_Diff = 1.0
        # JC69 for P=1.0 is Inf.
        self.assertIn("seq3 vs seq4\t4\t4\t1.0000\tInf\tN/A", output)


    # --- 3. Individual Distance Function Tests ---
    
    # JC69
    def test_calculate_jc69_defined(self):
        P_prop, Q_prop, valid_sites = 0.1, 0.0, 100
        d, se = calculate_dna_distances.calculate_jc69(P_prop, Q_prop, valid_sites)
        self.assertAlmostEqual(d, 0.1073, places=4)
        self.assertAlmostEqual(se, 0.0346, places=4)

    def test_calculate_jc69_undefined_log(self):
        d, se = calculate_dna_distances.calculate_jc69(0.75, 0.0, 100)
        self.assertEqual(d, "Inf")
        self.assertEqual(se, "N/A")

    def test_calculate_jc69_no_valid_sites(self):
        d, se = calculate_dna_distances.calculate_jc69(0.1, 0.0, 0)
        self.assertEqual(d, "N/A")
        self.assertEqual(se, "N/A")

    # K2P
    def test_calculate_k2p_defined(self):
        P_prop, Q_prop, valid_sites = 0.1, 0.05, 100 # P=transitions, Q=transversions
        d, se, R = calculate_dna_distances.calculate_k2p(P_prop, Q_prop, valid_sites)
        self.assertAlmostEqual(d, 0.1702, places=4)
        self.assertAlmostEqual(R, 2.0, places=4)
        self.assertAlmostEqual(se, 0.0497, places=4)

    def test_calculate_k2p_undefined_log(self):
        d1, _, _ = calculate_dna_distances.calculate_k2p(0.4, 0.3, 100) # 1-2P-Q <= 0
        self.assertEqual(d1, "Inf")
        d2, _, _ = calculate_dna_distances.calculate_k2p(0.1, 0.5, 100) # 1-2Q <= 0
        self.assertEqual(d2, "Inf")

    def test_calculate_k2p_ts_tv_ratio_edge_cases(self):
        _, _, R_inf = calculate_dna_distances.calculate_k2p(0.1, 0.0, 100) # Q=0, P>0
        self.assertEqual(R_inf, "Inf")
        _, _, R_na = calculate_dna_distances.calculate_k2p(0.0, 0.0, 100) # P=0, Q=0
        self.assertEqual(R_na, "N/A")
        _, _, R_zero = calculate_dna_distances.calculate_k2p(0.0, 0.1, 100) # P=0, Q>0
        self.assertAlmostEqual(R_zero, 0.0, places=4)

    def test_calculate_k2p_no_valid_sites(self):
        d, se, R = calculate_dna_distances.calculate_k2p(0.1, 0.05, 0)
        self.assertEqual(d, "N/A")
        self.assertEqual(se, "N/A")
        self.assertEqual(R, "N/A")

    # Tajima-Nei (simplified testing due to complexity)
    def test_calculate_tajima_nei_defined_symmetric_freqs(self):
        pair_counts_diff = {'AT': 1, 'GC': 1} # Transversions
        base_counts_all_seqs = {'A':1, 'C':1, 'G':1, 'T':1} # Gives equal freqs (0.25)
        valid_sites = 2
        P_prop_overall_diffs = 1.0 # All sites different
        # P_prop_transitions and Q_prop_transversions are for info print, not directly used in d_TN formula in script
        d, se = calculate_dna_distances.calculate_tajima_nei(
            pair_counts_diff, base_counts_all_seqs, valid_sites, 
            P_prop=0.0, Q_prop=1.0 # P_prop=transitions/valid, Q_prop=transversions/valid
        )
        self.assertAlmostEqual(d, 1.2978, places=4) # Value from previous manual calculation
        self.assertEqual(se, "N/A")

    def test_calculate_tajima_nei_h_zero(self):
        # h_tajima becomes zero if a required base for a pair type has zero frequency
        pair_counts_diff = {'AT': 1}
        base_counts_all_seqs = {'A':2, 'C':0, 'G':0, 'T':0} # Freq_T = 0
        valid_sites = 1
        d, se = calculate_dna_distances.calculate_tajima_nei(
            pair_counts_diff, base_counts_all_seqs, valid_sites, P_prop=0.0, Q_prop=1.0
        )
        self.assertEqual(d, "Inf")
        self.assertEqual(se, "N/A")
        
    def test_calculate_tajima_nei_no_valid_sites(self):
        d, se = calculate_dna_distances.calculate_tajima_nei({}, {}, 0, 0, 0)
        self.assertEqual(d, "N/A")
        self.assertEqual(se, "N/A")

    # Tamura-Nei (simplified testing)
    def test_calculate_tamura_nei_defined_example(self):
        # Example: seq1=AG, seq2=CT. P1=0, P2=0, Q=2/2=1.0
        # base_counts_diff_pair: A=1, C=1, G=1, T=1. Freqs = 0.25 each.
        # freq_purines_prop = 0.5, freq_pyrimidines_prop = 0.5
        # Term1 (P1): log(1 - 0/0.5 - 1.0/(2*0.5)) = log(1 - 0 - 1) = log(0) -> Inf
        # This implies that if any log argument is <=0, result is Inf.
        d, se = calculate_dna_distances.calculate_tamura_nei(
            p1_transitions_pur_prop=0.0, 
            p2_transitions_pyr_prop=0.0, 
            q_transversions_prop=1.0,
            base_counts_diff_pair={'A':1, 'C':1, 'G':1, 'T':1}, 
            base_counts_ident_pair={},
            valid_sites_pair=2
        )
        self.assertEqual(d, "Inf") # Because Q=1 and equal base freqs makes log args zero
        self.assertEqual(se, "N/A")

    def test_calculate_tamura_nei_no_diffs(self):
        # If P1=P2=Q=0, d should be 0.
        d, se = calculate_dna_distances.calculate_tamura_nei(
            0,0,0, {}, {'A':4}, 2 # Two identical sequences AA vs AA
        )
        self.assertAlmostEqual(d, 0.0, places=4)
        self.assertEqual(se, "N/A") # SE is N/A from script

    def test_calculate_tamura_nei_no_valid_sites(self):
        d, se = calculate_dna_distances.calculate_tamura_nei(0,0,0,{},{},0)
        self.assertEqual(d, "N/A")
        self.assertEqual(se, "N/A")

    # --- 4. Main Loop Logic ---
    @patch('pylib.utils.seq_parser.parse_fasta_file')
    def test_main_loop_p_q_aggregation(self, mock_parse_fasta):
        # Test how P (transitions) and Q (transversions) are counted.
        # Seq1: A G C T
        # Seq2: A A C C (Ts: G->A, Tv: T->C)
        mock_parse_fasta.return_value = [
            MockSeqRecord("AGCT", id="s1"),
            MockSeqRecord("AACC", id="s2")
        ]
        with patch.object(sys, 'argv', ['calculate_dna_distances.py', 'dummy.fasta']):
            with redirect_stdout(io.StringIO()) as stdout_capture:
                calculate_dna_distances.main()
        output = stdout_capture.getvalue()
        # Expected: 2 diffs, 4 valid sites. Prop_Diff = 0.5
        # G->A is transition (P=1). T->C is transition (P=1). Total P = 2.
        # No transversions (Q=0).
        # P_prop = 2/4 = 0.5. Q_prop = 0/4 = 0.0.
        # Check output line for s1 vs s2:
        # Raw_Diffs=2, Valid_Sites=4, Prop_Diff=0.5000
        # JC69 for P=0.5, Q=0: -0.75*ln(1-4/3*0.5) = -0.75*ln(1-0.6666) = -0.75*ln(0.3333) = -0.75*(-1.0986) = 0.824
        # K2P for P_ts=0.5, Q_tv=0: 
        #   d = -0.5*ln(1-2*0.5-0) -0.25*ln(1-2*0) = -0.5*ln(0) -0.25*ln(1) = Inf
        #   R = P/Q -> 0.5/0 -> Inf
        # Tajima-Nei: Needs base_counts_all_seqs. Here: A=3,C=3,G=1,T=1. Total=8.
        #   fA=3/8, fC=3/8, fG=1/8, fT=1/8.
        #   pair_counts_diff: {'GA':1, 'TC':1}
        # Tamura-Nei: P1(AG,GA)=1, P2(CT,TC)=1. Total P=2.
        #   p1_prop = 1/4=0.25, p2_prop = 1/4=0.25, q_prop = 0
        
        # Focus on what's easily verifiable in output: Raw_Diffs, Valid_Sites, Prop_Diff
        # And the P_prop, Q_prop that would be passed to the functions.
        # The printout includes P_prop (transitions) and Q_prop (transversions)
        # Header: ... Prop_Diff P_prop_Ts Q_prop_Tv JC69_d ... K2P_d K2P_SE K2P_R ...
        # Data: s1 vs s2  2  4  0.5000  0.5000  0.0000  0.8240  N/A  Inf  N/A  Inf ...
        self.assertIn("s1 vs s2\t2\t4\t0.5000\t0.5000\t0.0000", output) # Raw, Valid, Prop, P_ts, Q_tv
        self.assertIn("0.8240\tN/A", output) # JC69 d, SE (SE for P=0.5 is N/A by formula (1-4/3P)^2 in denom)
        self.assertIn("Inf\tN/A\tInf", output) # K2P d, SE, R

    @patch('pylib.utils.seq_parser.parse_fasta_file')
    def test_main_loop_all_gaps_no_valid_sites(self, mock_parse_fasta):
        mock_parse_fasta.return_value = [
            MockSeqRecord("----", id="s1"),
            MockSeqRecord("NNNN", id="s2")
        ]
        with patch.object(sys, 'argv', ['calculate_dna_distances.py', 'dummy.fasta']):
            with redirect_stdout(io.StringIO()) as stdout_capture:
                calculate_dna_distances.main()
        output = stdout_capture.getvalue()
        # Expected: 0 diffs, 0 valid sites. All N/A.
        self.assertIn("s1 vs s2\t0\t0\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A", output)

    # --- 5. Output Verification ---
    @patch('pylib.utils.seq_parser.parse_fasta_file')
    def test_output_formatting_precision_and_na_inf(self, mock_parse_fasta):
        mock_parse_fasta.return_value = [
            MockSeqRecord("AAAA", id="s1"), # P=0, Q=0
            MockSeqRecord("AAAA", id="s2"),
            MockSeqRecord("TTTT", id="s3")  # s1 vs s3: P=0, Q=1 (all diffs are transversions)
        ]
        with patch.object(sys, 'argv', ['calculate_dna_distances.py', 'dummy.fasta']):
            with redirect_stdout(io.StringIO()) as stdout_capture:
                calculate_dna_distances.main()
        output = stdout_capture.getvalue()
        
        # s1 vs s2: 0 diffs, 4 valid. All distances 0.0000. Ratios N/A.
        self.assertIn("s1 vs s2\t0\t4\t0.0000\t0.0000\t0.0000\t0.0000\t0.0000\t0.0000\t0.0000\tN/A", output)
        
        # s1 vs s3: 4 diffs, 4 valid. Prop_diff=1.0. P_ts=0, Q_tv=1.0
        # JC69 for P_overall=1.0 is Inf.
        # K2P for P_ts=0, Q_tv=1.0: d = -0.5*ln(1-0-1) -0.25*ln(1-2*1) = -0.5*ln(0) -0.25*ln(-1) = Inf
        # K2P R = P_ts/Q_tv = 0/1 = 0.0000
        self.assertIn("s1 vs s3\t4\t4\t1.0000\t0.0000\t1.0000\tInf\tN/A\tInf\tN/A\t0.0000", output)

    def test_output_header_presence_and_content(self):
        # Test with minimal valid sequences to just check the header
        with patch('pylib.utils.seq_parser.parse_fasta_file') as mock_parse, \
             patch.object(sys, 'argv', ['calculate_dna_distances.py', 'dummy.fasta']), \
             redirect_stdout(io.StringIO()) as stdout_capture:
            mock_parse.return_value = [MockSeqRecord("A", "s1"), MockSeqRecord("C", "s2")]
            calculate_dna_distances.main()
        
        output = stdout_capture.getvalue()
        expected_header_parts = [
            "Pair", "Raw_Diffs", "Valid_Sites", "Prop_Diff", 
            "P_prop_Ts", "Q_prop_Tv", 
            "JC69_d", "JC69_SE", 
            "K2P_d", "K2P_SE", "K2P_R",
            "TN93_d", "TN93_SE", # Tajima-Nei
            "TrNei_d", "TrNei_SE" # Tamura-Nei (TrN)
        ]
        # Check if all parts are in the first line of output
        first_line = output.splitlines()[0]
        for part in expected_header_parts:
            self.assertIn(part, first_line)

    # --- 6. Error Handling (already covered FileNotFoundError) ---
    # Other specific exceptions are not explicitly handled by the script apart from sys.exit calls.

if __name__ == '__main__':
    unittest.main(argv=['first-arg-is-ignored'], exit=False, verbosity=2)
