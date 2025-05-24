import unittest
from unittest.mock import patch, MagicMock, mock_open, call
import os
import sys
import io
from contextlib import redirect_stdout, redirect_stderr

# Add the directory containing the script to the Python path
# This ensures that 'blast_ncbi_tabular' can be imported
current_dir = os.path.dirname(os.path.abspath(__file__))
# project_root = os.path.dirname(os.path.dirname(os.path.dirname(current_dir))) # Assuming pylib is in project root
# sys.path.insert(0, project_root)
sys.path.insert(0, current_dir) # If script is in the same directory as test
sys.path.insert(0, os.path.dirname(current_dir)) # scripts directory
sys.path.insert(0, os.path.dirname(os.path.dirname(current_dir))) # pylib directory


# Import the script to be tested
try:
    from pylib.scripts.ncbi import blast_ncbi_tabular
except ModuleNotFoundError:
    # Fallback if running test directly and pylib structure isn't picked up
    import blast_ncbi_tabular


# Mock utility_functions early if it's imported at module level by the script
# Ensure this patch is active *before* blast_ncbi_tabular is fully imported and uses it,
# or patch it directly where it's used if that's cleaner.
# If blast_ncbi_tabular.utility_functions is the path, use that.
# The script uses: from utils import utility_functions
# This means we need to mock utils.utility_functions from the perspective of where blast_ncbi_tabular is.
# So, if blast_ncbi_tabular is in pylib.scripts.ncbi, then utils is pylib.utils
MOCK_USER_AGENT = "Test User Agent"

# Patching 'pylib.utils.utility_functions' if 'blast_ncbi_tabular' imports it as 'from utils import utility_functions'
# and 'utils' is a sibling of 'scripts' under 'pylib'.
# The script has a complex sys.path manipulation.
# Let's assume the script's import `from utils import utility_functions` works due to its sys.path changes.
# We need to patch `utils.utility_functions` or `blast_ncbi_tabular.utility_functions`
# depending on how it's resolved.
# The script itself tries to import `utils.utility_functions`.
# So we should patch that.
@patch('utils.utility_functions.get_user_agent', return_value=MOCK_USER_AGENT)
class TestBlastNcbiTabular(unittest.TestCase):

    def _mock_response(self, status_code=200, text_data="", content_data=b"", json_data=None, raise_for_status_exception=None):
        mock_resp = MagicMock(spec=requests.Response)
        mock_resp.status_code = status_code
        mock_resp.text = text_data
        mock_resp.content = content_data
        mock_resp.reason = "OK" if status_code < 400 else "Error"
        if json_data:
            mock_resp.json = MagicMock(return_value=json_data)
        
        if raise_for_status_exception:
            mock_resp.raise_for_status = MagicMock(side_effect=raise_for_status_exception)
        else:
            mock_resp.raise_for_status = MagicMock()
            if status_code >= 400:
                error_response_mock = MagicMock(spec=requests.Response)
                error_response_mock.status_code = status_code
                error_response_mock.reason = mock_resp.reason
                mock_resp.raise_for_status.side_effect = requests.exceptions.HTTPError(
                    f"{status_code} Client Error: {mock_resp.reason} for url: FAKE_URL", response=error_response_mock
                )
        return mock_resp

    # --- Test Argument Parsing ---
    def test_argument_parsing_defaults(self, mock_get_user_agent):
        # Test that default values are set correctly
        # We need to capture the output of parse_args
        with patch.object(sys, 'argv', ['blast_ncbi_tabular.py', '--query', 'ACGT']):
            args = blast_ncbi_tabular.main.__globals__['parser'].parse_args() # Access parser from main
            self.assertEqual(args.query, "ACGT")
            self.assertEqual(args.database, "nr")
            self.assertEqual(args.program, "blastp")
            self.assertEqual(args.filter, "L")
            self.assertEqual(args.hitlist_size, 20)
            self.assertEqual(args.evalue, 0.01)
            self.assertEqual(args.poll_interval, 10)
            self.assertEqual(args.max_poll_attempts, 120)
            self.assertEqual(args.ncbi_url, blast_ncbi_tabular.DEFAULT_NCBI_URL)
    
    def test_argument_parsing_required(self, mock_get_user_agent):
        # Test that required arguments are enforced
        with patch.object(sys, 'argv', ['blast_ncbi_tabular.py']):
            with self.assertRaises(SystemExit): # argparse exits on error
                with redirect_stderr(io.StringIO()) as _: # Suppress stderr
                    blast_ncbi_tabular.main.__globals__['parser'].parse_args()

    def test_argument_parsing_custom_values(self, mock_get_user_agent):
        custom_args = [
            'blast_ncbi_tabular.py',
            '--query', 'NNNN',
            '--database', 'refseq_protein',
            '--program', 'blastx',
            '--filter', 'F',
            '--hitlist_size', '50',
            '--evalue', '0.05',
            '--poll_interval', '5',
            '--max_poll_attempts', '60',
            '--ncbi_url', 'http://localhost/blast'
        ]
        with patch.object(sys, 'argv', custom_args):
            args = blast_ncbi_tabular.main.__globals__['parser'].parse_args()
            self.assertEqual(args.query, 'NNNN')
            self.assertEqual(args.database, 'refseq_protein')
            self.assertEqual(args.program, 'blastx')
            self.assertEqual(args.filter, 'F')
            self.assertEqual(args.hitlist_size, 50)
            self.assertEqual(args.evalue, 0.05)
            self.assertEqual(args.poll_interval, 5)
            self.assertEqual(args.max_poll_attempts, 60)
            self.assertEqual(args.ncbi_url, 'http://localhost/blast')

    # --- Test submit_blast_query ---
    @patch('requests.post')
    def test_submit_blast_query_success_rid_html(self, mock_post, mock_get_user_agent):
        mock_response_text = '<html><body><input type="hidden" name="RID" value="TEST_RID_123"/> ... QBlastInfoBegin RID = OTHER_RID QBlastInfoEnd</body></html>'
        mock_post.return_value = self._mock_response(text_data=mock_response_text)
        
        rid = blast_ncbi_tabular.submit_blast_query("ACGT", "nr", "blastp", "L", 20, 0.01, "fake_url")
        self.assertEqual(rid, "TEST_RID_123")
        mock_post.assert_called_once()
        self.assertEqual(mock_post.call_args[1]['headers']['User-Agent'], MOCK_USER_AGENT)


    @patch('requests.post')
    def test_submit_blast_query_success_rid_text(self, mock_post, mock_get_user_agent):
        mock_response_text = "\n    RID = TEST_RID_456\n"
        mock_post.return_value = self._mock_response(text_data=mock_response_text)
        
        rid = blast_ncbi_tabular.submit_blast_query("ACGT", "nr", "blastp", "L", 20, 0.01, "fake_url")
        self.assertEqual(rid, "TEST_RID_456")

    @patch('requests.post')
    def test_submit_blast_query_success_rid_qblastinfo(self, mock_post, mock_get_user_agent):
        mock_response_text = "Some other text\nQBlastInfoBegin\n     RID = TEST_RID_789\nQBlastInfoEnd\nMore text"
        mock_post.return_value = self._mock_response(text_data=mock_response_text)
        
        rid = blast_ncbi_tabular.submit_blast_query("ACGT", "nr", "blastp", "L", 20, 0.01, "fake_url")
        self.assertEqual(rid, "TEST_RID_789")

    @patch('requests.post')
    def test_submit_blast_query_no_rid_found(self, mock_post, mock_get_user_agent):
        mock_post.return_value = self._mock_response(text_data="No RID here.")
        with redirect_stderr(io.StringIO()) as stderr_capture:
            rid = blast_ncbi_tabular.submit_blast_query("ACGT", "nr", "blastp", "L", 20, 0.01, "fake_url")
        self.assertIsNone(rid)
        self.assertIn("Error: Could not find RID in NCBI response.", stderr_capture.getvalue())

    @patch('requests.post')
    def test_submit_blast_query_request_exception(self, mock_post, mock_get_user_agent):
        mock_post.side_effect = requests.exceptions.RequestException("Test network error")
        with redirect_stderr(io.StringIO()) as stderr_capture:
            rid = blast_ncbi_tabular.submit_blast_query("ACGT", "nr", "blastp", "L", 20, 0.01, "fake_url")
        self.assertIsNone(rid)
        self.assertIn("Error submitting BLAST request: Test network error", stderr_capture.getvalue())

    # --- Test poll_blast_results ---
    @patch('requests.post') # Changed from get to post based on script
    @patch('time.sleep', return_value=None) # Mock time.sleep
    def test_poll_blast_results_success_ready(self, mock_sleep, mock_post, mock_get_user_agent):
        sample_tabular_data = "# BLASTP 2.10.0+\n# Query: test_query\n# Database: nr\n# Fields: query id, subject id, ...\n# 0 hits found"
        mock_post.side_effect = [
            self._mock_response(text_data="Status=WAITING"),
            self._mock_response(text_data="Status=READY\n" + sample_tabular_data)
        ]
        results = blast_ncbi_tabular.poll_blast_results("TEST_RID", 1, 3, "fake_url")
        self.assertEqual(results, "Status=READY\n" + sample_tabular_data)
        self.assertEqual(mock_post.call_count, 2)
        self.assertEqual(mock_sleep.call_count, 2) # Sleep called before each poll attempt

    @patch('requests.post')
    @patch('time.sleep', return_value=None)
    def test_poll_blast_results_no_hits(self, mock_sleep, mock_post, mock_get_user_agent):
        ready_no_hits = "Status=READY\n\n# BLASTP 2.13.0+\n# Query: myseq\n# Database: nr\n# 0 hits found\n"
        mock_post.return_value = self._mock_response(text_data=ready_no_hits)
        with redirect_stdout(io.StringIO()) as stdout_capture:
             results = blast_ncbi_tabular.poll_blast_results("TEST_RID", 1, 3, "fake_url")
        self.assertEqual(results, "") # Returns empty string for "No hits found"
        self.assertIn("No hits found for the query.", stdout_capture.getvalue())


    @patch('requests.post')
    @patch('time.sleep', return_value=None)
    def test_poll_blast_results_failed(self, mock_sleep, mock_post, mock_get_user_agent):
        mock_post.return_value = self._mock_response(text_data="Status=FAILED\nSome error details.")
        with redirect_stderr(io.StringIO()) as stderr_capture:
            results = blast_ncbi_tabular.poll_blast_results("TEST_RID", 1, 3, "fake_url")
        self.assertIsNone(results)
        self.assertIn("Error: BLAST query failed.", stderr_capture.getvalue())

    @patch('requests.post')
    @patch('time.sleep', return_value=None)
    def test_poll_blast_results_timeout(self, mock_sleep, mock_post, mock_get_user_agent):
        mock_post.return_value = self._mock_response(text_data="Status=WAITING")
        with redirect_stderr(io.StringIO()) as stderr_capture:
            results = blast_ncbi_tabular.poll_blast_results("TEST_RID", 1, 3, "fake_url") # 3 attempts
        self.assertIsNone(results)
        self.assertIn("Error: Timed out waiting for BLAST results.", stderr_capture.getvalue())
        self.assertEqual(mock_post.call_count, 3)
        self.assertEqual(mock_sleep.call_count, 3)
        
    @patch('requests.post')
    @patch('time.sleep', return_value=None)
    def test_poll_blast_results_status_unknown(self, mock_sleep, mock_post, mock_get_user_agent):
        mock_post.return_value = self._mock_response(text_data="Status=UNKNOWN")
        with redirect_stderr(io.StringIO()) as stderr_capture:
            results = blast_ncbi_tabular.poll_blast_results("TEST_RID", 1, 3, "fake_url")
        self.assertIsNone(results)
        self.assertIn("Error: BLAST query RID expired or unknown.", stderr_capture.getvalue())

    @patch('requests.post')
    @patch('time.sleep', return_value=None)
    def test_poll_blast_results_assumed_ready(self, mock_sleep, mock_post, mock_get_user_agent):
        tabular_data_no_status = "# BLASTP 2.10.0+\n# Query: test_query\n# Database: nr\n# Fields: query id, subject id, ...\nseq1\tsubj1\t..."
        mock_post.return_value = self._mock_response(text_data=tabular_data_no_status)
        with redirect_stderr(io.StringIO()) as stderr_capture:
            results = blast_ncbi_tabular.poll_blast_results("TEST_RID", 1, 3, "fake_url")
        self.assertEqual(results, tabular_data_no_status)
        self.assertIn("Status: Assumed READY (tabular data found without explicit Status=READY).", stderr_capture.getvalue())

    # --- Test parse_and_print_tabular_results ---
    def test_parse_and_print_tabular_results_success(self, mock_get_user_agent):
        sample_report = """
# BLASTP 2.10.0+
# Query: myquery
# Database: nr
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 2 hits found
myquery	gb|AAA123.1|	100.00	150	0	0	1	150	10	160	1.0e-50	300
myquery	ref|XP_001.1|sp|P12345.1|	95.50	100	5	0	20	120	30	130	2.5e-30	200
"""
        expected_output_lines = [
            "#   Subject ID (Processed)         Identity  Length   E-value",
            "1   gb:AAA123.1                      100.00     150   1.0e-50",
            "2   ref:XP_001.1 sp:P12345.1          95.50     100   2.5e-30"
        ]
        with redirect_stdout(io.StringIO()) as stdout_capture:
            blast_ncbi_tabular.parse_and_print_tabular_results(sample_report)
        
        output = stdout_capture.getvalue().strip()
        # Normalize whitespace for comparison if needed, but f-string formatting should be consistent
        for expected_line in expected_output_lines:
            self.assertIn(expected_line, output)

    def test_parse_and_print_tabular_results_no_hits_report(self, mock_get_user_agent):
        sample_report_no_hits = """
# BLASTP 2.13.0+
# Query: Q00000
# Database: swissprot
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 0 hits found
"""
        with redirect_stdout(io.StringIO()) as stdout_capture, \
             redirect_stderr(io.StringIO()) as stderr_capture:
            blast_ncbi_tabular.parse_and_print_tabular_results(sample_report_no_hits)
        
        self.assertEqual(stdout_capture.getvalue().strip(), "") # No table printed
        # The script currently prints "No results to parse." to stderr if report is empty,
        # but if it's a "0 hits found" report, it might print nothing or the header.
        # Based on current parse_and_print, it should print the header then nothing else.
        # Let's adjust expectation: parse_and_print will print its header.
        # self.assertIn("#   Subject ID (Processed)         Identity  Length   E-value", stdout_capture.getvalue())
        # Actually, if "# 0 hits found" is present, it prints the header and then loops 0 times.
        # If the input `report_text` itself is empty (e.g. from poll_blast_results returning "" for "No hits found")
        # then "No results to parse." is printed to stderr.

    def test_parse_and_print_tabular_results_empty_input(self, mock_get_user_agent):
        with redirect_stderr(io.StringIO()) as stderr_capture:
            blast_ncbi_tabular.parse_and_print_tabular_results("")
        self.assertIn("No results to parse.", stderr_capture.getvalue())
        
    def test_parse_and_print_tabular_results_missing_fields_line(self, mock_get_user_agent):
        sample_report_no_fields = """
# BLASTP 2.10.0+
# Query: myquery
myquery	gb|AAA123.1|	100.00	150	0	0	1	150	10	160	1.0e-50	300
"""     # Fallback default headers are used
        with redirect_stderr(io.StringIO()) as stderr_capture, \
             redirect_stdout(io.StringIO()) as stdout_capture:
            blast_ncbi_tabular.parse_and_print_tabular_results(sample_report_no_fields)
        
        self.assertIn("Warning: '# Fields:' line not found. Attempting to parse data assuming standard columns.", stderr_capture.getvalue())
        self.assertIn("1   gb:AAA123.1                      100.00     150   1.0e-50", stdout_capture.getvalue())


    # --- Test Main Script Logic ---
    @patch('blast_ncbi_tabular.submit_blast_query')
    @patch('blast_ncbi_tabular.poll_blast_results')
    @patch('blast_ncbi_tabular.parse_and_print_tabular_results')
    @patch('time.sleep', return_value=None) # Mock initial sleep if any in main
    def test_main_successful_flow(self, mock_tsleep, mock_parse_print, mock_poll, mock_submit, mock_get_user_agent):
        mock_submit.return_value = "DUMMY_RID"
        mock_poll.return_value = "Sample BLAST report text"
        
        test_args = ['blast_ncbi_tabular.py', '--query', 'ACGT']
        with patch.object(sys, 'argv', test_args):
            with redirect_stderr(io.StringIO()) as stderr_capture:
                 with self.assertRaises(SystemExit) as cm: # Script calls sys.exit()
                    blast_ncbi_tabular.main()
                 self.assertEqual(cm.exception.code, None) # Should be sys.exit(0) or just return
        
        mock_submit.assert_called_once()
        mock_poll.assert_called_once_with("DUMMY_RID", 10, 120, blast_ncbi_tabular.DEFAULT_NCBI_URL)
        mock_parse_print.assert_called_once_with("Sample BLAST report text")
        self.assertIn("BLAST processing finished.", stderr_capture.getvalue())


    @patch('blast_ncbi_tabular.submit_blast_query')
    @patch('time.sleep', return_value=None)
    def test_main_submit_fails(self, mock_tsleep, mock_submit, mock_get_user_agent):
        mock_submit.return_value = None # Simulate submission failure
        
        test_args = ['blast_ncbi_tabular.py', '--query', 'ACGT']
        with patch.object(sys, 'argv', test_args):
            with redirect_stderr(io.StringIO()) as stderr_capture:
                with self.assertRaises(SystemExit) as cm:
                    blast_ncbi_tabular.main()
                self.assertEqual(cm.exception.code, 1) # Should exit with error
        
        self.assertIn("Exiting due to error in BLAST submission.", stderr_capture.getvalue())

    @patch('blast_ncbi_tabular.submit_blast_query')
    @patch('blast_ncbi_tabular.poll_blast_results')
    @patch('time.sleep', return_value=None)
    def test_main_poll_fails(self, mock_tsleep, mock_poll, mock_submit, mock_get_user_agent):
        mock_submit.return_value = "DUMMY_RID"
        mock_poll.return_value = None # Simulate polling failure
        
        test_args = ['blast_ncbi_tabular.py', '--query', 'ACGT']
        with patch.object(sys, 'argv', test_args):
            with redirect_stderr(io.StringIO()) as stderr_capture:
                with self.assertRaises(SystemExit) as cm:
                    blast_ncbi_tabular.main()
                self.assertEqual(cm.exception.code, 1)
        
        self.assertIn("Exiting due to error or timeout during polling.", stderr_capture.getvalue())

    def test_parse_subject_id(self, mock_get_user_agent):
        self.assertEqual(blast_ncbi_tabular.parse_subject_id("gb|AEG74000.1|emb|CAM39480.1|"), "gb:AEG74000.1 emb:CAM39480.1")
        self.assertEqual(blast_ncbi_tabular.parse_subject_id("ref|XP_001.2|"), "ref:XP_001.2")
        self.assertEqual(blast_ncbi_tabular.parse_subject_id("pdb|1ABC|A"), "pdb:1ABC A") # Handles odd number of parts
        self.assertEqual(blast_ncbi_tabular.parse_subject_id("XYZ123"), "XYZ123") # Single ID
        self.assertEqual(blast_ncbi_tabular.parse_subject_id(""), "")
        self.assertEqual(blast_ncbi_tabular.parse_subject_id("gi|12345|gb|AEG74000.1|"), "gi:12345 gb:AEG74000.1")


if __name__ == '__main__':
    # Need to ensure that the utility_functions mock is applied correctly when tests run.
    # One way is to ensure the class-level patch is correctly specified.
    # The path for patching should be 'blast_ncbi_tabular.utility_functions' if the script imports it as 'import utility_functions'
    # and utility_functions.py is effectively in its sys.path.
    # Given the script's structure: `from utils import utility_functions`
    # and the test setup adding `pylib` to sys.path, the module `utils` (i.e. `pylib.utils`) should be findable.
    # So, patching `utils.utility_functions.get_user_agent` at the class level is one way.
    # Or, if `blast_ncbi_tabular` stores `utility_functions` as a module attribute:
    # @patch('blast_ncbi_tabular.utility_functions.get_user_agent', return_value=MOCK_USER_AGENT)
    
    # The current class decorator @patch('utils.utility_functions.get_user_agent', ...) is likely correct.
    unittest.main(argv=['first-arg-is-ignored'], exit=False, verbosity=2)
