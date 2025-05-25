import unittest
from unittest.mock import patch, MagicMock, mock_open
import os
import sys
import io
from contextlib import redirect_stdout

# Add the directory of this test script to the Python path
# to ensure that the script being tested, 'blast_ncbi_seq_def', can be imported.
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_dir)

import blast_ncbi_seq_def

# Import requests for its exception types
import requests


class TestBlastNcbiSeqDef(unittest.TestCase):

    def setUp(self):
        """Set up for test methods."""
        self.test_fasta_content_single = ">seq1\nACGT"
        self.test_fasta_content_multiple = ">seq1\nACGT\n>seq2\nGCTA"
        self.test_fasta_content_comments_empty = "#comment\n\n>seq1\nACGT\n #another comment \n  \n>seq2\nGCTA  "  # Added trailing space to GCTA line
        self.test_fasta_content_empty_file = ""
        self.test_fasta_content_malformed = "ACGT\nTGCA"  # No header
        self.test_fasta_whitespace = ">seqWS\nACG T\n GCTA "  # Added trailing space

        self.dummy_fasta_path = os.path.join(current_dir, "dummy_test.fasta")

    def tearDown(self):
        """Tear down after test methods."""
        if os.path.exists(self.dummy_fasta_path):
            os.remove(self.dummy_fasta_path)

    def _create_dummy_fasta(self, content):
        with open(self.dummy_fasta_path, "w") as f:
            f.write(content)
        return self.dummy_fasta_path

    # --- Test FASTA Parsing ---
    def test_parse_fasta_single_sequence(self):
        path = self._create_dummy_fasta(self.test_fasta_content_single)
        sequences = blast_ncbi_seq_def.parse_fasta(path)
        self.assertEqual(len(sequences), 1)
        self.assertEqual(sequences["seq1"], "ACGT")

    def test_parse_fasta_multiple_sequences(self):
        path = self._create_dummy_fasta(self.test_fasta_content_multiple)
        sequences = blast_ncbi_seq_def.parse_fasta(path)
        self.assertEqual(len(sequences), 2)
        self.assertEqual(sequences["seq1"], "ACGT")
        self.assertEqual(sequences["seq2"], "GCTA")

    def test_parse_fasta_with_comments_and_empty_lines_and_trailing_spaces(self):
        path = self._create_dummy_fasta(self.test_fasta_content_comments_empty)
        sequences = blast_ncbi_seq_def.parse_fasta(path)
        self.assertEqual(len(sequences), 2)
        self.assertEqual(sequences["seq1"], "ACGT")
        self.assertEqual(sequences["seq2"], "GCTA")  # Whitespace removal check

    def test_parse_fasta_empty_file(self):
        path = self._create_dummy_fasta(self.test_fasta_content_empty_file)
        sequences = blast_ncbi_seq_def.parse_fasta(path)
        self.assertEqual(len(sequences), 0)

    def test_parse_fasta_malformed_no_header(self):
        path = self._create_dummy_fasta(self.test_fasta_content_malformed)
        sequences = blast_ncbi_seq_def.parse_fasta(path)
        self.assertEqual(len(sequences), 0)

    def test_parse_fasta_whitespace_in_sequence(self):
        path = self._create_dummy_fasta(self.test_fasta_whitespace)
        sequences = blast_ncbi_seq_def.parse_fasta(path)
        self.assertEqual(sequences["seqWS"], "ACGTGCTA")

    # --- Mock responses for requests.post ---
    def _mock_requests_post_response(
        self,
        status_code=200,
        text_data="",
        json_data=None,
        raise_for_status_exception=None,
    ):
        mock_resp = MagicMock(spec=requests.Response)
        mock_resp.status_code = status_code
        mock_resp.text = text_data
        mock_resp.reason = "OK" if status_code < 400 else "Error"
        if json_data:
            mock_resp.json = MagicMock(return_value=json_data)

        if raise_for_status_exception:
            mock_resp.raise_for_status = MagicMock(
                side_effect=raise_for_status_exception
            )
        else:
            mock_resp.raise_for_status = MagicMock()
            if status_code >= 400:
                # Create a mock response object to be the 'response' attribute of HTTPError
                error_response_mock = MagicMock(spec=requests.Response)
                error_response_mock.status_code = status_code
                error_response_mock.reason = mock_resp.reason
                mock_resp.raise_for_status.side_effect = requests.exceptions.HTTPError(
                    f"{status_code} Client Error: {mock_resp.reason} for url: FAKE_URL",
                    response=error_response_mock,
                )
        return mock_resp

    # --- Test BLAST Workflow ---
    # Note: The MockHSP class and its related comment `Biopython uses 1 for plus, -1 for minus`
    # were previously here but were removed as they are not used in this file.
    # This test suite focuses on the script's interaction with NCBI BLAST services.
    @patch("blast_ncbi_seq_def.requests.post")
    @patch("blast_ncbi_seq_def.time.sleep", return_value=None)
    def test_successful_blast_workflow(self, mock_sleep, mock_post):
        sample_report_text = "Status=READY\n\nSample BLAST Report Content"
        mock_post.side_effect = [
            self._mock_requests_post_response(text_data="RID = TESTRID123"),
            self._mock_requests_post_response(text_data="Status=WAITING"),
            self._mock_requests_post_response(text_data=sample_report_text),
        ]

        path = self._create_dummy_fasta(self.test_fasta_content_single)
        captured_output = io.StringIO()
        with redirect_stdout(captured_output):
            test_args = ["blast_ncbi_seq_def.py", path, "0.01"]
            with patch.object(sys, "argv", test_args):
                with self.assertRaises(SystemExit) as cm:
                    blast_ncbi_seq_def.main()
                self.assertIsNone(
                    cm.exception.code, "Script should exit cleanly on success."
                )

        output = captured_output.getvalue()
        self.assertIn(f"Input FASTA file: {path}", output)
        self.assertIn("E-value: 0.01", output)
        self.assertIn("Processing sequence: seq1", output)
        self.assertIn("1\tWAITING", output)
        self.assertIn(sample_report_text, output)
        self.assertNotIn("bad input", output)
        self.assertNotIn("timed out", output)
        self.assertEqual(mock_post.call_count, 3)
        self.assertEqual(mock_sleep.call_count, 3)  # 1 initial, 2 polling

    @patch("blast_ncbi_seq_def.requests.post")
    @patch("blast_ncbi_seq_def.time.sleep", return_value=None)
    def test_blast_rid_unknown(self, mock_sleep, mock_post):
        mock_post.return_value = self._mock_requests_post_response(
            text_data="RID = UNKNOWN"
        )
        path = self._create_dummy_fasta(self.test_fasta_content_single)

        captured_output = io.StringIO()
        with redirect_stdout(captured_output):
            test_args = ["blast_ncbi_seq_def.py", path, "0.01"]
            with patch.object(sys, "argv", test_args):
                with self.assertRaises(SystemExit) as cm:
                    blast_ncbi_seq_def.main()
                self.assertIsNone(
                    cm.exception.code, "Script should exit after 'bad input'."
                )

        output = captured_output.getvalue()
        self.assertIn("Processing sequence: seq1", output)
        self.assertIn("bad input", output)
        self.assertNotIn("WAITING", output)
        self.assertEqual(mock_post.call_count, 1)
        self.assertEqual(mock_sleep.call_count, 1)  # Only initial sleep

    @patch("blast_ncbi_seq_def.requests.post")
    @patch("blast_ncbi_seq_def.time.sleep", return_value=None)
    def test_blast_timeout(self, mock_sleep, mock_post):
        mock_post.side_effect = [
            self._mock_requests_post_response(text_data="RID = TESTRID456")
        ] + [
            self._mock_requests_post_response(text_data="Status=WAITING")
        ] * blast_ncbi_seq_def.MAX_ATTEMPTS

        path = self._create_dummy_fasta(self.test_fasta_content_single)
        captured_output = io.StringIO()
        with redirect_stdout(captured_output):
            test_args = ["blast_ncbi_seq_def.py", path, "0.01"]
            with patch.object(sys, "argv", test_args):
                with self.assertRaises(SystemExit) as cm:
                    blast_ncbi_seq_def.main()
                self.assertIsNone(
                    cm.exception.code, "Script should exit after timeout."
                )

        output = captured_output.getvalue()
        self.assertIn("Processing sequence: seq1", output)
        self.assertIn("1\tWAITING", output)
        self.assertIn(f"{blast_ncbi_seq_def.MAX_ATTEMPTS}\tWAITING", output)
        self.assertIn("timed out at 20 minutes", output)
        self.assertNotIn("READY", output)
        self.assertNotIn("bad input", output)
        self.assertEqual(mock_post.call_count, 1 + blast_ncbi_seq_def.MAX_ATTEMPTS)
        self.assertEqual(mock_sleep.call_count, 1 + blast_ncbi_seq_def.MAX_ATTEMPTS)

    @patch("blast_ncbi_seq_def.requests.post")
    @patch("blast_ncbi_seq_def.time.sleep", return_value=None)
    def test_blast_submit_http_error(self, mock_sleep, mock_post):
        http_error_exception = requests.exceptions.RequestException(
            "Test HTTP Error on Submit"
        )
        mock_post.return_value = self._mock_requests_post_response(
            status_code=500, raise_for_status_exception=http_error_exception
        )

        path = self._create_dummy_fasta(self.test_fasta_content_single)
        captured_output = io.StringIO()
        with redirect_stdout(captured_output):
            test_args = ["blast_ncbi_seq_def.py", path, "0.01"]
            with patch.object(sys, "argv", test_args):
                with self.assertRaises(SystemExit) as cm:
                    blast_ncbi_seq_def.main()
                self.assertIsNone(
                    cm.exception.code, "Script should exit on submit error."
                )

        output = captured_output.getvalue()
        self.assertIn(
            "Error submitting BLAST request: Test HTTP Error on Submit", output
        )
        self.assertEqual(mock_post.call_count, 1)
        self.assertEqual(mock_sleep.call_count, 1)

    @patch("blast_ncbi_seq_def.requests.post")
    @patch("blast_ncbi_seq_def.time.sleep", return_value=None)
    def test_blast_poll_http_error(self, mock_sleep, mock_post):
        http_error_exception = requests.exceptions.RequestException(
            "Test HTTP Error on Poll"
        )
        mock_post.side_effect = [
            self._mock_requests_post_response(text_data="RID = TESTRID789"),
            self._mock_requests_post_response(
                status_code=500, raise_for_status_exception=http_error_exception
            ),
        ]

        path = self._create_dummy_fasta(self.test_fasta_content_single)
        captured_output = io.StringIO()
        with redirect_stdout(captured_output):
            test_args = ["blast_ncbi_seq_def.py", path, "0.01"]
            with patch.object(sys, "argv", test_args):
                with self.assertRaises(SystemExit) as cm:
                    blast_ncbi_seq_def.main()
                self.assertIsNone(
                    cm.exception.code, "Script should exit on poll error."
                )

        output = captured_output.getvalue()
        self.assertIn("Error polling for results: Test HTTP Error on Poll", output)
        self.assertEqual(mock_post.call_count, 2)
        self.assertEqual(mock_sleep.call_count, 1 + 1)

    @patch("blast_ncbi_seq_def.time.sleep", return_value=None)
    @patch("blast_ncbi_seq_def.submit_blast_request")
    @patch("blast_ncbi_seq_def.poll_for_results")
    def test_argument_usage_in_main(self, mock_poll, mock_submit, mock_initial_sleep):
        mock_submit.return_value = "DUMMY_RID_ARG_TEST"
        path = self._create_dummy_fasta(self.test_fasta_content_single)
        custom_e_value = "0.005"
        test_args = ["blast_ncbi_seq_def.py", path, custom_e_value]

        with patch.object(sys, "argv", test_args):
            with self.assertRaises(SystemExit) as cm:
                blast_ncbi_seq_def.main()
            self.assertIsNone(
                cm.exception.code, "main() should exit cleanly with mocked functions."
            )

        mock_initial_sleep.assert_any_call(20)
        mock_submit.assert_called_once_with("ACGT", custom_e_value)
        mock_poll.assert_called_once_with("DUMMY_RID_ARG_TEST")

    @patch("blast_ncbi_seq_def.requests.post")
    @patch("blast_ncbi_seq_def.time.sleep", return_value=None)
    def test_blast_rid_not_found_general(self, mock_sleep, mock_post):
        mock_post.return_value = self._mock_requests_post_response(
            text_data="Some other response without RID"
        )
        path = self._create_dummy_fasta(self.test_fasta_content_single)

        captured_output = io.StringIO()
        with redirect_stdout(captured_output):
            test_args = ["blast_ncbi_seq_def.py", path, "0.01"]
            with patch.object(sys, "argv", test_args):
                with self.assertRaises(SystemExit) as cm:
                    blast_ncbi_seq_def.main()
                self.assertIsNone(
                    cm.exception.code, "Script should exit on RID not found."
                )

        output = captured_output.getvalue()
        self.assertIn("Could not find RID in response:", output)
        self.assertIn("Some other response without RID", output)
        self.assertEqual(mock_post.call_count, 1)
        self.assertEqual(mock_sleep.call_count, 1)

    @patch("blast_ncbi_seq_def.requests.post")
    @patch("blast_ncbi_seq_def.time.sleep", return_value=None)
    def test_blast_status_not_found(self, mock_sleep, mock_post):
        mock_post.side_effect = [
            self._mock_requests_post_response(text_data="RID = TESTRID_NO_STATUS"),
            self._mock_requests_post_response(
                text_data="Some response without Status="
            ),
        ]
        path = self._create_dummy_fasta(self.test_fasta_content_single)

        captured_output = io.StringIO()
        with redirect_stdout(captured_output):
            test_args = ["blast_ncbi_seq_def.py", path, "0.01"]
            with patch.object(sys, "argv", test_args):
                with self.assertRaises(SystemExit) as cm:
                    blast_ncbi_seq_def.main()
                self.assertIsNone(
                    cm.exception.code, "Script should exit on Status not found."
                )

        output = captured_output.getvalue()
        self.assertIn("Could not find Status in response:", output)
        self.assertIn("Some response without Status=", output)
        self.assertEqual(mock_post.call_count, 2)
        self.assertEqual(mock_sleep.call_count, 2)

    @patch("blast_ncbi_seq_def.requests.post")
    @patch("blast_ncbi_seq_def.time.sleep", return_value=None)
    def test_blast_status_not_found_but_ready_lenient(self, mock_sleep, mock_post):
        # This text should contain "READY" but NOT "Status=" to trigger the lenient path
        ready_response_text = "The job is READY. Results are available below."
        mock_post.side_effect = [
            self._mock_requests_post_response(text_data="RID = TESTRID_LENIENT"),
            self._mock_requests_post_response(text_data=ready_response_text),
        ]
        path = self._create_dummy_fasta(self.test_fasta_content_single)

        captured_output = io.StringIO()
        with redirect_stdout(captured_output):
            test_args = ["blast_ncbi_seq_def.py", path, "0.01"]
            with patch.object(sys, "argv", test_args):
                with self.assertRaises(SystemExit) as cm:
                    blast_ncbi_seq_def.main()
                self.assertIsNone(
                    cm.exception.code,
                    "Script should exit (even if Status not found but READY in text).",
                )

        output = captured_output.getvalue()
        self.assertIn("Could not find Status in response:", output)
        self.assertIn(
            "Assuming READY based on content, though Status line missing.", output
        )
        self.assertIn(ready_response_text, output)
        self.assertEqual(mock_post.call_count, 2)
        self.assertEqual(mock_sleep.call_count, 2)

    # Test initial 20s sleep
    @patch("blast_ncbi_seq_def.time.sleep")
    @patch("blast_ncbi_seq_def.submit_blast_request", return_value="RID_FOR_SLEEP_TEST")
    @patch("blast_ncbi_seq_def.poll_for_results") # Mocked to prevent further execution
    def test_initial_delay(self, mock_poll, mock_submit, mock_sleep):
        path = self._create_dummy_fasta(self.test_fasta_content_single)
        test_args = ["blast_ncbi_seq_def.py", path, "0.01"]

        with patch.object(sys, "argv", test_args):
            with self.assertRaises(SystemExit) as cm: # Capture context manager
                blast_ncbi_seq_def.main()
            self.assertIsNone( # Check for clean exit if main completes (even if mocked)
                cm.exception.code, 
                "main() should exit cleanly with mocked poll/submit."
            )

        # Check if time.sleep(20) was the first call to sleep
        mock_sleep.assert_any_call(20)


if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0]] + sys.argv[1:], verbosity=2, exit=False) # Ensure single newline at EOF
