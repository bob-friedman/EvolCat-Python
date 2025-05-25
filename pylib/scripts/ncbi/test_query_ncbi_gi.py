import unittest
from unittest.mock import patch, MagicMock, call
import os
import sys
import io
from contextlib import redirect_stdout

# Add the directory containing the script to the Python path
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_dir)

import query_ncbi_gi

# Import requests for its exception types
import requests


class TestQueryNcbiGi(unittest.TestCase):

    def setUp(self):
        """Set up for test methods."""
        self.search_term_ok = "BRCA1"
        self.search_term_no_results = "TERM_WITH_NO_RESULTS_XYZ"

        # Mocked NCBI search response HTML
        self.mock_search_html_multiple_uids = """
        <html><body>
        <input name="EntrezSystem2.PEntrez.Nuccore.Sequence_ResultsPanel.Sequence_RVDocSum.uid" value="12345">
        <input name="EntrezSystem2.PEntrez.Nuccore.Sequence_ResultsPanel.Sequence_RVDocSum.uid" value="67890">
        </body></html>
        """
        self.mock_search_html_single_uid = """
        <html><body>
        <input name="EntrezSystem2.PEntrez.Nuccore.Sequence_ResultsPanel.Sequence_RVDocSum.uid" value="12345">
        </body></html>
        """
        self.mock_search_html_no_uids = (
            "<html><body><p>No results found.</p></body></html>"
        )

        # Mocked GenBank record
        self.mock_genbank_record_12345 = """LOCUS       12345               100 bp    DNA     linear   SYN 15-JUL-2024
DEFINITION  Sample sequence 12345.
ACCESSION   12345
VERSION     12345.1
KEYWORDS    .
SOURCE      synthetic construct
  ORGANISM  synthetic construct
            other sequences; artificial sequences.
FEATURES             Location/Qualifiers
     source          1..100
                     /organism="synthetic construct"
                     /mol_type="genomic DNA"
ORIGIN
        1 gatcctccat atacaacggt atctccacct caggtttaga tctcaacaac ggaaccattg
       61 ccgacatgag acagttaggt atcgtcgaga gttacaagct
//
"""
        self.mock_genbank_record_67890 = """LOCUS       67890               100 bp    DNA     linear   SYN 15-JUL-2024
DEFINITION  Sample sequence 67890.
ACCESSION   67890
VERSION     67890.1
KEYWORDS    .
SOURCE      synthetic construct
  ORGANISM  synthetic construct
            other sequences; artificial sequences.
FEATURES             Location/Qualifiers
     source          1..100
                     /organism="synthetic construct"
                     /mol_type="genomic DNA"
ORIGIN
        1 agatcctccg atacaacggt atctccacct caggtttaga tctcaacaac ggaaccattg
       61 ccgacatgag acagttaggt atcgtcgaga gttacaagct
//
"""
        self.mock_invalid_genbank_record = "This is not a GenBank record."
        self.mock_empty_genbank_record = ""

    def tearDown(self):
        """Tear down after test methods."""
        pass  # No files created that need cleanup

    def _mock_response(
        self, status_code=200, text_data="", raise_for_status_exception=None
    ):
        mock_resp = MagicMock(spec=requests.Response)
        mock_resp.status_code = status_code
        mock_resp.text = text_data
        mock_resp.reason = "OK" if status_code < 400 else "Error"
        if raise_for_status_exception:
            mock_resp.raise_for_status = MagicMock(
                side_effect=raise_for_status_exception
            )
        else:
            mock_resp.raise_for_status = MagicMock()
            if status_code >= 400:
                mock_resp.raise_for_status.side_effect = requests.exceptions.HTTPError(
                    f"{status_code} Client Error: {mock_resp.reason} for url: FAKE_URL",
                    response=mock_resp,
                )
        return mock_resp

    # --- Test Argument Parsing (indirectly via main function) ---
    @patch(
        "query_ncbi_gi.search_ncbi", return_value=[]
    )  # Mock search to prevent further calls
    @patch("query_ncbi_gi.time.sleep")  # Mock all sleeps
    def test_argument_parsing(self, mock_sleep, mock_search_ncbi):
        test_term = "MyTestTerm"
        test_args = ["query_ncbi_gi.py", test_term]

        captured_output = io.StringIO()
        with redirect_stdout(captured_output):
            with patch.object(sys, "argv", test_args):
                query_ncbi_gi.main()

        output = captured_output.getvalue()
        self.assertIn(f"Search Term: {test_term}", output)
        mock_search_ncbi.assert_called_once_with(test_term)
        # Check final sleep
        self.assertIn(call(5), mock_sleep.call_args_list)

    # --- Test Search Workflow and UID Parsing ---
    @patch("query_ncbi_gi.requests.post")
    @patch("query_ncbi_gi.requests.get")  # Mock get as well for fetch
    @patch("query_ncbi_gi.time.sleep")
    def test_search_successful_multiple_uids(self, mock_sleep, mock_get, mock_post):
        # Search returns multiple UIDs, fetch returns valid GenBank for each
        mock_post.return_value = self._mock_response(
            text_data=self.mock_search_html_multiple_uids
        )
        mock_get.side_effect = [
            self._mock_response(text_data=self.mock_genbank_record_12345),
            self._mock_response(text_data=self.mock_genbank_record_67890),
        ]

        test_args = ["query_ncbi_gi.py", self.search_term_ok]
        captured_output = io.StringIO()
        with redirect_stdout(captured_output):
            with patch.object(sys, "argv", test_args):
                query_ncbi_gi.main()

        output = captured_output.getvalue()
        self.assertIn("Found UIDs: 12345, 67890", output)
        self.assertIn(self.mock_genbank_record_12345.strip(), output)
        self.assertIn(self.mock_genbank_record_67890.strip(), output)

        # Check search call
        mock_post.assert_called_once_with(
            query_ncbi_gi.SEARCH_URL,
            data={"cmd": "search", "term": self.search_term_ok, "db": "Nucleotide"},
            headers={"User-Agent": query_ncbi_gi.USER_AGENT},
        )
        # Check fetch calls
        self.assertEqual(mock_get.call_count, 2)
        expected_fetch_calls = [
            call(
                query_ncbi_gi.EFETCH_URL,
                params={
                    "db": "nuccore",
                    "id": "12345",
                    "rettype": "gb",
                    "retmode": "text",
                },
                headers={"User-Agent": query_ncbi_gi.USER_AGENT},
            ),
            call(
                query_ncbi_gi.EFETCH_URL,
                params={
                    "db": "nuccore",
                    "id": "67890",
                    "rettype": "gb",
                    "retmode": "text",
                },
                headers={"User-Agent": query_ncbi_gi.USER_AGENT},
            ),
        ]
        mock_get.assert_has_calls(expected_fetch_calls)

        # Check sleeps: 1s between fetches, 5s at the end
        expected_sleep_calls = [call(1), call(5)]
        self.assertEqual(mock_sleep.call_args_list, expected_sleep_calls)

    @patch("query_ncbi_gi.requests.post")
    @patch("query_ncbi_gi.requests.get")
    @patch("query_ncbi_gi.time.sleep")
    def test_search_successful_single_uid(self, mock_sleep, mock_get, mock_post):
        mock_post.return_value = self._mock_response(
            text_data=self.mock_search_html_single_uid
        )
        mock_get.return_value = self._mock_response(
            text_data=self.mock_genbank_record_12345
        )

        test_args = ["query_ncbi_gi.py", self.search_term_ok]
        captured_output = io.StringIO()
        with redirect_stdout(captured_output):
            with patch.object(sys, "argv", test_args):
                query_ncbi_gi.main()

        output = captured_output.getvalue()
        self.assertIn("Found UIDs: 12345", output)
        self.assertIn(self.mock_genbank_record_12345.strip(), output)
        mock_get.assert_called_once()
        # Only final sleep should be called as there's no loop for multiple fetches
        mock_sleep.assert_called_once_with(5)

    @patch("query_ncbi_gi.requests.post")
    @patch("query_ncbi_gi.requests.get")
    @patch("query_ncbi_gi.time.sleep")
    def test_search_no_uids_found(self, mock_sleep, mock_get, mock_post):
        mock_post.return_value = self._mock_response(
            text_data=self.mock_search_html_no_uids
        )

        test_args = ["query_ncbi_gi.py", self.search_term_no_results]
        captured_output = io.StringIO()
        with redirect_stdout(captured_output):
            with patch.object(sys, "argv", test_args):
                query_ncbi_gi.main()

        output = captured_output.getvalue()
        self.assertIn(
            f"No UIDs found for search term: {self.search_term_no_results}", output
        )
        mock_get.assert_not_called()  # Fetch should not be called
        mock_sleep.assert_called_once_with(5)  # Only final sleep

    # --- Test Fetch Workflow ---
    @patch("query_ncbi_gi.requests.post")
    @patch("query_ncbi_gi.requests.get")
    @patch("query_ncbi_gi.time.sleep")
    def test_fetch_invalid_genbank_record(self, mock_sleep, mock_get, mock_post):
        mock_post.return_value = self._mock_response(
            text_data=self.mock_search_html_single_uid
        )  # Search finds one UID
        mock_get.return_value = self._mock_response(
            text_data=self.mock_invalid_genbank_record
        )  # Fetch returns invalid data

        test_args = ["query_ncbi_gi.py", self.search_term_ok]
        captured_output = io.StringIO()
        with redirect_stdout(captured_output):
            with patch.object(sys, "argv", test_args):
                query_ncbi_gi.main()

        output = captured_output.getvalue()
        self.assertIn("Found UIDs: 12345", output)
        self.assertIn("Fetching GenBank record for UID: 12345", output)
        self.assertIn("LOCUS line not found in the record.", output)
        self.assertNotIn(
            "DEFINITION", output
        )  # Content of invalid record should not be printed as GenBank
        mock_sleep.assert_called_once_with(5)

    @patch("query_ncbi_gi.requests.post")
    @patch("query_ncbi_gi.requests.get")
    @patch("query_ncbi_gi.time.sleep")
    def test_fetch_empty_genbank_record(self, mock_sleep, mock_get, mock_post):
        mock_post.return_value = self._mock_response(
            text_data=self.mock_search_html_single_uid
        )
        mock_get.return_value = self._mock_response(
            text_data=self.mock_empty_genbank_record
        )

        test_args = ["query_ncbi_gi.py", self.search_term_ok]
        captured_output = io.StringIO()
        with redirect_stdout(captured_output):
            with patch.object(sys, "argv", test_args):
                query_ncbi_gi.main()

        output = captured_output.getvalue()
        # If record is empty, main() prints "Could not retrieve..."
        self.assertIn("Could not retrieve GenBank record for UID: 12345", output)
        mock_sleep.assert_called_once_with(5)

    # --- Test Error Handling ---
    @patch("query_ncbi_gi.requests.post")
    @patch("query_ncbi_gi.requests.get")
    @patch("query_ncbi_gi.time.sleep")
    def test_search_http_error(self, mock_sleep, mock_get, mock_post):
        http_error = requests.exceptions.RequestException("Search HTTP Error")
        mock_post.return_value = self._mock_response(
            status_code=500, raise_for_status_exception=http_error
        )

        test_args = ["query_ncbi_gi.py", self.search_term_ok]
        captured_output = io.StringIO()
        with redirect_stdout(captured_output):
            with patch.object(sys, "argv", test_args):
                query_ncbi_gi.main()

        output = captured_output.getvalue()
        self.assertIn("An HTTP error occurred: Search HTTP Error", output)
        mock_get.assert_not_called()  # Fetch should not be attempted
        mock_sleep.assert_called_once_with(5)  # Final sleep should still run

    @patch("query_ncbi_gi.requests.post")
    @patch("query_ncbi_gi.requests.get")
    @patch("query_ncbi_gi.time.sleep")
    def test_fetch_http_error_for_one_uid(self, mock_sleep, mock_get, mock_post):
        mock_post.return_value = self._mock_response(
            text_data=self.mock_search_html_multiple_uids
        )  # Two UIDs: 12345, 67890
        http_error_fetch = requests.exceptions.RequestException("Fetch HTTP Error")
        mock_get.side_effect = [
            self._mock_response(
                text_data=self.mock_genbank_record_12345
            ),  # First UID fetches successfully
            self._mock_response(
                status_code=500, raise_for_status_exception=http_error_fetch
            ),  # Second UID fetch fails
        ]

        test_args = ["query_ncbi_gi.py", self.search_term_ok]
        captured_output = io.StringIO()
        with redirect_stdout(captured_output):
            with patch.object(sys, "argv", test_args):
                query_ncbi_gi.main()

        output = captured_output.getvalue()
        self.assertIn("Found UIDs: 12345, 67890", output)
        self.assertIn(
            self.mock_genbank_record_12345.strip(), output
        )  # First record printed
        self.assertIn("Fetching GenBank record for UID: 67890", output)
        # The script currently exits on the first fetch error within the loop due to the outer try-except.
        self.assertIn("An HTTP error occurred: Fetch HTTP Error", output)
        self.assertNotIn(
            self.mock_genbank_record_67890.strip(), output
        )  # Second record not printed

        # Check sleeps: call(1) for the delay before the failing fetch, then call(5) in finally.
        expected_sleep_calls = [call(1), call(5)]
        self.assertEqual(mock_sleep.call_args_list, expected_sleep_calls)

    # --- Test Delays ---
    # test_successful_blast_workflow already tests delays for multiple UIDs.
    # test_argument_parsing, test_search_no_uids_found etc. test final delay for single/no UID cases.
    # Explicitly test the final sleep in all paths.
    @patch(
        "query_ncbi_gi.search_ncbi",
        side_effect=Exception("General error to stop processing"),
    )
    @patch("query_ncbi_gi.time.sleep")
    def test_final_sleep_runs_on_general_exception(self, mock_sleep, mock_search_ncbi):
        test_args = ["query_ncbi_gi.py", "any_term"]
        with redirect_stdout(io.StringIO()):  # Suppress other prints
            with patch.object(sys, "argv", test_args):
                query_ncbi_gi.main()
        mock_sleep.assert_any_call(5)


if __name__ == "__main__":
    unittest.main(argv=["first-arg-is-ignored"], exit=False) # Ensure single newline at EOF
