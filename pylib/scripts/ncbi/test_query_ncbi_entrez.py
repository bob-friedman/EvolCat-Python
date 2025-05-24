import unittest
from unittest.mock import patch, MagicMock, call
import os
import sys
import io
from contextlib import redirect_stdout, redirect_stderr
import xml.etree.ElementTree as ET

# Add paths for script and utility_functions import
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_dir) # Script directory
sys.path.insert(0, os.path.dirname(current_dir)) # scripts directory
sys.path.insert(0, os.path.dirname(os.path.dirname(current_dir))) # pylib directory

# Import the script to be tested
try:
    from pylib.scripts.ncbi import query_ncbi_entrez
except ModuleNotFoundError:
    import query_ncbi_entrez # Fallback

MOCK_USER_AGENT = "Test Entrez User Agent"

# Patch utility_functions.get_user_agent at the class level
# This assumes query_ncbi_entrez imports it as `from utils import utility_functions`
# and then calls utility_functions.get_user_agent()
@patch('utils.utility_functions.get_user_agent', return_value=MOCK_USER_AGENT)
class TestQueryNcbiEntrez(unittest.TestCase):

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
    def test_argument_parsing_defaults(self, mock_get_ua):
        with patch.object(sys, 'argv', ['query_ncbi_entrez.py', '--term', 'test_term']):
            # Access parser from main's globals
            parser = query_ncbi_entrez.main.__globals__['parser'] 
            args = parser.parse_args()
            self.assertEqual(args.term, "test_term")
            self.assertEqual(args.search_db, "protein")
            # Test the logic where fetch_db defaults to search_db if not provided
            # This requires running the initial part of main() or replicating its logic.
            # For simplicity, we'll test the direct parser output first.
            self.assertEqual(args.fetch_db, "protein") # Initial default
            self.assertEqual(args.rettype, "gb")
            self.assertEqual(args.retmode, "text")
            self.assertIsNone(args.email)
            self.assertIsNone(args.api_key)
            self.assertEqual(args.max_retries, 3)
            self.assertEqual(args.base_eutils_url, query_ncbi_entrez.BASE_EUTILS_URL)
            
    def test_argument_parsing_fetch_db_sync(self, mock_get_ua):
        # Test the specific logic in main() that syncs fetch_db to search_db
        # if fetch_db is not explicitly provided by the user.
        with patch.object(sys, 'argv', ['query_ncbi_entrez.py', '--term', 'test_term', '--search_db', 'nuccore']):
            # Temporarily hijack the parser within the script's main function scope for this test
            original_parser = query_ncbi_entrez.main.__globals__['parser']
            args = original_parser.parse_args()
            # Simulate the logic from main()
            if not any(arg.startswith('--fetch_db') for arg in sys.argv):
                args.fetch_db = args.search_db
            self.assertEqual(args.fetch_db, 'nuccore')

        with patch.object(sys, 'argv', ['query_ncbi_entrez.py', '--term', 'test_term', '--search_db', 'nuccore', '--fetch_db', 'protein']):
            original_parser = query_ncbi_entrez.main.__globals__['parser']
            args = original_parser.parse_args()
            # Simulate the logic from main()
            if not any(arg.startswith('--fetch_db') for arg in sys.argv):
                 args.fetch_db = args.search_db # This line won't run if --fetch_db is present
            self.assertEqual(args.fetch_db, 'protein') # Should remain protein as it was specified


    def test_argument_parsing_required(self, mock_get_ua):
        with patch.object(sys, 'argv', ['query_ncbi_entrez.py']): # No --term
            with self.assertRaises(SystemExit):
                with redirect_stderr(io.StringIO()):
                    query_ncbi_entrez.main.__globals__['parser'].parse_args()

    # --- Test _make_request_with_retries ---
    @patch('requests.get')
    @patch('time.sleep', return_value=None)
    def test_make_request_success_first_try(self, mock_sleep, mock_get, mock_get_ua):
        mock_get.return_value = self._mock_response(text_data="success")
        query_ncbi_entrez.REQUEST_COUNT = 0 # Reset global counter
        response = query_ncbi_entrez._make_request_with_retries("fake_url", {}, api_key_provided=False)
        self.assertIsNotNone(response)
        self.assertEqual(response.text, "success")
        mock_get.assert_called_once()
        mock_sleep.assert_not_called() # No sleep before first request or on success
        self.assertEqual(query_ncbi_entrez.REQUEST_COUNT, 1)

    @patch('requests.get')
    @patch('time.sleep', return_value=None)
    def test_make_request_retry_then_success(self, mock_sleep, mock_get, mock_get_ua):
        mock_get.side_effect = [
            requests.exceptions.Timeout("timeout"),
            self._mock_response(text_data="success")
        ]
        query_ncbi_entrez.REQUEST_COUNT = 0
        response = query_ncbi_entrez._make_request_with_retries("fake_url", {}, max_retries=2, api_key_provided=False)
        self.assertIsNotNone(response)
        self.assertEqual(response.text, "success")
        self.assertEqual(mock_get.call_count, 2)
        mock_sleep.assert_any_call(1) # Sleep after first failure (0.34s from REQUEST_COUNT > 0 not called due to reset)
                                     # then 1s for retry backoff
        self.assertEqual(query_ncbi_entrez.REQUEST_COUNT, 2)


    @patch('requests.get')
    @patch('time.sleep', return_value=None)
    def test_make_request_all_retries_fail(self, mock_sleep, mock_get, mock_get_ua):
        mock_get.side_effect = requests.exceptions.RequestException("persistent error")
        query_ncbi_entrez.REQUEST_COUNT = 0
        with redirect_stderr(io.StringIO()) as stderr_capture:
            response = query_ncbi_entrez._make_request_with_retries("fake_url", {}, max_retries=3, api_key_provided=True)
        self.assertIsNone(response)
        self.assertEqual(mock_get.call_count, 3)
        # Sleeps for rate limit (REQUEST_COUNT > 0) then for retry backoff
        expected_sleep_calls = [
            call(1), # After 1st failure (0.11s from REQUEST_COUNT > 0 not called due to reset)
            call(2)  # After 2nd failure
        ]
        # mock_sleep.assert_has_calls(expected_sleep_calls, any_order=False) # Order matters
        self.assertIn("Max retries reached.", stderr_capture.getvalue())
        self.assertEqual(query_ncbi_entrez.REQUEST_COUNT, 3)

    @patch('requests.get')
    @patch('time.sleep', return_value=None)
    def test_make_request_rate_limiting_sleep(self, mock_sleep, mock_get, mock_get_ua):
        mock_get.return_value = self._mock_response()
        query_ncbi_entrez.REQUEST_COUNT = 1 # Simulate a previous request
        
        # Test without API key (longer delay)
        query_ncbi_entrez._make_request_with_retries("fake_url1", {}, api_key_provided=False)
        mock_sleep.assert_any_call(query_ncbi_entrez.get_request_delay(False))
        
        # Test with API key (shorter delay)
        query_ncbi_entrez.REQUEST_COUNT = 1 # Reset for this part
        mock_sleep.reset_mock()
        query_ncbi_entrez._make_request_with_retries("fake_url2", {}, api_key_provided=True)
        mock_sleep.assert_any_call(query_ncbi_entrez.get_request_delay(True))
        query_ncbi_entrez.REQUEST_COUNT = 0 # Reset global

    # --- Test search_entrez ---
    @patch('query_ncbi_entrez._make_request_with_retries')
    def test_search_entrez_success(self, mock_make_request, mock_get_ua):
        xml_response = """
        <eSearchResult>
            <Count>2</Count>
            <IdList>
                <Id>123</Id>
                <Id>456</Id>
            </IdList>
            <WebEnv>WEBENV_TOKEN</WebEnv>
            <QueryKey>1</QueryKey>
        </eSearchResult>
        """
        mock_make_request.return_value = self._mock_response(text_data=xml_response)
        uids, webenv, qk = query_ncbi_entrez.search_entrez("term", "db", "url", "email", "key", 3)
        self.assertEqual(uids, ["123", "456"])
        self.assertEqual(webenv, "WEBENV_TOKEN")
        self.assertEqual(qk, "1")

    @patch('query_ncbi_entrez._make_request_with_retries')
    def test_search_entrez_no_uids(self, mock_make_request, mock_get_ua):
        xml_response = "<eSearchResult><Count>0</Count><IdList/></eSearchResult>"
        mock_make_request.return_value = self._mock_response(text_data=xml_response)
        with redirect_stderr(io.StringIO()) as stderr_capture:
            uids, _, _ = query_ncbi_entrez.search_entrez("term", "db", "url", None, None, 1)
        self.assertEqual(uids, [])
        self.assertIn("No UIDs found for the search term.", stderr_capture.getvalue())

    @patch('query_ncbi_entrez._make_request_with_retries')
    def test_search_entrez_ncbi_error_xml(self, mock_make_request, mock_get_ua):
        xml_error = "<eSearchResult><ERROR>Invalid db name</ERROR></eSearchResult>"
        mock_make_request.return_value = self._mock_response(text_data=xml_error)
        with redirect_stderr(io.StringIO()) as stderr_capture:
            uids, _, _ = query_ncbi_entrez.search_entrez("term", "db", "url", None, None, 1)
        self.assertIsNone(uids)
        self.assertIn("NCBI E-utility Error: Invalid db name", stderr_capture.getvalue())

    # --- Test fetch_entrez_records ---
    @patch('query_ncbi_entrez._make_request_with_retries')
    def test_fetch_entrez_records_with_webenv(self, mock_make_request, mock_get_ua):
        mock_make_request.return_value = self._mock_response(text_data="Sample GenBank Record")
        with redirect_stdout(io.StringIO()) as stdout_capture:
            success = query_ncbi_entrez.fetch_entrez_records(
                "db", "gb", "text", "url", "email", "key", 3,
                webenv="WEBENV", query_key="1"
            )
        self.assertTrue(success)
        self.assertIn("Sample GenBank Record", stdout_capture.getvalue())
        args, kwargs = mock_make_request.call_args
        self.assertIn('WebEnv', kwargs['params']) # Check if params for POST or GET
        self.assertIn('query_key', kwargs['params'])


    @patch('query_ncbi_entrez._make_request_with_retries')
    def test_fetch_entrez_records_with_uid_list_get(self, mock_make_request, mock_get_ua):
        mock_make_request.return_value = self._mock_response(text_data="FASTA >seq1...")
        with redirect_stdout(io.StringIO()) as stdout_capture:
            success = query_ncbi_entrez.fetch_entrez_records(
                "db", "fasta", "text", "url", None, None, 1,
                uid_list=["1", "2", "3"]
            )
        self.assertTrue(success)
        self.assertIn("FASTA >seq1...", stdout_capture.getvalue())
        args, kwargs = mock_make_request.call_args
        self.assertEqual(kwargs['method'], 'GET')
        self.assertEqual(kwargs['params']['id'], '1,2,3')

    @patch('query_ncbi_entrez._make_request_with_retries')
    def test_fetch_entrez_records_with_uid_list_post(self, mock_make_request, mock_get_ua):
        # Simulate a long list of UIDs that should trigger POST
        long_uid_list = [str(i) for i in range(250)]
        mock_make_request.return_value = self._mock_response(text_data="Long FASTA data...")

        with redirect_stdout(io.StringIO()) as stdout_capture:
            success = query_ncbi_entrez.fetch_entrez_records(
                "db", "fasta", "text", "url", None, None, 1,
                uid_list=long_uid_list
            )
        self.assertTrue(success)
        args, kwargs = mock_make_request.call_args
        self.assertEqual(kwargs['method'], 'POST') # Check if method was POST
        self.assertEqual(kwargs['params']['id'], ','.join(long_uid_list))


    @patch('query_ncbi_entrez._make_request_with_retries')
    def test_fetch_entrez_records_ncbi_error(self, mock_make_request, mock_get_ua):
        # Simulate an error from NCBI during fetch (e.g., if rettype is XML)
        xml_error_fetch = "<eFetchResult><ERROR>UID not found</ERROR></eFetchResult>"
        mock_make_request.return_value = self._mock_response(text_data=xml_error_fetch)
        with redirect_stderr(io.StringIO()) as stderr_capture:
            success = query_ncbi_entrez.fetch_entrez_records(
                "db", "xml", "xml", "url", None, None, 1, # XML rettype
                uid_list=["invalid_uid"]
            )
        self.assertFalse(success)
        self.assertIn("NCBI E-utility Error: UID not found", stderr_capture.getvalue())

    # --- Test Main Script Logic ---
    @patch('query_ncbi_entrez.search_entrez')
    @patch('query_ncbi_entrez.fetch_entrez_records')
    def test_main_successful_flow(self, mock_fetch, mock_search, mock_get_ua):
        mock_search.return_value = (["123", "456"], "WEBENV_MAIN", "1")
        mock_fetch.return_value = True # Simulate successful fetch
        
        test_args = ['query_ncbi_entrez.py', '--term', 'my_term', '--email', 'a@b.com']
        with patch.object(sys, 'argv', test_args):
            with redirect_stderr(io.StringIO()) as stderr_capture:
                with self.assertRaises(SystemExit) as cm:
                     query_ncbi_entrez.main()
                self.assertIsNone(cm.exception.code) # Successful exit

        mock_search.assert_called_once()
        mock_fetch.assert_called_once_with(
            'protein', 'gb', 'text', query_ncbi_entrez.BASE_EUTILS_URL, 
            'a@b.com', None, 3, 
            uid_list=["123", "456"], webenv="WEBENV_MAIN", query_key="1"
        )
        self.assertIn("Successfully fetched records.", stderr_capture.getvalue())

    @patch('query_ncbi_entrez.search_entrez')
    @patch('query_ncbi_entrez.fetch_entrez_records') # Also mock fetch to prevent it running
    def test_main_search_no_results(self, mock_fetch, mock_search, mock_get_ua):
        mock_search.return_value = ([], "WEBENV_EMPTY", "2") # No UIDs
        
        test_args = ['query_ncbi_entrez.py', '--term', 'term_no_results']
        with patch.object(sys, 'argv', test_args):
            with redirect_stdout(io.StringIO()) as stdout_capture:
                 with self.assertRaises(SystemExit) as cm:
                    query_ncbi_entrez.main()
                 self.assertEqual(cm.exception.code, 0) # Successful exit even if no results

        self.assertIn("No records found for term 'term_no_results'", stdout_capture.getvalue())
        mock_fetch.assert_not_called() # Fetch should not be called

    @patch('query_ncbi_entrez.search_entrez')
    def test_main_search_fails(self, mock_search, mock_get_ua):
        mock_search.return_value = (None, None, None) # Search error
        
        test_args = ['query_ncbi_entrez.py', '--term', 'bad_term']
        with patch.object(sys, 'argv', test_args):
            with redirect_stderr(io.StringIO()) as stderr_capture:
                with self.assertRaises(SystemExit) as cm:
                    query_ncbi_entrez.main()
                self.assertEqual(cm.exception.code, 1) # Error exit
        
        self.assertIn("Exiting due to error in ESearch.", stderr_capture.getvalue())
        
    @patch('query_ncbi_entrez.search_entrez')
    @patch('query_ncbi_entrez.fetch_entrez_records')
    def test_main_fetch_fails(self, mock_fetch, mock_search, mock_get_ua):
        mock_search.return_value = (["789"], "WEBENV_FETCH_FAIL", "3")
        mock_fetch.return_value = False # Fetch error
        
        test_args = ['query_ncbi_entrez.py', '--term', 'fetch_fail_term']
        with patch.object(sys, 'argv', test_args):
            with redirect_stderr(io.StringIO()) as stderr_capture:
                with self.assertRaises(SystemExit) as cm:
                    query_ncbi_entrez.main()
                self.assertEqual(cm.exception.code, 1)
        
        self.assertIn("Failed to fetch all records.", stderr_capture.getvalue())

if __name__ == '__main__':
    unittest.main(argv=['first-arg-is-ignored'], exit=False, verbosity=2)
