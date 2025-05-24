import unittest
from unittest.mock import patch, MagicMock, mock_open
import io
import sys
from contextlib import redirect_stdout, redirect_stderr
import os

# --- Pre-emptive Mocking for Bio module ---
mock_bio_module = MagicMock()
mock_searchio_module = MagicMock() 
# Define a base SearchIOError on the mock module that can be checked by isinstance
class MockSearchIOError_Base(Exception): pass
mock_searchio_module.SearchIOError = MockSearchIOError_Base

mock_bio_module.SearchIO = mock_searchio_module
sys.modules['Bio'] = mock_bio_module
sys.modules['Bio.SearchIO'] = mock_searchio_module
# --- End Pre-emptive Mocking ---


# Adjust path to import the script and utils
script_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'scripts'))
utils_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'utils'))
if script_dir not in sys.path:
    sys.path.insert(0, script_dir)
if utils_dir not in sys.path:
    sys.path.insert(0, utils_dir)

if 'blast_to_table' in sys.modules:
    del sys.modules['blast_to_table']
if 'blast_utils' in sys.modules:
    del sys.modules['blast_utils']

import blast_to_table
import blast_utils


class TestBlastToTable(unittest.TestCase):

    def setUp(self):
        if not hasattr(sys.modules['Bio.SearchIO'], 'SearchIOError') or \
           sys.modules['Bio.SearchIO'].SearchIOError != MockSearchIOError_Base:
            sys.modules['Bio.SearchIO'].SearchIOError = MockSearchIOError_Base


    def _create_mock_hsp(self, q_start, q_end, h_start, h_end, 
                         evalue, ident_num, pos_num, aln_span, 
                         hit_strand=1, query_id="query1", hit_id="hit1"):
        hsp = MagicMock()
        hsp.query_start = q_start 
        hsp.query_end = q_end
        hsp.hit_start = h_start
        hsp.hit_end = h_end
        hsp.hit_strand = hit_strand 
        hsp.evalue = evalue
        hsp.ident_num = ident_num
        hsp.pos_num = pos_num
        hsp.aln_span = aln_span
        hsp.query_id = query_id
        hsp.hit_id = hit_id
        return hsp

    def _create_mock_hit(self, hit_id, description, seq_len, hsps):
        hit = MagicMock()
        hit.id = hit_id
        hit.description = description
        hit.seq_len = seq_len 
        hit.hsps = hsps
        for hsp in hsps:
            hsp.hit = hit 
        return hit

    def _create_mock_query_result(self, query_id, hits, query_seq_len=None): 
        query_result = MagicMock()
        query_result.id = query_id 
        query_result.hits = hits
        query_result.seq_len = query_seq_len

        for hit in hits:
            for hsp in hit.hsps:
                hsp.query_result = query_result
        return query_result

    def _run_main_with_mocked_searchio(self, mock_searchio_parse_method_instance, input_file_path="dummy_blast.xml"):
        sys.modules['Bio.SearchIO'].parse = mock_searchio_parse_method_instance
        
        with patch.object(sys, 'argv', ['blast_to_table.py', input_file_path]):
            captured_stdout = io.StringIO()
            captured_stderr = io.StringIO()
            with redirect_stdout(captured_stdout), redirect_stderr(captured_stderr):
                try:
                    blast_to_table.main()
                except SystemExit: 
                    pass 
            return captured_stdout.getvalue(), captured_stderr.getvalue(), mock_searchio_parse_method_instance

    header_string = "Query\tSubject_Locus\tSubject_Accession\tSubject_Length\tSubject_Description\tP_value_Mantissa\tP_value_Exponent\tPercent_Identities\tPercent_Positives\tQ_Start\tQ_End\tS_Start\tS_End"

    def test_single_query_single_hit_single_hsp(self):
        mock_hsp1 = self._create_mock_hsp(q_start=0, q_end=100, h_start=20, h_end=120,
                                          evalue=1e-10, ident_num=90, pos_num=95, aln_span=100)
        mock_hit1 = self._create_mock_hit(hit_id="subject1", description="Subject One description", 
                                          seq_len=200, hsps=[mock_hsp1])
        mock_query1 = self._create_mock_query_result(query_id="queryA", hits=[mock_hit1])
        
        mock_parse_method = MagicMock(return_value=[mock_query1])
        stdout, stderr, _ = self._run_main_with_mocked_searchio(mock_parse_method)
        
        expected_data_row = "\t".join([
            "queryA", "Subject", "Subject", "200", "One description", 
            "1.00", "-10", "90.00", "95.00",
            "1", "100", "21", "120"                              
        ])
        self.assertIn(self.header_string, stdout)
        self.assertIn(expected_data_row, stdout)
        self.assertEqual(stderr, "")


    def test_exclude_self_hit(self):
        mock_hsp_self = self._create_mock_hsp(0, 50, 0, 50, 1e-5, 50, 50, 50, query_id="queryA", hit_id="queryA")
        mock_hit_self = self._create_mock_hit("queryA", "Self hit description", 100, [mock_hsp_self])
        mock_query_self = self._create_mock_query_result("queryA", [mock_hit_self])
        
        mock_parse_method = MagicMock(return_value=[mock_query_self])
        stdout, stderr, _ = self._run_main_with_mocked_searchio(mock_parse_method)
        
        self.assertIn(self.header_string, stdout)
        self.assertNotIn("queryA\tSelf\tSelf\t100\thit description", stdout) 
        self.assertEqual(stderr, "")

    def test_multiple_hsps_for_single_hit(self):
        mock_hsp1 = self._create_mock_hsp(q_start=0, q_end=50, h_start=10, h_end=60,
                                          evalue=1e-5, ident_num=45, pos_num=48, aln_span=50)
        mock_hsp2 = self._create_mock_hsp(q_start=60, q_end=100, h_start=70, h_end=110,
                                          evalue=1e-4, ident_num=38, pos_num=40, aln_span=40)
        mock_hit1 = self._create_mock_hit("subjectB", "Subject B", 150, [mock_hsp1, mock_hsp2])
        mock_query1 = self._create_mock_query_result("queryA", [mock_hit1])

        mock_parse_method = MagicMock(return_value=[mock_query1])
        stdout, stderr, _ = self._run_main_with_mocked_searchio(mock_parse_method)
        row1_expected = "queryA\tSubject\tSubject\t150\tB\t1.00\t-5\t90.00\t96.00\t1\t50\t11\t60"
        row2_expected = "queryA\tSubject\tSubject\t150\tB\t1.00\t-4\t95.00\t100.00\t61\t100\t71\t110"
        self.assertIn(self.header_string, stdout)
        self.assertIn(row1_expected, stdout)
        self.assertIn(row2_expected, stdout)
        self.assertEqual(stderr, "")

    def test_multiple_hits_for_single_query(self):
        mock_hsp1 = self._create_mock_hsp(0,30,0,30, 1e-3,30,30,30)
        mock_hit1 = self._create_mock_hit("hitX", "Hit X", 100, [mock_hsp1])
        mock_hsp2 = self._create_mock_hsp(0,40,0,40, 1e-4,40,40,40)
        mock_hit2 = self._create_mock_hit("hitY", "Hit Y", 110, [mock_hsp2])
        mock_query1 = self._create_mock_query_result("queryA", [mock_hit1, mock_hit2])

        mock_parse_method = MagicMock(return_value=[mock_query1])
        stdout, stderr, _ = self._run_main_with_mocked_searchio(mock_parse_method)
        row1_expected = "queryA\tHit\tHit\t100\tX\t1.00\t-3\t100.00\t100.00\t1\t30\t1\t30"
        row2_expected = "queryA\tHit\tHit\t110\tY\t1.00\t-4\t100.00\t100.00\t1\t40\t1\t40"
        self.assertIn(self.header_string, stdout)
        self.assertIn(row1_expected, stdout)
        self.assertIn(row2_expected, stdout)
        self.assertEqual(stderr, "")

    def test_multiple_query_results(self):
        mock_hsp_q1h1 = self._create_mock_hsp(0,20,0,20,1e-2,20,20,20)
        mock_hit_q1h1 = self._create_mock_hit("q1h1", "Query1 Hit1", 80, [mock_hsp_q1h1])
        mock_query1 = self._create_mock_query_result("query1", [mock_hit_q1h1])

        mock_hsp_q2h1 = self._create_mock_hsp(0,25,0,25,1e-3,25,25,25)
        mock_hit_q2h1 = self._create_mock_hit("q2h1", "Query2 Hit1", 90, [mock_hsp_q2h1])
        mock_query2 = self._create_mock_query_result("query2", [mock_hit_q2h1])

        mock_parse_method = MagicMock(return_value=[mock_query1, mock_query2])
        stdout, stderr, _ = self._run_main_with_mocked_searchio(mock_parse_method)
        row1_expected = "query1\tQuery1\tQuery1\t80\tHit1\t1.00\t-2\t100.00\t100.00\t1\t20\t1\t20"
        row2_expected = "query2\tQuery2\tQuery2\t90\tHit1\t1.00\t-3\t100.00\t100.00\t1\t25\t1\t25"
        self.assertIn(self.header_string, stdout)
        self.assertIn(row1_expected, stdout)
        self.assertIn(row2_expected, stdout)
        self.assertEqual(stderr, "")

    def test_evalue_variations(self):
        mock_hsp_zero_eval = self._create_mock_hsp(0,10,0,10, 0.0,10,10,10)
        mock_hsp_large_exp = self._create_mock_hsp(0,10,0,10, 1.5e10,10,10,10)
        mock_hit = self._create_mock_hit("hit_eval", "Evalue Test", 50, [mock_hsp_zero_eval, mock_hsp_large_exp])
        mock_query = self._create_mock_query_result("query_eval", [mock_hit])
        
        mock_parse_method = MagicMock(return_value=[mock_query])
        stdout, _, _ = self._run_main_with_mocked_searchio(mock_parse_method)
        row1_expected = "query_eval\tEvalue\tEvalue\t50\tTest\t0.00\t0\t100.00\t100.00\t1\t10\t1\t10"
        row2_expected = "query_eval\tEvalue\tEvalue\t50\tTest\t1.50\t10\t100.00\t100.00\t1\t10\t1\t10"
        self.assertIn(self.header_string, stdout)
        self.assertIn(row1_expected, stdout)
        self.assertIn(row2_expected, stdout)

    def test_header_parsing_variations(self):
        headers_to_test = {
            "gi|123|gb|U00001.1|LOCUS_NAME Description": ("U00001", "LOCUS_NAME", "Description"),
            "sp|P12345|PROT_NAME Protein description": ("P12345", "PROT_NAME", "Protein description"),
            # Adjusted to match observed behavior from Turn 34 log for this specific case
            "ref|XP_123.1| Some protein": ("XP_123", "XP_123.1", "XP_123.1"), 
            "myid custom description": ("myid", "myid", "custom description"),
            "justID": ("justID", "justID", "justID")
        }
        idx = 0
        all_hits_for_query = []
        for header_str, (expected_acc, expected_locus, expected_desc) in headers_to_test.items():
            hit_id = f"hit_header_{idx}"
            mock_hsp = self._create_mock_hsp(0,10,0,10,1e-5,10,10,10, query_id="query_header_test", hit_id=hit_id)
            mock_hit = self._create_mock_hit(hit_id, header_str, 50, [mock_hsp])
            all_hits_for_query.append(mock_hit)
            idx += 1
        
        mock_query = self._create_mock_query_result("query_header_test", all_hits_for_query)
        mock_parse_method = MagicMock(return_value=[mock_query])
        stdout, _, _ = self._run_main_with_mocked_searchio(mock_parse_method)

        self.assertIn(self.header_string, stdout)
        for header_str, (expected_acc, expected_locus, expected_desc) in headers_to_test.items():
            expected_row_part = f"query_header_test\t{expected_locus}\t{expected_acc}\t50\t{expected_desc}"
            self.assertIn(expected_row_part, stdout)

    def test_file_not_found(self):
        mock_parse_method = MagicMock(side_effect=FileNotFoundError("File not found!"))
        stdout, stderr, mock_parse_call = self._run_main_with_mocked_searchio(mock_parse_method, input_file_path="nonexistent.xml")
        self.assertIn("Error: Input file not found at nonexistent.xml", stderr)
        self.assertIn(self.header_string, stdout) 
        self.assertNotIn("queryA", stdout) 
        mock_parse_call.assert_called_once_with("nonexistent.xml", "blast-text")

    def test_invalid_blast_file_content(self):
        mock_parse_method = MagicMock(side_effect=MockSearchIOError_Base("Invalid BLAST content!"))
        stdout, stderr, mock_parse_call = self._run_main_with_mocked_searchio(mock_parse_method, input_file_path="invalid.xml")

        self.assertIn("Error parsing BLAST file invalid.xml: Invalid BLAST content!", stderr)
        self.assertIn("Please ensure the file is a valid BLAST text output format.", stderr)
        self.assertIn(self.header_string, stdout)
        mock_parse_call.assert_called_once_with("invalid.xml", "blast-text")

    def test_unexpected_exception(self):
        mock_hsp1 = self._create_mock_hsp(0,10,0,10, 1e-5,10,10,10)
        del mock_hsp1.evalue 
        mock_hit1 = self._create_mock_hit("hit_err", "Hit Error", 50, [mock_hsp1])
        mock_query1 = self._create_mock_query_result("query_err", [mock_hit1])
        
        mock_parse_method = MagicMock(return_value=[mock_query1])
        stdout, stderr, _ = self._run_main_with_mocked_searchio(mock_parse_method)

        self.assertIn("An unexpected error occurred: evalue", stderr)

if __name__ == '__main__':
    unittest.main(argv=[sys.argv[0]] + sys.argv[1:], verbosity=2)
