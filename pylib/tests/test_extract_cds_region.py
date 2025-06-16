import unittest
import subprocess
import os
import sys
import tempfile

# Ensure the script can find pylib
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.insert(0, project_root)

# Path to the script to be tested
SCRIPT_PATH = os.path.join(project_root, "pylib", "scripts", "extract_cds_region.py")

class TestExtractCDSRegion(unittest.TestCase):

    def run_script(self, gb_content, start_pos, end_pos):
        """Helper function to run the script and return its output."""
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".gb") as tmp_gb_file:
            tmp_gb_file.write(gb_content)
            tmp_gb_file_path = tmp_gb_file.name

        command = [
            sys.executable, SCRIPT_PATH,
            tmp_gb_file_path,
            str(start_pos),
            str(end_pos)
        ]

        process = subprocess.run(command, capture_output=True, text=True)

        os.remove(tmp_gb_file_path)
        return process.stdout, process.stderr

    def test_simple_cds_forward_strand_full_overlap(self):
        """Test a simple CDS on the forward strand, query fully overlaps CDS."""
        genbank_content = """
LOCUS       TestEntry                 200 bp    DNA     linear   SYN 01-JAN-2023
DEFINITION  Test sequence.
ACCESSION   TestEntry
VERSION     TestEntry.1
KEYWORDS    .
SOURCE      synthetic DNA construct
  ORGANISM  synthetic DNA construct
            .
FEATURES             Location/Qualifiers
     source          1..200
                     /organism="synthetic DNA construct"
                     /mol_type="genomic DNA"
     CDS             50..149
                     /codon_start=1
                     /transl_table=1
                     /protein_id="XYZ_001"
                     /product="hypothetical protein"
ORIGIN
        1 gatcagatca gatcagatca gatcagatca gatcagatca gatcagatca
       51 atgaagcttg aatctctatt cacaatttca tcaacatatc agtgtaagct
      101 agattcagta acaggctcag ctaagaaaca atatgagcaa taagatcaga
      151 tcaagatcag atcagatcag atcagatcag atcagatcag atcagatca
//
"""
        # Sequence for CDS (50-149): ATGAAGCTTGAATCTCTATTCACAATTTCATCAACATATCAGTGTAAGCTAGATTCAGTAACAGGCTCAGCTAAGAAACAATATGAGCAATAA
        # Translation: M K L N S L F T I S S T Y Q C K L D S V T G S A K K Q Y E Q Stop

        stdout, stderr = self.run_script(genbank_content, 50, 149) # Query exact CDS boundaries

        # Debugging: Print stdout and stderr if the test fails
        if stderr != "" or "Protein sequence from CDS XYZ_001 (frame adj: 0, table: 1): MKLNSLFTISSTYQCKLDSVTGSAKKQYEQ" not in stdout:
            print("\nDEBUG: test_simple_cds_forward_strand_full_overlap")
            print("STDOUT:\n", stdout)
            print("STDERR:\n", stderr)

        self.assertEqual(stderr, "") # Should be no errors
        self.assertIn("Processing record: TestEntry", stdout)
        # The full query range sequence
        self.assertIn("Nucleotide sequence in query range (50-149): ATGAAGCTTGAATCTCTATTCACAATTTCATCAACATATCAGTGTAAGCTAGATTCAGTAACAGGCTCAGCTAAGAAACAATATGAGCAATAA", stdout)
        self.assertIn("Found CDS: XYZ_001 located at [49:149]", stdout) # Biopython location [0-based_start:0-based_exclusive_end]
        # The sequence extracted from the CDS part that overlaps the query
        self.assertIn("Nucleotide sequence from CDS (XYZ_001) overlapping query: ATGAAGCTTGAATCTCTATTCACAATTTCATCAACATATCAGTGTAAGCTAGATTCAGTAACAGGCTCAGCTAAGAAACAATATGAGCAATAA", stdout)
        # The translated protein
        self.assertIn("Protein sequence from CDS XYZ_001 (frame adj: 0, table: 1): MKLNSLFTISSTYQCKLDSVTGSAKKQYEQ", stdout) # Stop is not included by to_stop=True

    def test_cds_on_reverse_strand(self):
        genbank_content = """
LOCUS       TestReverse              200 bp    DNA     linear   SYN 01-JAN-2023
DEFINITION  Test sequence with reverse strand CDS.
ACCESSION   TestReverse
VERSION     TestReverse.1
FEATURES             Location/Qualifiers
     source          1..200
     CDS             complement(50..149)
                     /codon_start=1
                     /transl_table=1
                     /protein_id="XYZ_002"
                     /product="hypothetical protein rev"
ORIGIN
        1 gatcagatca gatcagatca gatcagatca gatcagatca gatcagatca
       51 ttatgctcat attgtttctt agctgagcct gttactgaat ctagcttaca
      101 ctgatatgtt gatgaaatgt gaataagaga ttcaagcttc atgatcagat
      151 cagatcagat cagatcagat cagatcagatca gatcagatca gatcagatca
//
"""
        # Genomic 50..149 (0-based 49..148): TTATGCTCATATTGTTTCTTAGCTGAGCCTGTTACTGAATCTAGCTTACACTGATATGTTGATGAAATGTGAATAAGAGATTCAAGCTTCAT
        # Reverse Complement for 50-149: ATGAAGCTTGAATCTCTTATTCACATTTCATCAACATATCAGTGTAAGCTAGATTCAGTAACAGGCTCAGCTAAGAAACAATATGAGCATAA
        # Translation: M K L E S L I H I S S T Y Q C K L D S V T G S A K K Q Y E H Stop

        stdout, stderr = self.run_script(genbank_content, 50, 149)

        if stderr != "" or "Protein sequence from CDS XYZ_002 (frame adj: 0, table: 1): MKLESLIHISSTYQCKLDSVTGSAKKQYEH" not in stdout:
            print("\nDEBUG: test_cds_on_reverse_strand")
            print("STDOUT:\n", stdout)
            print("STDERR:\n", stderr)

        self.assertEqual(stderr, "")
        self.assertIn("Processing record: TestReverse", stdout)
        self.assertIn("Nucleotide sequence in query range (50-149): TTATGCTCATATTGTTTCTTAGCTGAGCCTGTTACTGAATCTAGCTTACACTGATATGTTGATGAAATGTGAATAAGAGATTCAAGCTTCAT", stdout)
        self.assertIn("Found CDS: XYZ_002 located at complement([49:149])", stdout)
        self.assertIn("Nucleotide sequence from CDS (XYZ_002) overlapping query: ATGAAGCTTGAATCTCTTATTCACATTTCATCAACATATCAGTGTAAGCTAGATTCAGTAACAGGCTCAGCTAAGAAACAATATGAGCATAA", stdout)
        self.assertIn("Protein sequence from CDS XYZ_002 (frame adj: 0, table: 1): MKLESLIHISSTYQCKLDSVTGSAKKQYEH", stdout)

    def test_query_partial_overlap_cds_start(self):
        """Query overlaps the beginning of a CDS."""
        genbank_content = """
LOCUS       TestPartialStart          200 bp    DNA     linear   SYN 01-JAN-2023
DEFINITION  Test sequence.
ACCESSION   TestPartialStart
FEATURES             Location/Qualifiers
     source          1..200
     CDS             50..149
                     /codon_start=1
                     /transl_table=1
                     /protein_id="XYZ_003"
ORIGIN
        1 gatcagatca gatcagatca gatcagatca gatcagatca gatcagatca
       51 atgaagcttg aatctctatt cacaatttca tcaacatatc agtgtaagct
      101 agattcagta acaggctcag ctaagaaaca atatgagcaa taagatcaga
      151 tcaagatcag atcagatcag atcagatcag atcagatcag atcagatca
// CDS (50-149): ATGAAGCTTGAATCTCTATTCACAATTTCATCAACATATCAGTGTAAGCTAGATTCAGTAACAGGCTCAGCTAAGAAACAATATGAGCAATAA
// Protein: M K L N S L F T I S S T Y Q C K L D S V T G S A K K Q Y E Q Stop
"""
        # Query 30-74 (1-based). This is 0-based 29 to 74(exclusive for query_end_0based).
        # CDS is 50-149 (1-based). This is 0-based 49 to 149(exclusive for feature.location.end in some contexts, but inclusive for BioPython FeatureLocation.end).
        # The script converts query to 0-based: query_start_0based = 29, query_end_0based = 74.
        # Record sequence from 29 to 73: gatcagatcaATGAAGCTTGAATCTCTATTCACAATTTCATCAAC (length 45)
        # CDS feature runs from record index 49 to 148.
        # Overlap on record: indices 49 to 73. (start=max(29,49)=49, end=min(74,149)=74).
        # This corresponds to bases 50 through 74 on the 1-based record.
        # Sequence of this overlap on record: ATGAAGCTTGAATCTCTATTCACAATTTCATCAAC (25 bases)
        # This is the first 25 bases of the CDS.
        # full_cds_sequence starts with ATGAAGCTT...
        # indices_in_full_cds_for_query will be [0, 1, ..., 24]
        # first_base_offset_in_cds = 0. codon_start_qualifier = 1.
        # effective_frame_start = ((1-1) + 0) % 3 = 0.
        # Sequence for translation: ATGAAGCTTGAATCTCTATTCACAATTTCATCAAC (25 bases)
        # Protein: M K L N S L F T I (25 bases // 3 = 8 codons, plus one base 'I' from 'ATC') -> MKLNSLFTI

        stdout, stderr = self.run_script(genbank_content, 30, 74)

        if stderr != "" or "Protein sequence from CDS XYZ_003 (frame adj: 0, table: 1): MKLNSLFTI" not in stdout:
            print("\nDEBUG: test_query_partial_overlap_cds_start")
            print("STDOUT:\n", stdout)
            print("STDERR:\n", stderr)

        self.assertEqual(stderr, "")
        self.assertIn("Nucleotide sequence in query range (30-74): gatcagatcaATGAAGCTTGAATCTCTATTCACAATTTCATCAAC", stdout)
        self.assertIn("Nucleotide sequence from CDS (XYZ_003) overlapping query: ATGAAGCTTGAATCTCTATTCACAATTTCATCAAC", stdout)
        self.assertIn("Protein sequence from CDS XYZ_003 (frame adj: 0, table: 1): MKLNSLFTI", stdout)

    def test_no_cds_in_region(self):
        """Test a query region that does not overlap with any CDS."""
        genbank_content = """
LOCUS       TestNoCDS                 200 bp    DNA     linear   SYN 01-JAN-2023
DEFINITION  Test sequence with no CDS in query.
ACCESSION   TestNoCDS
FEATURES             Location/Qualifiers
     source          1..200
     CDS             150..180
                     /codon_start=1
                     /transl_table=1
                     /protein_id="XYZ_004"
ORIGIN
        1 gatcagatca gatcagatca gatcagatca gatcagatca gatcagatca
       51 gatcagatca gatcagatca gatcagatca gatcagatca gatcagatca
      101 gatcagatca gatcagatca gatcagatca gatcagatca gatcagatca
      151 atgagcagatcagatcagatcagatcagatcagatcagatcagatcagatca
//
"""
        stdout, stderr = self.run_script(genbank_content, 10, 40)

        if "No CDS features found whose exonic parts overlap with the query range 10-40." not in stdout:
            print("\nDEBUG: test_no_cds_in_region")
            print("STDOUT:\n", stdout)
            print("STDERR:\n", stderr)

        self.assertEqual(stderr, "")
        self.assertIn("Nucleotide sequence in query range (10-40): gatcagatcagatcagatcagatcagatca", stdout)
        # This message changed in the script slightly. Let's check for the key part.
        self.assertIn("No CDS features found whose exonic parts overlap", stdout)
        self.assertNotIn("Protein sequence from CDS", stdout)

    def test_invalid_range_start_greater_than_end(self):
        """Test invalid input: start position is greater than end position."""
        genbank_content = """
LOCUS       TestInvalidRange          100 bp    DNA     linear   SYN 01-JAN-2023
DEFINITION  Minimal test sequence.
ACCESSION   TestInvalidRange
FEATURES             Location/Qualifiers
     source          1..100
ORIGIN
        1 gatcagatca gatcagatca gatcagatca gatcagatca gatcagatca
       51 gatcagatca gatcagatca gatcagatca gatcagatca gatcagatca
//
"""
        stdout, stderr = self.run_script(genbank_content, 50, 20)
        self.assertIn("Error: Start position (50) cannot be greater than end position (20).", stderr)
        self.assertEqual(stdout, "")

    def test_invalid_range_out_of_bounds(self):
        """Test invalid input: range is out of sequence bounds."""
        genbank_content = """
LOCUS       TestOutOfBounds           100 bp    DNA     linear   SYN 01-JAN-2023
DEFINITION  Minimal test sequence.
ACCESSION   TestOutOfBounds
FEATURES             Location/Qualifiers
     source          1..100
ORIGIN
        1 gatcagatca gatcagatca gatcagatca gatcagatca gatcagatca
       51 gatcagatca gatcagatca gatcagatca gatcagatca gatcagatca
//
"""
        stdout, stderr = self.run_script(genbank_content, 80, 120)
        # The script prints record info, then the error.
        self.assertIn("Processing record: TestOutOfBounds", stdout)
        self.assertIn("Error: Query range 80-120 is out of bounds for record TestOutOfBounds (length 100).", stderr)


if __name__ == '__main__':
    unittest.main()
