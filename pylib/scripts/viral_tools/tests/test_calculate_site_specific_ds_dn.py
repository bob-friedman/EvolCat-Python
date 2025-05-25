import unittest
import subprocess
import os
import shutil
import tempfile
import sys

# Path to the script to be tested
SCRIPT_PATH = os.path.join(os.path.dirname(__file__), "..", "calculate_site_specific_ds_dn.py")
CODEML_AVAILABLE = shutil.which("codeml") is not None

# Example FASTA content
FASTA_CONTENT = """>Seq1
ATGCGTAGCATT
>Seq2
ATGCGTACCATT
>Seq3
ATGCGTTCGATT
""" # 12 bases = 4 codons

# Example Newick tree content
NEWICK_CONTENT = "(Seq1:0.1,Seq2:0.1,Seq3:0.1);"

class TestCalculateSiteSpecificDsDn(unittest.TestCase):

    def setUp(self):
        self.temp_dir_obj = tempfile.TemporaryDirectory()
        self.temp_dir_path = self.temp_dir_obj.name

        self.alignment_file = os.path.join(self.temp_dir_path, "test_alignment.fasta")
        with open(self.alignment_file, "w") as f:
            f.write(FASTA_CONTENT)

        self.tree_file = os.path.join(self.temp_dir_path, "test_tree.nwk")
        with open(self.tree_file, "w") as f:
            f.write(NEWICK_CONTENT)
            
        # Create a dummy non-executable file for testing codeml path failure
        self.dummy_non_exec_codeml = os.path.join(self.temp_dir_path, "dummy_codeml")
        with open(self.dummy_non_exec_codeml, "w") as f:
            f.write("#!/bin/sh\nexit 1") # Not executable by default
        # os.chmod(self.dummy_non_exec_codeml, 0o644) # Ensure it's not executable

    def tearDown(self):
        self.temp_dir_obj.cleanup()

    def _run_script(self, args_list):
        """Helper method to run the script and return its output."""
        full_args = [sys.executable, SCRIPT_PATH] + args_list # Use sys.executable for python interpreter
        return subprocess.run(full_args, capture_output=True, text=True, cwd=self.temp_dir_path)

    @unittest.skipIf(CODEML_AVAILABLE, "codeml is available, skipping test for its absence")
    def test_codeml_not_found_graceful_exit_dummy_path(self):
        """Test graceful exit if codeml is not found at the specified dummy path."""
        # This test runs if shutil.which("codeml") is None.
        # We provide a path to a non-executable file.
        args = [
            "--alignment", self.alignment_file,
            "--tree", self.tree_file,
            "--model", "M0",
            "--outfile_prefix", "test_run_m0_codeml_fail",
            "--paml_path", self.dummy_non_exec_codeml 
        ]
        process = self._run_script(args)
        self.assertNotEqual(process.returncode, 0, "Script should exit with non-zero status if dummy codeml is not executable.")
        self.assertIn(f"Error: codeml executable not found at '{self.dummy_non_exec_codeml}'", process.stdout)

    # A variation for when codeml is not in PATH and --paml_path is not given
    @unittest.skipIf(CODEML_AVAILABLE, "codeml is available, skipping test for its absence from PATH")
    def test_codeml_not_in_path_graceful_exit(self):
        """Test graceful exit if codeml is not in PATH and --paml_path is not given."""
        args = [
            "--alignment", self.alignment_file,
            "--tree", self.tree_file,
            "--model", "M0",
            "--outfile_prefix", "test_run_m0_codeml_path_fail"
            # No --paml_path here
        ]
        process = self._run_script(args)
        self.assertNotEqual(process.returncode, 0, "Script should exit with non-zero status if codeml is not in PATH.")
        # shutil.which uses os.environ.get("PATH"), so this should trigger the "not found"
        self.assertIn("Error: codeml executable not found at 'codeml'.", process.stdout)


    def test_input_file_not_found(self):
        """Test script behavior when input files are not found."""
        non_existent_alignment = os.path.join(self.temp_dir_path, "no_align.fasta")
        non_existent_tree = os.path.join(self.temp_dir_path, "no_tree.nwk")

        # Test non-existent alignment
        args_no_align = [
            "--alignment", non_existent_alignment,
            "--tree", self.tree_file,
            "--model", "M0",
            "--outfile_prefix", "test_no_align"
        ]
        process_no_align = self._run_script(args_no_align)
        self.assertNotEqual(process_no_align.returncode, 0)
        self.assertIn(f"Error: Alignment file not found: {non_existent_alignment}", process_no_align.stdout)

        # Test non-existent tree
        args_no_tree = [
            "--alignment", self.alignment_file,
            "--tree", non_existent_tree,
            "--model", "M0",
            "--outfile_prefix", "test_no_tree"
        ]
        process_no_tree = self._run_script(args_no_tree)
        self.assertNotEqual(process_no_tree.returncode, 0)
        self.assertIn(f"Error: Tree file not found: {non_existent_tree}", process_no_tree.stdout)

    def test_invalid_model_name(self):
        """Test script behavior with an invalid model name."""
        args = [
            "--alignment", self.alignment_file,
            "--tree", self.tree_file,
            "--model", "M99", # Invalid model
            "--outfile_prefix", "test_invalid_model"
        ]
        process = self._run_script(args)
        self.assertNotEqual(process.returncode, 0)
        # argparse error messages go to stderr
        self.assertIn("invalid choice: 'M99'", process.stderr)

    @unittest.skipUnless(CODEML_AVAILABLE, "codeml not found in PATH, skipping codeml execution test.")
    def test_successful_run_M0(self):
        """Test a successful run with model M0."""
        prefix = "m0_success"
        outfile_prefix_path = os.path.join(self.temp_dir_path, prefix)
        args = [
            "--alignment", self.alignment_file,
            "--tree", self.tree_file,
            "--model", "M0",
            "--outfile_prefix", outfile_prefix_path,
            "--cleandata", "1" # Ensure it runs with this common option
        ]
        process = self._run_script(args)
        
        # PAML can print warnings to stderr even on success, so check stdout for key success messages
        # and primarily rely on return code and file existence.
        if process.returncode != 0:
            print(f"M0 run STDOUT:\n{process.stdout}")
            print(f"M0 run STDERR:\n{process.stderr}")
        self.assertEqual(process.returncode, 0, "Script failed for M0 model.")

        mlc_file = outfile_prefix_path + ".mlc"
        tsv_file = outfile_prefix_path + "_site_analysis.tsv"
        self.assertTrue(os.path.exists(mlc_file), f"PAML output file {mlc_file} not created.")
        self.assertTrue(os.path.exists(tsv_file), f"Summary TSV file {tsv_file} not created.")

        with open(tsv_file, "r") as f:
            lines = f.readlines()
        
        self.assertTrue(len(lines) >= 1, "TSV file is empty.")
        expected_header = "Site\tAminoAcid\tdN_dS\tPosteriorProbability_PositiveSelection\tNote\n"
        self.assertEqual(lines[0], expected_header, "TSV header is incorrect.")
        
        # For M0, site-specific dN/dS is not really applicable, the script writes a note.
        # Let's check that.
        if len(lines) > 1:
             self.assertIn("Site-specific dN/dS is not applicable for model M0", lines[1])


    @unittest.skipUnless(CODEML_AVAILABLE, "codeml not found in PATH, skipping codeml execution test.")
    def test_successful_run_M8_beb(self):
        """Test a successful run with model M8 for BEB results."""
        prefix = "m8_beb_success"
        outfile_prefix_path = os.path.join(self.temp_dir_path, prefix)
        args = [
            "--alignment", self.alignment_file,
            "--tree", self.tree_file,
            "--model", "M8",
            "--outfile_prefix", outfile_prefix_path,
            "--cleandata", "1"
        ]
        process = self._run_script(args)

        if process.returncode != 0:
            print(f"M8 run STDOUT:\n{process.stdout}")
            print(f"M8 run STDERR:\n{process.stderr}")
        self.assertEqual(process.returncode, 0, "Script failed for M8 model.")

        mlc_file = outfile_prefix_path + ".mlc"
        tsv_file = outfile_prefix_path + "_site_analysis.tsv"
        self.assertTrue(os.path.exists(mlc_file), f"PAML output file {mlc_file} not created.")
        self.assertTrue(os.path.exists(tsv_file), f"Summary TSV file {tsv_file} not created.")

        with open(tsv_file, "r") as f:
            lines = f.readlines()
        
        self.assertTrue(len(lines) >= 1, "TSV file is empty.")
        expected_header = "Site\tAminoAcid\tdN_dS\tPosteriorProbability_PositiveSelection\tNote\n"
        self.assertEqual(lines[0], expected_header, "TSV header is incorrect for M8.")
        
        num_codons = 4 # From FASTA_CONTENT
        # Expect header + num_codons lines
        self.assertEqual(len(lines), num_codons + 1, f"TSV file should have {num_codons+1} lines (header + data).")

        # Check that data lines are populated for M8
        for i in range(1, len(lines)):
            parts = lines[i].strip().split('\t')
            self.assertEqual(len(parts), 5, f"TSV data line {i} does not have 5 columns.")
            self.assertTrue(parts[0].isdigit(), f"Site column in line {i} is not a digit.")
            self.assertTrue(len(parts[1]) > 0, f"AminoAcid column in line {i} is empty.") # AA can be X if cleandata=0
            try:
                float(parts[2]) # dN_dS
                float(parts[3]) # PosteriorProbability_PositiveSelection
            except ValueError:
                self.fail(f"dN_dS or PosteriorProbability in line {i} is not a float: {lines[i]}")
            # Note column can be empty if no positive selection detected with high probability
            # self.assertTrue(len(parts[4]) > 0, f"Note column in line {i} is empty.")


if __name__ == "__main__":
    unittest.main()
