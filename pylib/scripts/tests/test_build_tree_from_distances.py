#!/usr/bin/env python3

import unittest
import subprocess
import os
import sys
from Bio import Phylo

class TestBuildTreeFromDistances(unittest.TestCase):

    def setUp(self):
        """Set up test environment."""
        self.script_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) # pylib/scripts/
        self.script_path = os.path.join(self.script_dir, "build_tree_from_distances.py")
        self.data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
        self.sample_matrix_path = os.path.join(self.data_dir, "sample_distance_matrix.phy")
        
        self.output_tree_file_nj = "test_tree_nj.nwk"
        self.output_tree_file_upgma = "test_tree_upgma.nwk"

        # Ensure no old output files are present
        if os.path.exists(self.output_tree_file_nj):
            os.remove(self.output_tree_file_nj)
        if os.path.exists(self.output_tree_file_upgma):
            os.remove(self.output_tree_file_upgma)

    def tearDown(self):
        """Clean up after tests."""
        if os.path.exists(self.output_tree_file_nj):
            os.remove(self.output_tree_file_nj)
        if os.path.exists(self.output_tree_file_upgma):
            os.remove(self.output_tree_file_upgma)

    def test_nj_tree_creation(self):
        """Test tree creation using Neighbor-Joining (NJ) method."""
        command_nj = [
            sys.executable, self.script_path,
            self.sample_matrix_path,
            "--method", "nj",
            "--outfile", self.output_tree_file_nj,
            "--informat", "nexus" # Script uses Nexus parser for PHYLIP-like matrices
        ]
        
        result_nj = subprocess.run(command_nj, capture_output=True, text=True)
        
        # Check for successful execution
        # Bio.Nexus can sometimes print PAML citation to stderr, which is not an error.
        # We check if stderr contains "Error:" or if returncode is non-zero
        self.assertEqual(result_nj.returncode, 0, f"NJ Script execution failed: {result_nj.stderr}")
        if "Error:" in result_nj.stderr: # Stricter check for actual errors
            self.fail(f"NJ Script execution produced an error: {result_nj.stderr}")

        # Verify output file was created
        self.assertTrue(os.path.exists(self.output_tree_file_nj), "NJ output tree file was not created.")
        
        # Try to parse the tree and check terminals
        try:
            tree_nj = Phylo.read(self.output_tree_file_nj, "newick")
            self.assertEqual(len(tree_nj.get_terminals()), 4, "NJ tree does not have 4 terminals.")
        except Exception as e:
            self.fail(f"Could not parse NJ output tree: {e}")

    def test_upgma_tree_creation(self):
        """Test tree creation using UPGMA method."""
        command_upgma = [
            sys.executable, self.script_path,
            self.sample_matrix_path,
            "--method", "upgma",
            "--outfile", self.output_tree_file_upgma,
            "--informat", "nexus" # Script uses Nexus parser for PHYLIP-like matrices
        ]
        
        result_upgma = subprocess.run(command_upgma, capture_output=True, text=True)
        
        self.assertEqual(result_upgma.returncode, 0, f"UPGMA Script execution failed: {result_upgma.stderr}")
        if "Error:" in result_upgma.stderr: # Stricter check for actual errors
            self.fail(f"UPGMA Script execution produced an error: {result_upgma.stderr}")

        self.assertTrue(os.path.exists(self.output_tree_file_upgma), "UPGMA output tree file was not created.")
        
        try:
            tree_upgma = Phylo.read(self.output_tree_file_upgma, "newick")
            self.assertEqual(len(tree_upgma.get_terminals()), 4, "UPGMA tree does not have 4 terminals.")
        except Exception as e:
            self.fail(f"Could not parse UPGMA output tree: {e}")

if __name__ == '__main__':
    unittest.main()
