import unittest
import subprocess
import os
import tempfile

# Path to the script to be tested
SCRIPT_PATH = os.path.join(os.path.dirname(__file__), "..", "calculate_nucleotide_diversity.py")

class TestCalculateNucleotideDiversity(unittest.TestCase):

    def run_script(self, fasta_content=None, file_path=None):
        """Helper method to run the script and return its output."""
        if fasta_content is not None:
            with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fasta") as tmp_fasta:
                tmp_fasta.write(fasta_content)
                script_args = ["python3", SCRIPT_PATH, tmp_fasta.name]
                temp_file_name = tmp_fasta.name
        elif file_path is not None:
            script_args = ["python3", SCRIPT_PATH, file_path]
            temp_file_name = None # No temporary file created by this helper in this case
        else:
            raise ValueError("Either fasta_content or file_path must be provided")

        process = subprocess.run(script_args, capture_output=True, text=True)
        
        if temp_file_name and fasta_content is not None: # Only delete if this helper created it
            os.remove(temp_file_name)
            
        return process

    def test_valid_alignment(self):
        """Test with a valid FASTA alignment."""
        fasta_content = """>Seq1
ATGC
>Seq2
ATGT
>Seq3
AGGC
"""
        # Expected pi calculation:
        # Seq1: ATGC
        # Seq2: ATGT
        # Seq3: AGGC
        # Alignment length = 4
        # Pairs:
        # (Seq1, Seq2): diff = 1 (G vs T at pos 3)
        # (Seq1, Seq3): diff = 1 (T vs G at pos 1)
        # (Seq2, Seq3): diff = 2 (T vs G at pos 1, T vs C at pos 3)
        # Sum of diffs = 1 + 1 + 2 = 4
        # Number of sequences (n) = 3
        # Number of pairs = 3 * (3-1) / 2 = 3
        # pi = (Sum of diffs) / (Number of pairs) / (Alignment length)
        # pi = 4 / 3 / 4 = 1/3 = 0.333333...
        expected_pi_value = 1/3
        
        process = self.run_script(fasta_content=fasta_content)
        
        self.assertEqual(process.returncode, 0, f"Script failed with error: {process.stderr}")
        self.assertIn("Nucleotide diversity (pi):", process.stdout)
        
        # Extract pi value from output
        output_pi_str = None
        for line in process.stdout.splitlines():
            if "Nucleotide diversity (pi):" in line:
                output_pi_str = line.split(":")[1].strip()
                break
        
        self.assertIsNotNone(output_pi_str, "Pi value not found in script output.")
        try:
            output_pi_value = float(output_pi_str)
        except ValueError:
            self.fail(f"Could not convert pi value '{output_pi_str}' to float.")
            
        self.assertAlmostEqual(output_pi_value, expected_pi_value, places=5, 
                               msg=f"Expected pi {expected_pi_value}, got {output_pi_value}")

    def test_unaligned_sequences(self):
        """Test with unaligned sequences (different lengths)."""
        fasta_content = """>Seq1
ATGC
>Seq2
ATGCA
>Seq3
AGGC
"""
        process = self.run_script(fasta_content=fasta_content)
        self.assertNotEqual(process.returncode, 0, "Script should exit with an error for unaligned sequences.")
        self.assertIn("Error: Sequences in the FASTA file are not all of the same length.", process.stdout, # Script prints to stdout for now
                      "Error message for unaligned sequences not found.")

    def test_empty_file(self):
        """Test with an empty FASTA file."""
        fasta_content = ""
        process = self.run_script(fasta_content=fasta_content)
        # The script might exit with 0 if no sequences are found but prints a message.
        # Or it might exit with non-zero. Let's check for the message.
        self.assertIn("No sequences found in the FASTA file.", process.stdout,
                      "Error message for empty file not found.")

    def test_single_sequence(self):
        """Test with a FASTA file containing only one sequence."""
        fasta_content = """>Seq1
ATGCGTAGCAT
"""
        process = self.run_script(fasta_content=fasta_content)
        self.assertIn("Nucleotide diversity calculation requires at least two sequences.", process.stdout,
                      "Error message for single sequence not found.")

    def test_file_not_found(self):
        """Test with a non-existent file path."""
        non_existent_file = "non_existent_file.fasta"
        # Ensure the file does not exist before running the test
        if os.path.exists(non_existent_file):
            os.remove(non_existent_file)
            
        process = self.run_script(file_path=non_existent_file)
        self.assertNotEqual(process.returncode, 0, "Script should exit with an error for file not found.")
        self.assertIn(f"Error: File not found at {non_existent_file}", process.stdout, # Script prints to stdout for now
                      "Error message for file not found not found.")

if __name__ == "__main__":
    unittest.main()
