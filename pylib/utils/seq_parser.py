from Bio import SeqIO

def parse_fasta_file(fasta_file_path):
    """Parses a FASTA file and yields SeqRecord objects."""
    with open(fasta_file_path, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            yield record

def parse_genbank_file(genbank_file_path):
    """Parses a GenBank file and yields SeqRecord objects."""
    with open(genbank_file_path, "r") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            yield record

def write_fasta_file(records, fasta_file_path):
    """Writes SeqRecord objects to a FASTA file."""
    with open(fasta_file_path, "w") as handle:
        SeqIO.write(records, handle, "fasta")

# Add more functions as needed for other formats (Phylip, etc.)
# based on the Perl scripts' capabilities.
