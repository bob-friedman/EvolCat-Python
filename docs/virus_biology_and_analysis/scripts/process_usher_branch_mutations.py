# 1. Mount Google MyDrive in Google Colab
from google.colab import drive
import os

# Mount Google MyDrive for use in Google Colab
# In Colab: had to activate all allowable permissions for access to MyDrive
#  to avoid authentication error
drive.mount('/content/drive')
os.chdir('/content/drive/MyDrive')

# 2. Installation of Conda in Google Colab
!pip install -q condacolab
import condacolab; condacolab.install()
# Initialize the shell and restart the kernel
# Colab typically and expectedly restarts with a log/report on `crash reported`
!conda init bash
# Verify installation
!conda --version

# 3. Installation of the UShER toolkit via Conda
# Otherwise, see documentation for other options:
# https://usher-wiki.readthedocs.io/en/latest/Installation.html
# Create a new environment for UShER
!conda create -n usher-env # python=3.10 to support BTE library, if installed
# Activate the newly created environment
!conda activate usher-env

# Set up channels
!conda config --add channels defaults
!conda config --add channels bioconda
!conda config --add channels conda-forge
# Install the UShER package
!conda install -q usher

# 4. Download the latest UShER Mutation-Annotated Tree (MAT) data
!wget http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.all.masked.pb.gz

# Uncompress the MAT data file (.pb extension name; -f parameter will force a file overwrite)
# Creates public-latest.all.masked.pb
!gunzip -f public-latest.all.masked.pb.gz

# Export summary data associated with the MAT file (e.g., --clades, --node-stats,
# --mutations, --samples, --get-all)
!matUtils summary --input-mat public-latest.all.masked.pb --clades clades.tsv
# !matUtils summary --input-mat public-latest.all.masked.pb --samples samples.txt

# 5. Obtain mutation data for each node in the subtree

# If any issues arise, verify that public-latest.all.masked.pb is in the current working directory
# Replace "YOUR_CLADE_OF_INTEREST" with the actual clade name, e.g., "20H (Beta)"
# May replace "mutations_for_clade.txt" with another output filename

# Tested with clade `20H (Beta)` clade of SARS-CoV-2 with 10179 samples:
# If scaling to a larger file size, note the full SARS-CoV-2 dataset is ~800x as large
# as clade `20H (Beta)`

!matUtils extract \
    --input-mat public-latest.all.masked.pb \
    --clade "YOUR_CLADE_OF_INTEREST" \
    --all-paths mutations_for_clade.txt

# Explanation of the command:
# `--input-mat public-latest.all.masked.pb`: Specifies the input MAT file.
# `--clade "YOUR_CLADE_OF_INTEREST"`: Focuses the extraction on the members of the named
#  clade. This name must exactly match a clade name present in the MAT file's metadata.
#  May specify multiple clade names as a comma-delimited list. Add double quotes to the
#  names with spaces.
# `--all-paths mutations_for_clade.txt`: This crucial option tells `matUtils` to output the mutations
#  along each path from the clade's common ancestor to every sample and internal node within
#  that clade. The output is saved to ` mutations_for_clade.txt`. The list is created by a depth-first
#  traversal order.

# Output Format:
# The output file (`mutations_for_clade.txt`) will typically list each node (internal nodes often
# labeled like `node_X:` or sample (e.g., `Country/SampleID/Date|Accession|Date:`) followed by
# the mutations inferred to have occurred on the branch immediately leading to it. For example:
#  node_1: G15910T
#  Sample/ID/Date|Accession|Date: C1191T,C11674T
#  node_2: T13090C

# This detailed mutation information is invaluable for understanding the specific evolutionary
# changes within a lineage and can serve as input for further analyses, including preparing data for # training predictive models like Transformers.

"""
# (Optional) Convert VCF formatted file to Fasta formatted sequence data
# The vcf2fasta binary fails to run in Colab or WSL of Windows 11. Untested in other environments.

# Installation of VCF library
!conda install -q vcflib
!echo "Current working directory: $PWD"

# Download reference sequence for reconstruction of SARS-CoV-2 Fasta sequences from vcf file
!wget https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/chromosomes/NC_045512v2.fa.gz
!gunzip -f NC_045512v2.fa.gz

# The vcf file lists each nucleotide variant of genotypes as compared
# to a reference Fasta formatted sequence (specified by --reference).
!vcfindex my_clade.vcf > my_clade_idx.vcf
!vcf2fasta --reference NC_045512v2.fa my_clade_idx.vcf
"""
