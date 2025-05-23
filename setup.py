from setuptools import setup, find_packages

setup(
    name="bio_perl_conversion",
    version="0.1.0",
    author="Original Perl Authors & AI Conversion",
    description="A Python library converted from Perl bioinformatics scripts.",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    packages=find_packages(where=".", include=['pylib', 'pylib.*']),
    install_requires=[
        "biopython>=1.79",  # Specify a reasonable version
        "matplotlib>=3.0"   # Specify a reasonable version
    ],
    entry_points={
        'console_scripts': [
            'gb2fasta_py=pylib.scripts.gb2fasta:main',
            'translate_seq_py=pylib.scripts.translate_seq:main',
            'join_tables_py=pylib.scripts.join_tables:main',
            'dot_plot_py=pylib.scripts.dot_plot:main',
            'fas2csv_py=pylib.scripts.fas2csv:main',
            'fas2meg_py=pylib.scripts.fas2meg:main',
            'fas2phy_py=pylib.scripts.fas2phy:main',
            'phy2fas_py=pylib.scripts.phy2fas:main',
            'phy2meg_py=pylib.scripts.phy2meg:main',
            'gbCDS_py=pylib.scripts.gbCDS:main',
            'clean_fasta_name_py=pylib.scripts.clean_fasta_name:main',
            'extract_region_py=pylib.scripts.extract_region:main',
            'merge_fastas_py=pylib.scripts.merge_fastas:main',
            'rev_comp_py=pylib.scripts.rev_comp:main',
            'nogaps_py=pylib.scripts.nogaps:main',
            'print_dup_col_vals_py=pylib.scripts.print_duplicate_column_values:main',
            'sort_numbers_py=pylib.scripts.sort_numbers:main',
            'table_to_binary_py=pylib.scripts.table_to_binary:main',
            'transpose_tsv_py=pylib.scripts.transpose_tsv:main',
            'transpose_text_py=pylib.scripts.transpose_text_matrix:main',
            'iupac_to_regexp_py=pylib.scripts.iupac_to_regexp:main',
            'count_kmers_py=pylib.scripts.count_kmers:main',
            'approx_match_py=pylib.scripts.approximate_string_match:main',
            'parse_blast_text_py=pylib.scripts.parse_blast_text:main',
            'blast_to_table_py=pylib.scripts.blast_to_table:main',
            'find_blast_top_pairs_py=pylib.scripts.find_blast_top_pairs:main',
            'calculate_k2p_py=pylib.scripts.calculate_k2p:main',
            'calculate_dna_distances_py=pylib.scripts.calculate_dna_distances:main',
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License", # Assuming MIT based on typical open source bio projects, adjust if known otherwise
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires='>=3.7', # Specify a reasonable minimum Python version
)
