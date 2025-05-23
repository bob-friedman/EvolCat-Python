# Perl to Python Bioinformatics Library

This project is a Python-based bioinformatics library converted from an existing collection of Perl scripts. It provides a set of tools for common bioinformatics tasks such as sequence manipulation, format conversion, and analysis.

## Overview

The library is organized into utility modules and executable scripts.
- `pylib/utils/`: Contains core utility modules for tasks like sequence parsing.
- `pylib/scripts/`: Contains executable Python scripts that replicate the functionality of the original Perl command-line tools.

## User Responsibility

The scripts and library components provided here are for research and informational purposes. Users are responsible for validating the results obtained using this software, interpreting them correctly, and ensuring they are appropriate for their specific application. The original authors and the converters of this code disclaim any liability for its use or misuse. It is recommended to test the tools with known datasets and compare results with other established bioinformatics software where appropriate. Users may need to adapt or modify the code to suit their specific research needs and computational environment.

## Usage

See `docs/USAGE.md` for details on how to use individual scripts.

## Dependencies

- Python 3.x
- Biopython
- Matplotlib
- (Other dependencies may be added as more scripts are converted)

## Development and Contributions

This library was primarily converted from its original Perl source using AI-assisted tooling (Model: Gemini, Agent: Jules). The process involved:
*   Analyzing original Perl scripts.
*   Translating Perl logic to Python, utilizing standard libraries like Biopython and Matplotlib.
*   Structuring the code into a Python package with scripts and utility modules.
*   Creating command-line interfaces using `argparse`.
*   Setting up packaging with `setup.py`.
*   Developing initial documentation.

Human oversight and review are crucial for ensuring the accuracy and robustness of the converted code.

*This library is currently under development.*
