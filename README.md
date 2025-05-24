# Perl to Python Bioinformatics Library

This project is a Python-based bioinformatics library converted from an existing collection of Perl scripts. It provides a set of tools for common bioinformatics tasks such as sequence manipulation, format conversion, and analysis.

## Overview

The library is organized into utility modules and executable scripts.
- `pylib/utils/`: Contains core utility modules for tasks like sequence parsing.
- `pylib/scripts/`: Contains executable Python scripts that replicate the functionality of the original Perl command-line tools.

## User Responsibility

The scripts and library components provided here are for research and informational purposes. Users are responsible for validating the results obtained using this software, interpreting them correctly, and ensuring they are appropriate for their specific application. The original authors and the converters of this code disclaim any liability for its use or misuse. It is recommended to test the tools with known datasets and compare results with other established bioinformatics software where appropriate. Users may need to adapt or modify the code to suit their specific research needs and computational environment.

## Usage

See `docs/USAGE.md` for details on how to use individual scripts. Generally, many of the scripts depend on the pylib core utility modules. These scripts will find the modules by default when they are located in EvolCat-Python with pylib/utils as a subdirectory.

## Access in a Windows OS with WSL

To access these scripts from a linux environment in supported versions of Windows, the first step to verify the installation of WSL. The next steps involve use of "pip" to manage packages in Python, but the default version with Ubuntu may not be compatible with it in WSL. To workaround the issue, it is possible to follow third-party procedures. These are not recommended, but I will describe steps below that worked for my system.

1. In the Ubuntu shell: sudo apt install -y gcc make build-essential libssl-dev libffi-dev zlib1g-dev libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm libncurses5-dev xz-utils tk-dev libxml2-dev libxmlsec1-dev liblzma-dev
2. curl https://pyenv.run | bash
3. Add to .bashrc file:
export PYENV_ROOT="$HOME/.pyenv"
command -v pyenv >/dev/null || export PATH="$PYENV_ROOT/bin:$PATH"
eval "$(pyenv init -)"
4. pyenv install 3.10.0
5. Choose the default version. I chose global. This version is not too far from the one installed.
pyenv global 3.10.0  # Set 3.10.0 as the default version of Python
pyenv local 3.10.0   # Set 3.10.0 as the version when running within the current folder

## Dependencies

- Python 3.x
- Biopython
- Matplotlib
- (Other dependencies may be added as more scripts are converted)

- EvolCat/Python/pylib/scripts/ncbi/ has experimental scripts to access NCBI sequence data, but also depends on a module Requests.

## Development and Contributions

This library was primarily converted from its original Perl source using AI-assisted tooling (Model: Gemini 2.5 Pro, Agent: Jules AI). The process involved:
*   Analyzing original Perl scripts.
*   Translating Perl logic to Python, utilizing standard libraries like Biopython and Matplotlib.
*   Structuring the code into a Python package with scripts and utility modules.
*   Creating command-line interfaces using `argparse`.
*   Setting up packaging with `setup.py`.
*   Developing initial documentation.

Human oversight and review are crucial for ensuring the accuracy and robustness of the converted code.

*This library is currently under development.*
