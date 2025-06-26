# Solvent Accessible Surface Area (SASA)

FreeSASA is an open-source tool for calculating the Solvent Accessible Surface Area (SASA) of biomolecules. It can be used as a command-line tool, a C library, or a Python module.

## Applications in Bioinformatics

SASA is a fundamental property of proteins and is used in various bioinformatics applications:

*   **Understanding Protein Structure and Function:** SASA helps identify surface and core residues, providing insights into protein folding, stability, and interactions.
*   **Analyzing Molecular Interactions:** Changes in SASA upon complex formation can identify interaction interfaces and estimate binding affinity.
*   **Predicting Protein Solubility:** The ratio of polar to non-polar surface area can indicate protein solubility.
*   **Drug Design:** SASA calculations can help identify potential binding pockets on target proteins.

## Algorithms

FreeSASA implements two common algorithms for SASA calculation:

*   **Lee & Richards' Algorithm (Default):** This algorithm slices the molecule and calculates accessible arc lengths for each atom in each slice.
*   **Shrake & Rupley's Algorithm:** This algorithm places points on each atom's surface and checks their accessibility to the solvent.

Users can choose the algorithm and adjust parameters to balance accuracy and speed.

## Studying Virus Evolution

SASA calculations are relevant for studying virus evolution:

*   **Analyzing Mutation Impacts:** Comparing SASA of wild-type and mutated proteins can quantify the structural impact of mutations, such as exposing new binding sites.
*   **Understanding Host-Virus Interactions:** Changes in SASA of viral surface proteins can affect binding to host cell receptors, influencing infectivity.
*   **Immune Evasion:** SASA calculations can help identify surface changes that allow viruses to evade the host immune system.

## Examples

### Command-Line Interface

1.  **Basic Calculation (Lee & Richards):**
    ```bash
    freesasa protein.pdb
    ```
2.  **Shrake & Rupley Algorithm with Custom Parameters:**
    ```bash
    freesasa --shrake-rupley -n 200 --probe-radius 1.2 --n-threads 4 protein.pdb
    ```
3.  **JSON Output with Residue-Level Detail:**
    ```bash
    freesasa --format=json --output-depth=residue protein.pdb
    ```

### Python Interface

1.  **Basic Calculation:**
    ```python
    import freesasa

    structure = freesasa.Structure("protein.pdb")
    result = freesasa.calc(structure)
    area_classes = freesasa.classifyResults(result, structure)

    print("Total : %.2f A2" % result.totalArea())
    for key in area_classes:
        print(key, ": %.2f A2" % area_classes[key])
    ```
2.  **Using Selections:**
    ```python
    import freesasa

    structure = freesasa.Structure("protein.pdb")
    result = freesasa.calc(structure)
    # Calculate SASA for Alanine residues and residues 1 to 10
    selections = freesasa.selectArea(('alanine, resn ala', 'r1_10, resi 1-10'),
                                     structure, result)
    for key in selections:
        print(key, ": %.2f A2" % selections[key])
    ```

For more details, refer to the [official FreeSASA documentation](https://freesasa.github.io/).
