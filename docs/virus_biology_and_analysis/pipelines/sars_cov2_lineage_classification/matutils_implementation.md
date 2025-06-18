# `matUtils`: Implementation Details and Advanced Usage

This document provides further details on how `matUtils`, a key component of the UShER toolkit, is utilized within the SARS-CoV-2 lineage classification pipeline. It also outlines a conceptual strategy for extracting data for all clades from the comprehensive Mutation-Annotated Tree (MAT) `.pb` file.

## Current `matUtils` Usage for Single Clade Extraction

As detailed in the main [SARS-CoV-2 Lineage Classification Pipeline README](./README.md) and the accompanying Python script (`sars_cov2_lineage_classifier.py`), `matUtils` is employed to extract variant data for specific clades of interest.

The primary command used for this purpose is:

```bash
matUtils extract -i public-latest.all.masked.pb -c "CLADE_NAME" -v output_clade_specific.vcf
```

Where:
*   `--input-mat public-latest.all.masked.pb`: Specifies the input MAT file (e.g., `public-latest.all.masked.pb`).
*   `--clade "CLADE_NAME"`: Focuses the extraction on the members of the specified clade. The `CLADE_NAME` must exactly match a name present in the MAT file's metadata (e.g., 'XBB.2.3', "20H (Beta)"). Note the use of quotes for clade names containing spaces or special characters.
*   `-v output_clade_specific.vcf`: Instructs `matUtils` to output the data for the specified clade in VCF (Variant Call Format) to the designated file.

This targeted extraction is crucial for manageable, focused analyses as demonstrated in the main pipeline, allowing for detailed feature engineering and model training on a per-clade or per-group-of-clades basis.

## Conceptual Strategy for Extracting All Clades (Divide and Conquer)

While the main pipeline focuses on specific clades, there is significant interest in analyzing data across *all* available clades within the MAT `.pb` file. Given the immense size of the full dataset (over 8 million samples), a direct extraction of all data into a single file is impractical due to prohibitive memory and time requirements.

A more feasible "divide-and-conquer" strategy can be implemented as follows:

1.  **Retrieve the Full List of Clades:**
    The UShER toolkit provides a way to list all clades annotated within the MAT. This can be achieved using the `matUtils summary` command:
    ```bash
    matUtils summary --input-mat public-latest.all.masked.pb --clades clades.tsv
    ```
    This command will generate a tab-separated values (TSV) file (e.g., `clades.tsv`) containing the names of all clades in the `.pb` file.

2.  **Iterative Extraction for Each Clade:**
    With the list of all clades, an iterative process can be scripted (e.g., using a bash loop or a Python script) to extract data for each clade individually. The core of this iteration would be:

    *   **Read Clade Names:** The script would parse the `clades.tsv` file to get each clade name.
    *   **Execute `matUtils extract`:** For each clade name, the script would dynamically construct and execute the `matUtils extract` command. It's crucial to output each clade's data to a unique file to prevent overwriting. For example:
        ```bash
        # Example in a bash script loop
        # while read CLADE_NAME; do
        #   echo "Extracting clade: ${CLADE_NAME}"
        #   matUtils extract -i public-latest.all.masked.pb -c "${CLADE_NAME}" -v "${CLADE_NAME}_variants.vcf"
        # done < clades.tsv
        ```
        Or, using Python with the `subprocess` module:
        ```python
        # import subprocess
        # import pandas as pd
        #
        # clades_df = pd.read_csv("clades.tsv", sep="\t") # Assuming header is 'clade'
        # for clade_name in clades_df["clade"]:
        #     print(f"Extracting clade: {clade_name}")
        #     output_filename = f"{clade_name.replace('/', '_').replace(' ', '_')}_variants.vcf" # Sanitize filename
        #     command = [
        #         "matUtils", "extract",
        #         "-i", "public-latest.all.masked.pb",
        #         "-c", clade_name,
        #         "-v", output_filename
        #     ]
        #     subprocess.run(command, check=True)
        ```

3.  **Managing Resources and Output:**
    *   **Time and Memory:** This iterative process will still be time-consuming, but it breaks the problem into manageable chunks, significantly reducing the peak memory required compared to a bulk extraction.
    *   **Storage:** Ensure sufficient disk space is available, as a VCF file will be generated for each clade.
    *   **Error Handling:** The script should include robust error handling to manage cases where a clade name might be problematic or an extraction fails.

4.  **Downstream Analysis:**
    Once all clades have their individual VCF files, these can be processed further. For example, scripts can be developed to parse these VCFs, consolidate relevant information (variants, sample IDs per variant/clade), and build a comprehensive table suitable for large-scale comparative analysis, population genetics studies, or identifying unique/shared mutations across the entire SARS-CoV-2 phylogeny.

This divide-and-conquer approach, while requiring careful scripting and resource management, provides a viable path to systematically access and analyze the rich dataset contained within the UShER MAT `.pb` files.
