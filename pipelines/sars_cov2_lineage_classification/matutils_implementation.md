# `matUtils`: A Deep Dive

This document provides advanced details on how `matUtils` is used in this pipeline and outlines a robust strategy for extracting data for *all* clades from a Mutation-Annotated Tree (MAT).

## Part 1: Single-Clade Extraction (Current Usage)

The main pipeline uses `matUtils` to perform targeted extractions of specific clades. This allows for focused analysis on a per-clade basis.

The core command for this operation is:

```bash
matUtils extract -i public-latest.all.masked.pb -c "CLADE_NAME" -v output_clade_specific.vcf
```

**Command Breakdown:**

*   `-i, --input-mat`: Specifies the input MAT protobuf file (e.g., `public-latest.all.masked.pb`).
*   `-c, --clade`: Filters the tree to include only samples belonging to `CLADE_NAME`. The name must be an exact match from the tree's metadata.
    > **Tip:** Use quotes for clade names containing spaces or special characters (e.g., `"20H (Beta, V2)"`).
*   `-v, --vcf`: Specifies the output VCF (Variant Call Format) file for the extracted clade data.

This approach is ideal for the pipeline's goal: manageable feature engineering and model training on a limited set of clades.

---

## Part 2: Strategy for Full-Scale Data Extraction

A common goal is to analyze data across all clades in the MAT. However, a direct, single extraction is often computationally infeasible.

> **⚠️ The Challenge: Scale**
> The full UShER MAT contains millions of samples. Attempting to extract all data into a single file would likely exhaust system memory and take an impractical amount of time.

The recommended solution is a **"divide-and-conquer"** strategy, which breaks the problem into manageable steps.

### The Divide-and-Conquer Workflow

#### **Step 1: Get All Clade Names**

First, generate a complete list of all clades present in the MAT file using the `matUtils summary` command.

```bash
matUtils summary -i public-latest.all.masked.pb --clades clades.tsv
```

This creates a clean, tab-separated file (`clades.tsv`) containing one clade name per line, ready for iteration.

#### **Step 2: Iteratively Extract Each Clade**

Next, use a script to loop through `clades.tsv` and run the `matUtils extract` command for each clade, saving each one to a unique file.

**Example: Bash Loop**
```bash
#!/bin/bash
INPUT_MAT="public-latest.all.masked.pb"
CLADE_LIST="clades.tsv"

while read -r CLADE_NAME; do
  # Sanitize clade name for use as a filename
  SANITIZED_NAME=$(echo "${CLADE_NAME}" | tr '/' '_')
  echo "Extracting clade: ${CLADE_NAME}"
  
  matUtils extract -i "${INPUT_MAT}" -c "${CLADE_NAME}" -v "${SANITIZED_NAME}_variants.vcf"
done < "${CLADE_LIST}"
```

**Example: Python Script**
```python
import subprocess
import pandas as pd

INPUT_MAT = "public-latest.all.masked.pb"
CLADE_LIST = "clades.tsv"

clades_df = pd.read_csv(CLADE_LIST, sep="\t", header=None, names=["clade"])

for clade_name in clades_df["clade"]:
    print(f"Extracting clade: {clade_name}")
    
    # Sanitize for safe filenames (e.g., "B.1/Q.1" -> "B.1_Q.1")
    output_filename = f"{clade_name.replace('/', '_').replace(' ', '_')}_variants.vcf"
    
    command = [
        "matUtils", "extract",
        "-i", INPUT_MAT,
        "-c", clade_name,
        "-v", output_filename
    ]
    subprocess.run(command, check=True)
```

#### **Step 3: Manage the Large-Scale Job**

This iterative process is resource-intensive. Plan accordingly:

*   **Time & CPU:** The full extraction will take many hours or even days, depending on your hardware.
*   **Storage:** Ensure you have sufficient disk space. The total size of the individual VCF files can be substantial.
*   **Error Handling:** A robust script should include error handling to log failed extractions and continue the process, rather than halting on the first error.

#### **Step 4: Consolidate and Analyze**

Once all individual VCF files are generated, they become the foundation for comprehensive downstream analysis. You can now develop scripts to parse these files, consolidate key information, and build a master dataset for comparative genomics, identifying unique mutations, or tracking variant prevalence across the entire phylogeny.

### Summary

This divide-and-conquer workflow provides a scalable and robust path to systematically process the entire UShER MAT, transforming an intractably large dataset into a collection of manageable files ready for deep analysis.
