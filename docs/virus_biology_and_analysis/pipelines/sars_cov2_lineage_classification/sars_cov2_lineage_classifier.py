# Processing SARS-CoV-2 Mutation Data

# 1. Mount Google MyDrive in Google Colab
from google.colab import drive
import os

# Mount Google MyDrive for use in Google Colab
# For Google Colab: had to activate all allowable permissions for access to MyDrive to avoid authentication error
drive.mount('/content/drive')
os.chdir('/content/drive/MyDrive')

# 2. Installation of Conda for Google Colab

# Installation of Conda for Google Colab
!pip install -q condacolab > /dev/null # suppress standard output
# Install Conda
import condacolab; condacolab.install()
# Initialize the shell and restart the kernel
# Google Colab usually restarts by itself and reports a warning message
#  however, the installation procedure is expected to continue
!conda init bash
# Verify installation
!conda --version

# 3. Installation of the UShER toolkit via Conda

# Installation of the UShER toolkit for Conda
# Otherwise, see documentation for other options:
#  https://usher-wiki.readthedocs.io/en/latest/Installation.html
# Create a new environment for UShER
!conda create -n usher-env # python=3.10 # to support BTE library, if installed
# Activate the newly created environment
!conda activate usher-env

# Set up channels
!conda config --add channels defaults
!conda config --add channels bioconda
!conda config --add channels conda-forge
# Install the UShER package
!conda install -q usher

# 4. Download the UShER Mutation-Annotated Tree (MAT) data

# Download the latest UShER Mutation-Annotated Tree (MAT) data (.pb extension name)
file_path = "/content/drive/MyDrive/public-latest.all.masked.pb"
if os.path.isfile(file_path):
    print(f"The file '{file_path}' exists.")
else:
    print(f"The file '{file_path}' does not exist.")
    !wget http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.all.masked.pb.gz

# Uncompress the MAT data file (-f parameter will force a file overwrite)
# Creates public-latest.all.masked.pb
!gunzip -f public-latest.all.masked.pb.gz

# Optional: Export a Newick tree
# !matUtils extract --input-mat public-latest.all.masked.pb --write-tree usher_tree.nwk

# Export summary data associated with the MAT file (e.g., --clades, --node-stats, --mutations, --samples, --get-all)
# !matUtils summary --input-mat public-latest.all.masked.pb --clades clades.tsv
# !matUtils summary --input-mat public-latest.all.masked.pb --samples samples.txt

# 5. Obtain mutation data for subtree

# Obtain mutation data for each node in the subtree

# Verify public-latest.all.masked.pb is in the current working directory
# Replace "YOUR_CLADE_OF_INTEREST" with the actual clade name, e.g., "20H (Beta)"
# May replace "mutations_for_clade.txt" with another output filename

# Tested with clade `20H (Beta)` clade of SARS-CoV-2 with 10179 samples:
# If scaling to larger file size, note the full SARS-CoV-2 dataset is
# ~800x as many samples as clade `20H (Beta)`

# Command below to extract mutations for clade "20H (Beta)"
# leads to file output size of ~500 Mb
# !matUtils extract --input-mat public-latest.all.masked.pb --clade "20H (Beta)" --all-paths mutations_for_clade.txt

# Clade         Samples
# 20H (Beta)    10179
# 20J (Gamma)   32210
# 23E (XBB.2.3) 16362

# !cat public-latest.metadata.tsv | grep -nr "21H" | wc -l
# 6362
# !cat public-latest.metadata.tsv | grep -nr "20J" | wc -l
# 32317
# !cat public-latest.metadata.tsv | grep -nr "23E" | wc -l
# 21550
# !cat public-latest.metadata.tsv | grep -nr "XBB.2.3" | wc -l
# 16393

# !matUtils extract -i public-latest.all.masked.pb -c "20H (Beta)" -v my_clade.vcf
# !matUtils extract -i public-latest.all.masked.pb -c "20J (Gamma)" -v my_clade.vcf
# Try single quotes to access clade names with the '.' characters
# Try of double quotes for the clade name fails to properly identify the clade
!matUtils extract -i public-latest.all.masked.pb -c 'XBB.2.3' -v my_clade.vcf

#!matUtils extract #    --input-mat public-latest.all.masked.pb #    --clade "YOUR_CLADE_OF_INTEREST" #    --all-paths mutations_for_clade.txt

# Explanation of the command:
# `--input-mat public-latest.all.masked.pb`: Specifies the input MAT file.
# `--clade "YOUR_CLADE_OF_INTEREST"`: Focuses the extraction on the members of the named clade. This name must exactly
#  match a clade name present in the MAT file's metadata. May specify multiple clade names as a comma-delimited list.
#  Add double quotes to the names with spaces.
# `--all-paths mutations_for_clade.txt`: This crucial option tells `matUtils` to output the mutations along each path
#  from the clade's common ancestor to every sample and internal node within that clade. The output is saved to
#` mutations_for_clade.txt`. The list is created by a depth-first traversal order.

# Output Format:
# The output file (`mutations_for_clade.txt`) will typically list each node (internal nodes often labeled like
#  `node_X:`) or sample (e.g., `Country/SampleID/Date|Accession|Date:`) followed by the mutations inferred to have
#  occurred on the branch immediately leading to it. For example:
#  node_1: G15910T
#  Sample/ID/Date|Accession|Date: C1191T,C11674T
#  node_2: T13090C

# This detailed mutation information is invaluable for understanding the specific evolutionary changes within a
# lineage and can serve as input for further analyses, including preparing data for training predictive models like
# Transformers.

# 6. Create Fasta sequence from a VCF file in Colab

# Install bcftools
!conda install bcftools

# Download reference sequence for SARS-CoV-2 vcf data
file_path = "/content/drive/MyDrive/NC_045512v2.fa"
if os.path.isfile(file_path):
    print(f"The file '{file_path}' exists.")
else:
    print(f"The file '{file_path}' does not exist.")
    !wget https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/chromosomes/NC_045512v2.fa.gz

!gunzip -f NC_045512v2.fa.gz

# Compress vcf file by use of bgzip
!bgzip --force my_clade.vcf

# Index compressed file
!echo "Current working directory: $PWD"
!bcftools index --tbi my_clade.vcf.gz

# Construct consensus sequence from the VCF file and a reference sequence
!bcftools consensus -f NC_045512v2.fa my_clade.vcf.gz > consensus.fa

with open('consensus.fa', 'r') as file:
    lines=file.readline()
    print(lines[0], "
", lines[1])

# 7. Iterate over VCF file samples

# Install PySam
!conda install -q -y pysam

import pysam
from pysam import VariantFile

# The VCF file as created and indexed earlier
vcf_filepath = 'my_clade.vcf.gz'
vcf = VariantFile(vcf_filepath)

print(f"Successfully opened {vcf_filepath}")
print("-" * 30)
print(f"Reference sequence in VCF: {vcf.header.contigs}")
print(f"Samples in VCF: {list(vcf.header.samples)}")
print("-" * 30)

# Iterate over each variant record in the VCF file
# .fetch() is the standard way to iterate over an indexed VCF/BCF

variant_count = 0
log_interval = 500  # Only log detailed info every N variants

for record in vcf.fetch():
    variant_count += 1

    # Now, let's see which samples have this variant
    # The 'GT' field in the FORMAT column tells us the genotype (0=ref, 1=alt)
    samples_with_variant = []
    for sample in record.samples.values():
        # sample['GT'] returns a tuple, e.g., (1, 1) for a homozygous alt
        # We check if the first alternate allele (1) is present in the genotype
        if 1 in sample['GT']:
            samples_with_variant.append(sample.name)

    # Only log detailed info periodically
    if variant_count % log_interval == 0:
        chrom = record.chrom
        pos = record.pos
        ref_allele = record.ref
        alt_alleles = record.alts

        print(f"Processed {variant_count} variants...")
        print(f"
Variant at {chrom}:{pos} | REF: {ref_allele} | ALT: {alt_alleles[0] if alt_alleles else ''}")

        if samples_with_variant:
            print(f"  > Found in {len(samples_with_variant)} samples: {samples_with_variant[:3]}...") # Print first 3
        else:
            print("  > Not found in any sample in this VCF (might be an ancestral variant).")

# Final summary
print(f"
Processing complete. Analyzed {variant_count} total variants.")

# Close the file
vcf.close()

# 8. Feature Engineering with Pandas

# Convert the raw data from the variants output and our metadata file into a clean,
# numerical matrix for a machine learning model. This matrix `X` will have:
#   Rows: One for each sample.
#   Columns: One for each unique genetic mutation found in the dataset.
#   Values: `1` if a sample has a mutation, and `0` if it does not.

# Mount Google MyDrive in Google Colab
from google.colab import drive
import os

# Mount Google MyDrive
# For Google Colab: had to activate all allowable permissions for access to MyDrive to avoid authentication error
drive.mount('/content/drive')
os.chdir('/content/drive/MyDrive')

# Installation of PySam (recommended method in Colab)
!pip install -q pysam

import pandas as pd
import pysam
from collections import defaultdict
import gc # Import garbage collector

# --- Configuration ---
# The file with sample IDs and lineages (~1 GB file size)
METADATA_FILE = 'public-latest.metadata.tsv'
# METADATA_FILE = 'filtered_clade_metadata.csv' # load grep-filtered file if original file too large
VCF_FILE = 'my_clade.vcf.gz'

# --- Step 1: Load Metadata (Potentially Filtered) and Get Samples from VCF ---
print(f"[1/7] Loading metadata from '{METADATA_FILE}'...")
# The DtypeWarning in the screen output is common with large, mixed-type files.
# Pandas warns that it had to guess the data type for some columns, which can use more memory.
# It is not an error. For this script, it does not lead to any problems.
metadata_df = pd.read_csv(METADATA_FILE, sep='	') # Or ',' if grep output is CSV
# Verify all rows have a sample identifier.
metadata_df.dropna(subset=['strain'], inplace=True) # remove rows with NaN (missing) values for 'strain' column

print(f"
[1b/7] Extracting sample list from VCF '{VCF_FILE}'...")
vcf_in_for_samples = pysam.VariantFile(VCF_FILE)
samples_in_vcf = list(vcf_in_for_samples.header.samples)
vcf_in_for_samples.close() # Close the file
print(f"Found {len(samples_in_vcf)} samples in the VCF file. These will be the rows for matrix X.")

# The following step is critical for performance and memory management.
# By filtering the massive metadata DataFrame to *only* the samples we care about (those in our VCF),
# we create a much smaller, more manageable DataFrame `clade_metadata_df`.
clade_metadata_df = metadata_df[metadata_df['strain'].isin(samples_in_vcf)].copy()

## Using .copy() is important here. It tells pandas to create a new, independent
## DataFrame, which prevents a "SettingWithCopyWarning" later on.

# It is good practice to free RAM by deleting unused large DataFrames
del metadata_df
gc.collect()

# Detect for lack of useful data samples in the VCF file
if clade_metadata_df.empty and samples_in_vcf:
    print("WARNING: No samples from the VCF were found in the provided metadata!")
elif not samples_in_vcf:
    print("ERROR: No samples found in the VCF header. Cannot proceed.")
    exit()

# --- Step 2: Identify All Unique Variants in the VCF ---
print(f"
[2/7] Identifying all unique variants in '{VCF_FILE}'...")
vcf_in = pysam.VariantFile(VCF_FILE)
all_variants = []
# fetch is preferred for higher efficiency
for record in vcf_in.fetch():
    # This creates a new unique string identifier for each variant.
    # It combines the variant's position (record.pos), reference base (record.ref),
    # and its first alternate allele (record.alts[0]) into a format like "490_C>T".
    variant_id = f"{record.pos}_{record.ref}>{record.alts[0]}"
    all_variants.append(variant_id)

# This is a standard Python idiom to create a sorted list of unique items.
# set(all_variants) creates a collection of only the unique variant IDs.
# list(...) converts it back to a list.
# sorted(...) sorts the list, ensuring the columns in our final matrix have a consistent order.
all_variants = sorted(list(set(all_variants)))
vcf_in.close() # Close the file
print(f"Found {len(all_variants)} unique variants to use as features.")

# --- Step 3: Create a Map of Variant -> Samples ---
print("
[3/7] Mapping variants to the samples that contain them...")
# defaultdict is a specialized dictionary. If you try to access a key
# that does not exist, it automatically creates it with a default value (here, an empty list `[]`),
# preventing errors and simplifying the code in the loop below.
variant_to_samples_map = defaultdict(list)
vcf_in = pysam.VariantFile(VCF_FILE)
for record in vcf_in.fetch():
    variant_id = f"{record.pos}_{record.ref}>{record.alts[0]}"
    # Renamed 'sample' to 'sample_obj' to avoid naming conflict
    # 'sample' is a common variable name, and `pysam` itself
    # is sometimes imported as `sample`. Renaming it to `sample_obj` is safer and clearer.
    for sample_obj in record.samples.values():
        # Alternate Allele data
        # `sample_obj['GT']` retrieves the Genotype for this sample.
        # In a VCF, 0 = reference allele, 1 = first alternate allele, 2 = second, etc.
        # `if 1 in sample_obj['GT']:` checks if this sample's genotype contains the alternate allele,
        # meaning the sample *has* this mutation.
        if 1 in sample_obj['GT']:
            variant_to_samples_map[variant_id].append(sample_obj.name) # Build database of variants/samples
vcf_in.close() # Close the file
print("Finished mapping.")

print(f"Size of variant_to_samples_map: {len(variant_to_samples_map)}")
if variant_to_samples_map:
    print(f"Example entry in map (first one): {list(variant_to_samples_map.items())[0]}")

# --- Step 4: Build the Feature Matrix (X) with Pandas ---
# Create an empty DataFrame with VCF samples as rows and variants as columns.
print("
[4/7] Initializing the feature matrix X with zeros...")
# Use samples_in_vcf (derived from your clade VCF) for the index!
# Using `dtype='int8'` is a smart memory optimization. Since our values
# will only be 0 or 1, we can use the smallest integer type (`int8`), which uses
# 8 times less memory than the default (`int64`).
X = pd.DataFrame(0, index=samples_in_vcf, columns=all_variants, dtype='int8')
print(f"Matrix created with shape: {X.shape}")

# --- Step 5: Populate the Feature Matrix ---
print("
[5/7] Populating the feature matrix with variant data...")
for variant, samples_with_variant in variant_to_samples_map.items():
    # This is a "list comprehension", a concise way to build a list.
    # It iterates through all sample names in `samples_with_variant` and keeps only
    # those (`s`) that are also present in the matrix's index (`X.index`).
    # This is a safety check; in this workflow, all samples should be valid.
    valid_samples = [s for s in samples_with_variant if s in X.index]
    if valid_samples:
        # This is the most efficient way to update the matrix. Instead of looping
        # through each sample one by one, we give pandas' `.loc` accessor the entire list
        # of `valid_samples` and the `variant` column name. It then sets all corresponding cells to 1 at once.
        X.loc[valid_samples, variant] = 1
print("Matrix population complete.")

print(f"Number of samples with at least one variant: {(X.sum(axis=1) > 0).sum()}")
print(f"Number of variants present in at least one sample: {(X.sum(axis=0) > 0).sum()}")

# --- Step 6: Prepare and Filter the Labels (y) and Features (X) ---
print("
[6/7] Preparing and filtering the label vector y and feature matrix X...")
# Use the filtered clade_metadata_df from Step 1
labels_df = clade_metadata_df.set_index('strain')

# This is code for identifying singletons.
print("Identifying and removing singleton classes for stratification...")
class_counts = labels_df['pango_lineage_usher'].value_counts()
classes_to_keep = class_counts[class_counts >= 2].index
df_filtered = labels_df[labels_df['pango_lineage_usher'].isin(classes_to_keep)]

print(f"Original number of samples: {len(labels_df)}")
print(f"Number of samples after filtering singletons: {len(df_filtered)}")
print(f"Number of lineages dropped: {len(class_counts) - len(classes_to_keep)}")

# --- This is an important correction step ---
# Filter the feature matrix X to match the filtered labels.
# We use .loc to select rows from X based on the index (sample names) of our filtered dataframe.
X_filtered = X.loc[df_filtered.index]

# Now, create the final label vector from the filtered dataframe.
y_filtered_labels = df_filtered['pango_lineage_usher']

# Convert the filtered string labels into numerical labels for the model.
# This is the 'y' we will actually use for training.
y_numerical, lineage_map = pd.factorize(y_filtered_labels)
print(f"Created {len(lineage_map)} numerical labels for the lineages.")
print("Example mapping:", list(enumerate(lineage_map[:5])))

# --- Step 7: Final Verification ---
print("
[7/7] --- Feature Engineering Complete! ---")
print(f"Shape of FINAL feature matrix X: {X_filtered.shape}")
print(f"Shape of FINAL label vector y: {y_numerical.shape}")

# This check should now pass perfectly!
if X_filtered.shape[0] != len(y_numerical):
    print("ERROR: Mismatch between number of samples in X and y.")
else:
    print("SUCCESS: Feature matrix X and label vector y are perfectly aligned.")

print("
Preview of the final feature matrix (X_filtered):")
print(X_filtered.head())

print("
Preview of the final numerical labels (y_numerical):")
print(y_numerical[:10])

# Free RAM
# del X, y_numerical, labels_df, clade_metadata_df, variant_to_samples_map, all_variants, samples_in_vcf
del labels_df, clade_metadata_df, variant_to_samples_map, samples_in_vcf
gc.collect()
print(gc.collect()) # Calling it twice can sometimes be more effective.
# The `0` you see in the output means 0 unreachable objects were collected in the above function.

# --- 9. Section 9: Building and Training a Neural Network ---

# This block assumes the previous script (Feature Engineering) has been run,
# and the variables `X` (pandas DataFrame) and `y_numerical` (numpy array)
# are available in memory.

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
from sklearn.model_selection import train_test_split
from sklearn.feature_extraction.text import CountVectorizer
import numpy as np

# print("
--- Vectorizer ---")

# Create a collection of the mutation names
document = set(all_variants)

# Create a Vectorizer Object
vectorizer = CountVectorizer()
vectorizer.fit(document)

# Printing the identified unique words along with their indices
# Disable this statement unless debugging code.
# print("Vocabulary: ", vectorizer.vocabulary_)

print("
--- Starting Keras Model Building ---")

# --- Step 1: Prepare Data for Keras ---

# The feature matrix X is already in a good format (a pandas DataFrame of 0s and 1s).
# The label vector y_numerical is also ready (a numpy array of integers).

# Convert the filtered pandas DataFrame to a numpy array.
X_np = X_filtered.to_numpy()
# y_numerical is already a numpy array, so we can just rename it for clarity.
y_np = y_numerical

# --- Step 2: Split Data into Training and Testing Sets ---
# This was the source of the error. We now use the correctly aligned X_np and y_np.
# 'stratify' needs the numerical labels (y_np) to work correctly.

print(f"
Splitting data with shapes X: {X_np.shape} and y: {y_np.shape}")

# test_size=0.2 means we'll use 20% of the data for testing and 80% for training.
# stratify=y_np ensures that the proportion of each lineage is the same in both the
# training and testing sets, which is crucial for imbalanced datasets.
# random_state=42 makes the split reproducible.
X_train, X_test, y_train, y_test = train_test_split(
    X_np, y_np, # Use the aligned numpy arrays
    test_size=0.2,
    random_state=42,
    stratify=y_np # Stratify using the numerical labels
)

print("
Successfully split the filtered data!")
print(f"Shape of X_train: {X_train.shape}")
print(f"Shape of y_train: {y_train.shape}")
print(f"Shape of X_test: {X_test.shape}")
print(f"Shape of y_test: {y_test.shape}")

# --- Step 3: Define the Neural Network Architecture ---
# We will use the Keras Sequential API, which is a simple stack of layers.

# Get the number of features (from the columns of X) and the number of classes (lineages).
num_features = X_train.shape[1]
num_classes = len(np.unique(y_np)) # This should match len(lineage_map)

model = Sequential([
    # Input Layer: Dense layer with 128 neurons.
    # `input_dim` tells the model the shape of the input data (number of features).
    # `relu` (Rectified Linear Unit) is a standard, effective activation function for hidden layers.
    # Note: 'input_layer' is already the name of a layer in this model.
    Dense(128, activation='relu', input_dim=num_features, name='input_layer_'),

    # Hidden Layer: Another dense layer to learn more complex patterns.
    Dense(64, activation='relu', name='hidden_layer_1'),

    # Output Layer: The final layer.
    # The number of neurons must equal `num_classes`.
    # `softmax` activation is used for multi-class classification. It converts the
    # layer's outputs into a probability distribution, where the sum of all outputs is 1.
    # The neuron with the highest probability is the model's prediction.
    Dense(num_classes, activation='softmax', name='output_layer')
])

# Print a summary of the model's architecture
print("
Model Architecture:")
model.summary()

# --- Step 4: Compile the Model ---
# This step configures the model for training.

model.compile(
    # Optimizer: 'adam' is a popular and effective optimization algorithm. It adjusts the
    # model's internal parameters (weights) to minimize the loss.
    optimizer='adam',

    # Loss Function: 'sparse_categorical_crossentropy' is the standard choice for
    # multi-class classification when the labels (`y`) are provided as integers (0, 1, 2...).
    # If our labels were one-hot encoded, we would use 'categorical_crossentropy'.
    loss='sparse_categorical_crossentropy',

    # Metrics: 'accuracy' is what we want to monitor. It tells us the fraction of
    # samples that are correctly classified.
    metrics=['accuracy']
)

# --- Step 5: Train the Model ---
print("
--- Training the Model ---")
# `fit` is the command to start the training process.
history = model.fit(
    X_train, y_train,
    validation_data=(X_test, y_test), # Provide test data to monitor performance on unseen data during training.
    epochs=10,                        # An epoch is one complete pass through the entire training dataset.
    batch_size=32,                    # The model will update its weights after processing 32 samples at a time.
    verbose=1                         # Set to 1 to show a progress bar for each epoch.
)

# --- Step 6: Evaluate the Model ---
print("
--- Evaluating Model Performance ---")
# Use the `.evaluate` method on the test set to get the final loss and accuracy.
# This gives us the definitive measure of how well our model performs on data it has never seen before.
loss, accuracy = model.evaluate(X_test, y_test, verbose=0)
print(f"
Final Test Accuracy: {accuracy * 100:.2f}%")
print(f"Final Test Loss: {loss:.4f}")

### Code to Generate the Confusion Matrix

# --- 10. Visualizing Performance with a Confusion Matrix ---

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import confusion_matrix, classification_report
import pandas as pd # Make sure pandas is imported

print("
--- Generating Confusion Matrix and Classification Report ---")

# --- Step 1: Get Model's Predictions for the Test Set ---
y_pred_probs = model.predict(X_test)
y_pred_numerical = np.argmax(y_pred_probs, axis=1)

# --- Step 2: Identify Top N Most Frequent Classes in the Test Set ---
N_TOP_CLASSES = 15 # You can change this number
class_counts = pd.Series(y_test).value_counts()
top_n_classes_indices = class_counts.head(N_TOP_CLASSES).index.tolist()

# Get the corresponding string names for these top classes
top_n_class_names = [lineage_map[i] for i in top_n_classes_indices]

# --- Step 3: Filter the Confusion Matrix for Only Top N Classes ---
# Calculate the full confusion matrix first
cm_full = confusion_matrix(y_test, y_pred_numerical)

# Create a filtered version
cm_filtered = cm_full[top_n_classes_indices, :][:, top_n_classes_indices]

# --- Step 4: Visualize the Filtered Confusion Matrix ---
plt.figure(figsize=(12, 10)) # A good size for 15x15
sns.heatmap(cm_filtered, annot=True, fmt='d', cmap='Blues',
            xticklabels=top_n_class_names, yticklabels=top_n_class_names)
plt.title(f'Confusion Matrix (Top {N_TOP_CLASSES} Most Frequent Lineages in Test Set)', fontsize=16)
plt.ylabel('True Lineage', fontsize=12)
plt.xlabel('Predicted Lineage', fontsize=12)
plt.xticks(rotation=45, ha="right")
plt.tight_layout()

# --- SAVE THE FIGURE ---
plt.savefig(f'confusion_matrix_top_{N_TOP_CLASSES}.png', dpi=150) # Lower DPI for smaller file
plt.show()

# --- Step 5: Print a Detailed Classification Report (This can remain the same) ---
# It's okay for the report to be long, as it's text.
print("
Full Classification Report:")
present_labels = np.unique(np.concatenate((y_test, y_pred_numerical)))
filtered_target_names = [lineage_map[i] for i in present_labels]
print(classification_report(
    y_test,
    y_pred_numerical,
    labels=present_labels,
    target_names=filtered_target_names,
    zero_division=0
))

### Code to Generate Saliency Maps for Model Interpretability

# --- 11. Verifying the Model's "Reasoning" with a Saliency Table ---

# This block assumes the following are in memory from previous steps:
# model, X_test, y_test, lineage_map, and the 'vectorizer' from data prep.

# Vectorizer defined in Step 9 above (neural network model):
# from sklearn.feature_extraction.text import CountVectorizer
# # Create a collection of the mutation names
# document = set(all_variants)
# # Create a Vectorizer Object
# vectorizer = CountVectorizer()
# vectorizer.fit(document)

import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

print("
--- Running Saliency Map Analysis for Model Interpretability ---")

# --- Step 1: Ensure we have the mutation names (features) ---
# The 'vectorizer' object holds the mapping from column index to mutation name.
try:
    feature_names = vectorizer.get_feature_names_out()
    print("Successfully retrieved feature names from vectorizer.")
except NameError:
    print("Error: 'vectorizer' object not found.")
    # You would need to re-run the cell with the CountVectorizer if this error occurs.
    # For now, we will exit if it's not available.
    feature_names = None

if feature_names is not None:
    # --- Step 2: Find one example index for each class in the test set ---
    unique_classes_in_test = np.unique(y_test)
    sample_indices = [np.where(y_test == cls)[0][0] for cls in unique_classes_in_test]

    ### --- START: MODIFICATION FOR TABLE OUTPUT --- ###

    # Create an empty list to store the results from each class
    saliency_results_list = []

    print(f"Calculating top 10 influential mutations for {len(sample_indices)} classes...")

    # --- Step 3: Loop through each sample and calculate saliency ---
    for i, sample_index in enumerate(sample_indices):
        sample_x = X_test[sample_index:sample_index+1]
        sample_y_numerical = y_test[sample_index]
        sample_y_name = lineage_map[sample_y_numerical]

        sample_tensor = tf.convert_to_tensor(sample_x, dtype=tf.float32)

        with tf.GradientTape() as tape:
            tape.watch(sample_tensor)
            predictions = model(sample_tensor)
            top_class_output = predictions[:, sample_y_numerical]

        saliency = tape.gradient(top_class_output, sample_tensor)
        saliency_scores = np.abs(saliency.numpy().flatten())

        # Create a DataFrame for this sample's results
        df_saliency = pd.DataFrame({
            'mutation': feature_names,
            'present': sample_x.flatten(),
            'importance': saliency_scores
        })

        # Filter for mutations present in this sequence and get the top 10
        top_10_features = df_saliency[df_saliency['present'] == 1].sort_values(
            by='importance', ascending=False
        ).head(10)

        # Add the lineage name to the DataFrame
        top_10_features['Lineage'] = sample_y_name

        # Append the results for this class to our master list
        saliency_results_list.append(top_10_features)

        # No plotting in the loop anymore
        # plt.figure(...)
        # plt.show()

    # --- Step 4: Combine all results into a single DataFrame ---
    if saliency_results_list:
        final_saliency_table = pd.concat(saliency_results_list, ignore_index=True)

        # Reorder columns for clarity
        final_saliency_table = final_saliency_table[['Lineage', 'mutation', 'importance']]

        print("
--- Saliency Analysis Complete ---")
        print("Top 10 most influential mutations for each class:")
        # Display the first 20 rows of the final table
        print(final_saliency_table.head(20))

        # --- Step 5: Save the table to a file ---
        csv_filename = 'saliency_report_top10.csv'
        final_saliency_table.to_csv(csv_filename, index=False)
        print(f"
Full saliency table saved to '{csv_filename}'")

    else:
        print("No saliency results were generated.")
