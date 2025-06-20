"""
End-to-End SARS-CoV-2 Lineage Classification Pipeline

This script implements a full workflow for:
1. Setting up a Conda environment with necessary bioinformatics tools in Google Colab.
2. Downloading and preparing SARS-CoV-2 data (UShER MAT, reference genome).
3. Extracting data for a specific clade of interest into a VCF file.
4. Performing feature engineering to create a (samples x mutations) matrix.
5. Training and evaluating a neural network to classify fine-grained lineages.
6. Generating interpretability reports (saliency maps) to identify key mutations.

Designed for execution in Google Colab, but can be adapted for other environments.
"""

import os
import gc
import subprocess
import pandas as pd
import numpy as np
import pysam
import tensorflow as tf
from collections import defaultdict
from google.colab import drive
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, classification_report
from sklearn.feature_extraction.text import CountVectorizer
import matplotlib.pyplot as plt
import seaborn as sns

# ==============================================================================
# --- 1. CONFIGURATION ---
# All user-configurable parameters and file paths are here.
# ==============================================================================

# --- Key Parameters ---
CLADE_OF_INTEREST = 'XBB.2.3' # The target clade for analysis
N_TOP_CLASSES_VIZ = 15        # Number of classes to show in the confusion matrix

# --- File Paths (relative to Google Drive root) ---
# Assumes the script is run from '/content/drive/MyDrive'
DRIVE_MOUNT_POINT = '/content/drive'
WORKING_DIRECTORY = os.path.join(DRIVE_MOUNT_POINT, 'MyDrive')

# Input data URLs
MAT_URL = "http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.all.masked.pb.gz"
REF_SEQ_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/chromosomes/NC_045512v2.fa.gz"

# Local file names
MAT_GZ = "public-latest.all.masked.pb.gz"
MAT_PB = "public-latest.all.masked.pb"
METADATA_TSV = 'public-latest.metadata.tsv'
REF_SEQ_GZ = "NC_045512v2.fa.gz"
REF_SEQ_FA = "NC_045512v2.fa"

# Generated files
CLADE_VCF = 'my_clade.vcf'
CLADE_VCF_GZ = f'{CLADE_VCF}.gz'
CONSENSUS_FA = 'consensus.fa'
CONFUSION_MATRIX_PNG = f'confusion_matrix_top_{N_TOP_CLASSES_VIZ}.png'
SALIENCY_REPORT_CSV = 'saliency_report_top10.csv'


def run_command(command):
    """Executes a shell command and raises an error if it fails."""
    print(f"Executing: {' '.join(command)}")
    subprocess.run(command, check=True)


def setup_environment():
    """Installs all necessary dependencies using Conda and Pip."""
    print("\n--- [STAGE 1/6] SETTING UP ENVIRONMENT ---")
    
    # Mount Google Drive
    drive.mount(DRIVE_MOUNT_POINT)
    if not os.path.exists(WORKING_DIRECTORY):
        os.makedirs(WORKING_DIRECTORY)
    os.chdir(WORKING_DIRECTORY)
    print(f"Working directory set to: {os.getcwd()}")
    
    # Install Conda
    print("Installing Conda...")
    run_command(["pip", "install", "-q", "condacolab"])
    import condacolab
    condacolab.install()
    run_command(["conda", "init", "bash"])

    # Install tools
    print("Installing UShER, bcftools, and PySam via Conda...")
    run_command(["conda", "config", "--add", "channels", "defaults"])
    run_command(["conda", "config", "--add", "channels", "bioconda"])
    run_command(["conda", "config", "--add", "channels", "conda-forge"])
    run_command(["conda", "install", "-q", "-y", "usher", "bcftools", "pysam"])
    print("Environment setup complete.")


def prepare_data_files():
    """Downloads, unzips, and prepares all required data files."""
    print("\n--- [STAGE 2/6] PREPARING DATA FILES ---")

    # Download and unzip MAT file
    if not os.path.isfile(MAT_PB):
        print(f"Downloading MAT file from {MAT_URL}...")
        run_command(["wget", "-O", MAT_GZ, MAT_URL])
        print("Unzipping MAT file...")
        run_command(["gunzip", "-f", MAT_GZ])
    else:
        print("MAT file already exists.")

    # Download and unzip reference sequence
    if not os.path.isfile(REF_SEQ_FA):
        print(f"Downloading reference sequence from {REF_SEQ_URL}...")
        run_command(["wget", "-O", REF_SEQ_GZ, REF_SEQ_URL])
        print("Unzipping reference sequence...")
        run_command(["gunzip", "-f", REF_SEQ_GZ])
    else:
        print("Reference sequence already exists.")
        
    # Extract clade-specific VCF
    print(f"Extracting clade '{CLADE_OF_INTEREST}' to VCF...")
    run_command([
        "matUtils", "extract",
        "-i", MAT_PB,
        "-c", CLADE_OF_INTEREST,
        "-v", CLADE_VCF
    ])
    
    # Compress and index VCF
    print("Compressing and indexing VCF file...")
    run_command(["bgzip", "--force", CLADE_VCF])
    run_command(["bcftools", "index", "--tbi", CLADE_VCF_GZ])
    print("Data file preparation complete.")


def feature_engineering():
    """Loads data, builds, and returns the feature matrix (X) and labels (y)."""
    print("\n--- [STAGE 3/6] PERFORMING FEATURE ENGINEERING ---")
    
    # Step 1: Load metadata and filter by samples in our VCF
    print("[1/5] Loading metadata and VCF sample list...")
    try:
        metadata_df = pd.read_csv(METADATA_TSV, sep='\t', low_memory=False)
        metadata_df.dropna(subset=['strain'], inplace=True)
    except FileNotFoundError:
        print(f"ERROR: Metadata file '{METADATA_TSV}' not found in {os.getcwd()}.")
        raise

    vcf_in = pysam.VariantFile(CLADE_VCF_GZ)
    samples_in_vcf = list(vcf_in.header.samples)
    vcf_in.close()
    print(f"Found {len(samples_in_vcf)} samples in the VCF file.")
    
    clade_metadata_df = metadata_df[metadata_df['strain'].isin(samples_in_vcf)].copy()
    del metadata_df
    gc.collect()

    if clade_metadata_df.empty:
        raise ValueError("No samples from the VCF were found in the metadata.")

    # Step 2: Identify all unique variants to create features
    print("[2/5] Identifying all unique variants...")
    vcf_in = pysam.VariantFile(CLADE_VCF_GZ)
    all_variants = sorted(list(set(
        f"{rec.pos}_{rec.ref}>{rec.alts[0]}" for rec in vcf_in.fetch() if rec.alts
    )))
    vcf_in.close()
    print(f"Found {len(all_variants)} unique variants to use as features.")

    # Step 3: Map variants to samples that contain them (for efficiency)
    print("[3/5] Mapping variants to samples...")
    variant_to_samples_map = defaultdict(list)
    vcf_in = pysam.VariantFile(CLADE_VCF_GZ)
    for record in vcf_in.fetch():
        if not record.alts: continue
        variant_id = f"{record.pos}_{record.ref}>{record.alts[0]}"
        for sample in record.samples.values():
            if 1 in sample['GT']:
                variant_to_samples_map[variant_id].append(sample.name)
    vcf_in.close()

    # Step 4: Build the feature matrix X
    print("[4/5] Building the feature matrix X...")
    X = pd.DataFrame(0, index=samples_in_vcf, columns=all_variants, dtype='int8')
    for variant, samples in variant_to_samples_map.items():
        valid_samples = [s for s in samples if s in X.index]
        if valid_samples:
            X.loc[valid_samples, variant] = 1

    # Step 5: Prepare and filter labels (y), removing singletons
    print("[5/5] Preparing and filtering labels y...")
    labels_df = clade_metadata_df.set_index('strain')
    class_counts = labels_df['pango_lineage_usher'].value_counts()
    classes_to_keep = class_counts[class_counts >= 2].index
    df_filtered = labels_df[labels_df['pango_lineage_usher'].isin(classes_to_keep)]
    
    # Align X and y perfectly
    X_filtered = X.loc[df_filtered.index]
    y_filtered_labels = df_filtered['pango_lineage_usher']
    y_numerical, lineage_map = pd.factorize(y_filtered_labels)
    
    print("--- Feature Engineering Complete! ---")
    print(f"Final matrix shape: {X_filtered.shape}")
    print(f"Final labels shape: {y_numerical.shape}")
    
    return X_filtered, y_numerical, lineage_map, all_variants


def train_and_evaluate_model(X_filtered, y_numerical, lineage_map, all_variants):
    """Trains a neural network, evaluates it, and visualizes performance."""
    print("\n--- [STAGE 4/6] BUILDING AND TRAINING NEURAL NETWORK ---")

    # Step 1: Prepare data for Keras
    X_np = X_filtered.to_numpy()
    y_np = y_numerical
    
    vectorizer = CountVectorizer()
    vectorizer.fit(all_variants)

    # Step 2: Stratified train-test split
    X_train, X_test, y_train, y_test = train_test_split(
        X_np, y_np, test_size=0.2, random_state=42, stratify=y_np
    )

    # Step 3: Define, compile, and train the model
    num_features = X_train.shape[1]
    num_classes = len(np.unique(y_np))

    model = tf.keras.models.Sequential([
        tf.keras.layers.Dense(128, activation='relu', input_dim=num_features, name='input_layer'),
        tf.keras.layers.Dense(64, activation='relu', name='hidden_layer_1'),
        tf.keras.layers.Dense(num_classes, activation='softmax', name='output_layer')
    ])
    
    model.compile(
        optimizer='adam',
        loss='sparse_categorical_crossentropy',
        metrics=['accuracy']
    )
    
    print("Model Architecture:")
    model.summary()

    print("\nTraining the model...")
    history = model.fit(
        X_train, y_train,
        validation_data=(X_test, y_test),
        epochs=10,
        batch_size=32,
        verbose=1
    )

    # Step 4: Evaluate and visualize
    print("\n--- [STAGE 5/6] EVALUATING MODEL AND VISUALIZING RESULTS ---")
    loss, accuracy = model.evaluate(X_test, y_test, verbose=0)
    print(f"\nFinal Test Accuracy: {accuracy * 100:.2f}%")
    
    # Confusion Matrix
    y_pred_probs = model.predict(X_test)
    y_pred_numerical = np.argmax(y_pred_probs, axis=1)
    
    class_counts = pd.Series(y_test).value_counts()
    top_n_indices = class_counts.head(N_TOP_CLASSES_VIZ).index.tolist()
    top_n_names = [lineage_map[i] for i in top_n_indices]

    cm_full = confusion_matrix(y_test, y_pred_numerical)
    cm_filtered = cm_full[top_n_indices, :][:, top_n_indices]

    plt.figure(figsize=(12, 10))
    sns.heatmap(cm_filtered, annot=True, fmt='d', cmap='Blues',
                xticklabels=top_n_names, yticklabels=top_n_names)
    plt.title(f'Confusion Matrix (Top {N_TOP_CLASSES_VIZ} Most Frequent Lineages)', fontsize=16)
    plt.ylabel('True Lineage', fontsize=12)
    plt.xlabel('Predicted Lineage', fontsize=12)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(CONFUSION_MATRIX_PNG, dpi=150)
    plt.show()

    # Classification Report
    print("\nFull Classification Report:")
    present_labels = np.unique(np.concatenate((y_test, y_pred_numerical)))
    filtered_target_names = [lineage_map[i] for i in present_labels]
    print(classification_report(
        y_test, y_pred_numerical,
        labels=present_labels, target_names=filtered_target_names, zero_division=0
    ))

    return model, X_test, y_test, vectorizer, lineage_map


def generate_saliency_report(model, X_test, y_test, vectorizer, lineage_map):
    """Calculates saliency scores and saves a report of influential mutations."""
    print("\n--- [STAGE 6/6] GENERATING SALIENCY REPORT FOR INTERPRETABILITY ---")

    feature_names = vectorizer.get_feature_names_out()
    saliency_results = []
    
    unique_classes_in_test = np.unique(y_test)
    sample_indices = [np.where(y_test == cls)[0][0] for cls in unique_classes_in_test]
    
    print(f"Calculating top 10 influential mutations for {len(sample_indices)} classes...")

    for sample_index in sample_indices:
        sample_x = X_test[sample_index:sample_index+1]
        sample_y_num = y_test[sample_index]
        sample_y_name = lineage_map[sample_y_num]
        
        sample_tensor = tf.convert_to_tensor(sample_x, dtype=tf.float32)

        with tf.GradientTape() as tape:
            tape.watch(sample_tensor)
            predictions = model(sample_tensor)
            top_class_output = predictions[:, sample_y_num]
        
        saliency_scores = np.abs(tape.gradient(top_class_output, sample_tensor).numpy().flatten())
        
        df_saliency = pd.DataFrame({
            'mutation': feature_names,
            'present_in_sample': sample_x.flatten(),
            'importance': saliency_scores
        })
        
        top_10 = df_saliency[df_saliency['present_in_sample'] == 1].sort_values(
            by='importance', ascending=False
        ).head(10)
        
        top_10['Lineage'] = sample_y_name
        saliency_results.append(top_10)
        
    if saliency_results:
        final_report = pd.concat(saliency_results, ignore_index=True)
        final_report = final_report[['Lineage', 'mutation', 'importance']]
        final_report.to_csv(SALIENCY_REPORT_CSV, index=False)
        print(f"\nSaliency analysis complete. Report saved to '{SALIENCY_REPORT_CSV}'")
        print("Top 20 rows of the report:")
        print(final_report.head(20))
    else:
        print("Could not generate saliency report.")


def main():
    """Main function to run the entire pipeline."""
    setup_environment()
    prepare_data_files()
    X_filtered, y_numerical, lineage_map, all_variants = feature_engineering()
    model, X_test, y_test, vectorizer, lineage_map = train_and_evaluate_model(
        X_filtered, y_numerical, lineage_map, all_variants
    )
    generate_saliency_report(model, X_test, y_test, vectorizer, lineage_map)
    print("\n--- PIPELINE EXECUTION COMPLETE ---")


if __name__ == "__main__":
    main()
