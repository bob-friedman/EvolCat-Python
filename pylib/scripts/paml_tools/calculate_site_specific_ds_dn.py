"""
Wraps PAML's codeml to perform site-specific dN/dS analysis.

This script automates the process of:
1. Preparing input files (converting alignment to PHYLIP).
2. Generating a PAML codeml control file (.ctl) based on user inputs.
3. Running codeml.
4. Parsing the codeml output.
5. Summarizing key results and site-specific dN/dS estimates,
   including Bayes Empirical Bayes (BEB) probabilities for positive selection.

Note: This script is focused on site-specific dN/dS analysis using PAML's `codeml`.
For simple pairwise dN/dS calculations between sequences, consider using the
`calculate_dn_ds.py` script in this directory, which utilizes PAML's `yn00` program.

Requires PAML (specifically the 'codeml' executable) to be installed and
either in the system PATH or its location provided via --paml_path.
"""

import argparse
import os
import sys
import subprocess
import shutil
from Bio import AlignIO, Phylo
from Bio.Phylo.PAML import codeml, CodemlError

def setup_codeml_options(cml, args):
    """Sets codeml options based on command-line arguments."""
    cml.set_options(seqtype=1)  # 1 for codons
    cml.set_options(noisy=9 if args.verbose else 3)
    cml.set_options(cleandata=args.cleandata) # 1 to remove sites with gaps/ambiguities
    cml.set_options(RateAncestor=1) # Required for BEB analysis
    
    # Common model parameters
    cml.set_options(fix_kappa=0) # Estimate kappa
    cml.set_options(kappa=2.0)   # Initial kappa
    cml.set_options(fix_omega=0) # Estimate omega
    cml.set_options(omega=0.5)   # Initial omega (background)

    # Model-specific parameters
    if args.model == "M0":
        cml.set_options(model=0)
        cml.set_options(NSsites=[0])
    elif args.model == "M1a":
        cml.set_options(model=0)
        cml.set_options(NSsites=[1])
    elif args.model == "M2a":
        cml.set_options(model=0)
        cml.set_options(NSsites=[2])
    elif args.model == "M3":
        cml.set_options(model=0)
        cml.set_options(NSsites=[3])
        cml.set_options(ncatG=args.num_site_categories)
    elif args.model == "M7":
        cml.set_options(model=0)
        cml.set_options(NSsites=[7])
    elif args.model == "M8":
        cml.set_options(model=0)
        cml.set_options(NSsites=[8])
    else:
        # This should ideally be caught by argparse choices, but as a safeguard:
        raise ValueError(f"Unsupported PAML model: {args.model}")

def main():
    parser = argparse.ArgumentParser(
        description="Run PAML codeml for site-specific dN/dS analysis.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
Requires PAML (codeml) to be installed and in PATH or specified via --paml_path.

Example Usage (using the provided sample data):
  python calculate_site_specific_ds_dn.py \\
    --alignment text_files/sample_aligned.fasta \\
    --tree your_tree_file.newick \\
    --model M2a \\
    --outfile_prefix M2a_results \\
    --verbose

Supported models:
  M0: one ratio for all sites
  M1a (neutral): nearly neutral (0 < w0 < 1) and conserved (w1 = 1) sites
  M2a (selection): adds a class of sites with w > 1 to M1a
  M3 (discrete): discrete distribution of w values (k categories)
  M7 (beta): beta distribution for w (0 < w < 1)
  M8 (beta&w>1): M7 + a class of sites with w >= 1
"""
    )
    parser.add_argument("--alignment", required=True, help="Path to the coding sequence alignment file (FASTA format).")
    parser.add_argument("--tree", required=True, help="Path to the phylogenetic tree file (Newick format).")
    parser.add_argument("--model", required=True, choices=["M0", "M1a", "M2a", "M3", "M7", "M8"],
                        help="PAML model for site-specific analysis.")
    parser.add_argument("--outfile_prefix", required=True, help="Prefix for PAML output files and the summary TSV.")
    parser.add_argument("--paml_path", help="Optional path to the codeml executable. If not provided, assumes codeml is in system PATH.")
    parser.add_argument("--verbose", action="store_true", help="Enable detailed PAML output (noisy = 9). Default is less verbose (noisy = 3).")
    parser.add_argument("--cleandata", type=int, default=1, choices=[0, 1],
                        help="PAML cleandata option: 1 to remove sites with gaps/ambiguities (default), 0 to keep.")
    parser.add_argument("--num_site_categories", type=int, default=3,
                        help="Number of site categories for M3 (discrete) model. Default: 3.")

    args = parser.parse_args()

    # --- Validate inputs ---
    if not os.path.exists(args.alignment):
        print(f"Error: Alignment file not found: {args.alignment}")
        sys.exit(1)
    if not os.path.exists(args.tree):
        print(f"Error: Tree file not found: {args.tree}")
        sys.exit(1)
    
    paml_exec = args.paml_path if args.paml_path else "codeml"
    if shutil.which(paml_exec) is None:
        print(f"Error: codeml executable not found at '{paml_exec}'. "
              "Please ensure PAML is installed and in your PATH, or provide the correct path using --paml_path.")
        sys.exit(1)

    # --- Prepare input files ---
    phylip_alignment_file = args.outfile_prefix + ".phy"
    try:
        AlignIO.convert(args.alignment, "fasta", phylip_alignment_file, "phylip-relaxed")
        print(f"Converted FASTA alignment '{args.alignment}' to PHYLIP '{phylip_alignment_file}'.")
    except Exception as e:
        print(f"Error converting alignment to PHYLIP: {e}")
        sys.exit(1)

    # --- Setup PAML codeml ---
    cml = codeml.Codeml()
    cml.alignment = phylip_alignment_file
    cml.tree = args.tree
    cml.out_file = args.outfile_prefix + ".mlc" # Main results file
    cml.working_dir = "." # Run in the current directory

    try:
        setup_codeml_options(cml, args)
    except ValueError as e:
        print(f"Error setting up PAML options: {e}")
        sys.exit(1)
    
    # Print control file content for debugging if verbose
    if args.verbose:
        print("\n--- PAML Control File (codeml.ctl preview) ---")
        # Temporarily set ctl_file to generate content for preview
        original_ctl_file_option = cml.ctl_file # Store original if set
        cml.ctl_file = "tmp_codeml_preview.ctl" 
        cml.write_ctl_file()
        with open(cml.ctl_file, "r") as f:
            print(f.read())
        os.remove(cml.ctl_file) # Clean up temp ctl
        cml.ctl_file = original_ctl_file_option # Restore original ctl_file option
        print("---------------------------------------------\n")


    # --- Run codeml ---
    print(f"Running PAML codeml with model {args.model}...")
    try:
        # Biopython's run method uses self._ctl_file attribute for the ctl filename.
        # It defaults to "codeml.ctl". Let's ensure it's unique if multiple runs occur.
        # However, cml.run() generates it based on working_dir and a standard name.
        # The main concern is if cml.run() doesn't clean it up or if it conflicts.
        # For now, assume Biopython handles ctl file naming within its execution context.
        results = cml.run(command=paml_exec, verbose=args.verbose, parse=False) 
        
        print(f"codeml execution finished. Main output: {cml.out_file}")
        
        if not os.path.exists(cml.out_file):
            print(f"Error: codeml main output file '{cml.out_file}' not found. Check for errors from codeml.")
            rst_file_path = os.path.join(cml.working_dir, "rst")
            if os.path.exists(rst_file_path):
                with open(rst_file_path, "r") as rst_f:
                    print(f"\n--- Content of '{rst_file_path}' file (may contain error details) ---")
                    print(rst_f.read(2000)) 
                    print("----------------------------------------------------------\n")
            sys.exit(1)
            
        results = codeml.read(cml.out_file)
        print("Successfully parsed codeml results.")

    except CodemlError as e:
        print(f"Error running PAML codeml: {e}")
        if hasattr(e, 'stderr') and e.stderr:
            print(f"codeml stderr:\n{e.stderr}")
        
        rst_file_path = os.path.join(cml.working_dir, "rst") 
        rub_file_path = os.path.join(cml.working_dir, "rub")
        if os.path.exists(rst_file_path):
            with open(rst_file_path, "r") as rst_f:
                print(f"\n--- Content of '{rst_file_path}' file (may contain error details) ---")
                print(rst_f.read(2000))
                print("----------------------------------------------------------\n")
        if os.path.exists(rub_file_path):
             with open(rub_file_path, "r") as rub_f:
                print(f"\n--- Content of '{rub_file_path}' file (may contain error details) ---")
                print(rub_f.read(2000))
                print("----------------------------------------------------------\n")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred during codeml execution or parsing: {e}")
        sys.exit(1)

    # --- Parse and Summarize Output ---
    print("\n--- PAML codeml Results Summary ---")
    print(f"Model: {args.model}")
    print(f"Log-likelihood (lnL): {results.get('lnL', 'N/A')}")

    if 'parameters' in results:
        params = results['parameters']
        if 'kappa' in params:
            print(f"kappa (ts/tv ratio): {params['kappa']:.4f}")
        
        ns_sites_key = str(cml.get_option("NSsites")) # e.g., "[0]", "[1]", "[2]"
        
        if ns_sites_key in results.get("NSsites", {}):
            ns_params = results["NSsites"][ns_sites_key].get("parameters", {})
            if "site classes" in ns_params: # For M1a, M2a, M3, M7, M8
                print("Site classes:")
                for sc in ns_params["site classes"]:
                    proportion = sc.get('proportion', {}).get('value', 'N/A')
                    # Omega value key can vary, sometimes 'omega', sometimes 'M0' (for the omega value itself)
                    omega_val = sc.get('omega', {}).get('M0', sc.get('omega', 'N/A')) 
                    if isinstance(omega_val, dict): # If omega is still a dict, try common keys
                        omega_val = omega_val.get('M0', omega_val.get('omega', 'N/A'))

                    print(f"  - Proportion: {proportion:.4f}, omega: {omega_val:.4f}")
            elif 'omega' in ns_params: # For M0, if parsed into NSsites structure
                 print(f"omega: {ns_params['omega']:.4f}")
        elif args.model == "M0" and 'omega' in params: # M0 omega might be directly in parameters
             print(f"omega: {params['omega']:.4f}")


    site_analysis_file = args.outfile_prefix + "_site_analysis.tsv"
    print(f"\nWriting site-specific analysis to: {site_analysis_file}")

    with open(site_analysis_file, "w") as tsv_out:
        header = ["Site", "AminoAcid", "dN_dS", "PosteriorProbability_PositiveSelection", "Note"]
        tsv_out.write("\t".join(header) + "\n")
        
        beb_sites_data = []
        all_sites_posterior_data = []
        
        ns_sites_key_str = str(cml.get_option("NSsites"))
        
        if ns_sites_key_str in results.get("NSsites", {}):
            model_specific_params = results["NSsites"][ns_sites_key_str].get("parameters", {})

            # Extract site class omegas
            site_class_omegas = {}
            if "site classes" in model_specific_params:
                for i, sc_detail in enumerate(model_specific_params["site classes"]):
                    omega_val = sc_detail.get('omega', {}).get('M0', sc_detail.get('omega', float('nan')))
                    if isinstance(omega_val, dict): # If omega is still a dict, try common keys
                        omega_val = omega_val.get('M0', omega_val.get('omega', float('nan')))
                    site_class_omegas[i] = omega_val
            
            positive_selection_class_idx = -1
            positive_selection_omega_value = float('nan')
            for class_idx, omega_val in site_class_omegas.items():
                if omega_val > 1.0:
                    positive_selection_class_idx = class_idx
                    positive_selection_omega_value = omega_val
                    break # Assuming only one such class for M2a/M8

            # For M2a and M8, BEB results are explicitly listed for positively selected sites
            if args.model in ["M2a", "M8"] and "BEB positive sites" in model_specific_params:
                for site_info in model_specific_params["BEB positive sites"]:
                    # Structure: (site_index_1_based, amino_acid, prob_of_belonging_to_pos_sel_class, omega_of_pos_sel_class)
                    # Or sometimes: (site_index_1_based, amino_acid, prob_class_k, omega_k) where class k has w > 1
                    site_idx_1based, aa, prob_pos_sel_class = site_info[0], site_info[1], site_info[2]
                    # The dN/dS value might be the specific omega for that class, or we use the identified positive_selection_omega_value
                    dnds_val_beb = positive_selection_omega_value if not (len(site_info) > 3 and isinstance(site_info[3], float)) else site_info[3]
                    
                    note = ""
                    if prob_pos_sel_class > 0.95: note = "Positive Selection (BEB > 0.95)"
                    elif prob_pos_sel_class > 0.90: note = "Positive Selection (BEB > 0.90)"
                    elif prob_pos_sel_class > 0.50: note = "Potential Positive Selection (BEB > 0.50)"
                    
                    beb_sites_data.append({
                        "site": site_idx_1based, "aa": aa, "dnds": dnds_val_beb,
                        "prob_positive": prob_pos_sel_class, "note": note
                    })

            # General site posterior probabilities (for all models with NSsites, if RateAncestor=1)
            # Structure: (site_num_1_based, aa, [prob_class0, prob_class1, ...])
            if "site_posterior_prob" in model_specific_params:
                for site_info in model_specific_params["site_posterior_prob"]:
                    site_idx_1based, aa = site_info[0], site_info[1]
                    probs_per_class = site_info[2:] # This should be a list of probabilities

                    # Determine the most probable class and its omega
                    most_probable_class_idx = probs_per_class.index(max(probs_per_class))
                    dnds_for_site = site_class_omegas.get(most_probable_class_idx, float('nan'))
                    
                    prob_positive_selection_class = 0.0
                    note_for_site = ""
                    if positive_selection_class_idx != -1 and positive_selection_class_idx < len(probs_per_class):
                        prob_positive_selection_class = probs_per_class[positive_selection_class_idx]
                        if prob_positive_selection_class > 0.95: note_for_site = "Positive Selection (BEB > 0.95)"
                        elif prob_positive_selection_class > 0.90: note_for_site = "Positive Selection (BEB > 0.90)"
                        elif prob_positive_selection_class > 0.50: note_for_site = "Potential Positive Selection (BEB > 0.50)"
                    
                    all_sites_posterior_data.append({
                        "site": site_idx_1based, "aa": aa, "dnds": dnds_for_site,
                        "prob_positive": prob_positive_selection_class, "note": note_for_site
                    })

        # Output logic:
        # For M2a/M8, if beb_sites_data is populated, it's preferred as it's specific.
        # Otherwise, or for other models, use all_sites_posterior_data if available.
        
        final_output_data = []
        if args.model in ["M2a", "M8"] and beb_sites_data:
            final_output_data = beb_sites_data
        elif all_sites_posterior_data:
            # For models not M2a/M8, or if BEB list was empty for M2a/M8,
            # use the comprehensive list. The "prob_positive" will be 0 if no class has w > 1.
            final_output_data = all_sites_posterior_data
        
        if final_output_data:
            final_output_data.sort(key=lambda x: x['site'])
            for site_entry in final_output_data:
                 tsv_out.write(f"{site_entry['site']}\t{site_entry['aa']}\t"
                               f"{site_entry['dnds']:.4f}\t{site_entry['prob_positive']:.4f}\t{site_entry['note']}\n")
        else:
            print("Warning: Could not extract detailed site-specific dN/dS data or BEB probabilities. "
                  "The TSV file will be empty or incomplete. This might be expected for model M0, "
                  "or indicates an issue with PAML output parsing or lack of positive selection findings.")
            if args.model == "M0":
                 tsv_out.write("Site-specific dN/dS is not applicable for model M0 (single ratio for all sites).\n")


    # --- Cleanup temporary files ---
    # phylip_alignment_file is <outfile_prefix>.phy, PAML output is <outfile_prefix>.mlc
    # These are named based on user prefix, so generally should be kept.
    # PAML also creates files like 2NG.dN, 2NG.dS, rst, rub in the working_dir.
    # Users might want to inspect these. We will not delete them by default.
    
    print(f"\nScript finished. Check '{cml.out_file}' and '{site_analysis_file}' for results.")
    print(f"Other PAML output files (e.g., rst, rub, 2NG.dN, 2NG.dS) may also be present in '{cml.working_dir}'.")

if __name__ == "__main__":
    main()
