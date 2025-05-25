#!/usr/bin/env python3

import argparse
import re
from Bio.Seq import Seq
import os
import matplotlib
matplotlib.use('Agg') # Use Agg backend for non-interactive plotting
import matplotlib.pyplot as plt

# Add the project root directory to the Python path
# This allows imports like 'from pylib.scripts...'
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.insert(0, project_root)
from pylib.utils import seq_parser

def find_all_occurrences(sub_string, main_string):
    """Finds all starting positions of sub_string in main_string."""
    return [m.start() for m in re.finditer(re.escape(sub_string), main_string)]

def main():
    parser = argparse.ArgumentParser(description="Generate a dot plot from two FASTA sequences.")
    parser.add_argument("--seqfile1", "-1", required=True, help="Path to the first input FASTA file.")
    parser.add_argument("--seqfile2", "-2", required=True, help="Path to the second input FASTA file.")
    parser.add_argument("--wordlen", "-w", type=int, required=True, help="Word size for matching.")
    parser.add_argument("--step", "-s", type=int, required=True, help="Step length for moving the word.")
    parser.add_argument("--outfile", "-o", required=True, help="Output PNG file name.")
    parser.add_argument("--dotfile", "-d", required=True, help="Path to the intermediate file to store match coordinates.")
    parser.add_argument("--title", "-t", default="", help="Title for the plot.")
    args = parser.parse_args()

    # Read sequences
    seq1_record = None
    seq2_record = None

    try:
        for record in seq_parser.parse_fasta_file(args.seqfile1):
            seq1_record = record
            break # Take the first sequence
        if seq1_record is None:
            print(f"Error: No sequences found in {args.seqfile1}")
            return
        seq1_str = str(seq1_record.seq).upper() # Ensure uppercase for matching
    except FileNotFoundError:
        print(f"Error: File not found {args.seqfile1}")
        return
    except Exception as e:
        print(f"Error reading {args.seqfile1}: {e}")
        return

    try:
        for record in seq_parser.parse_fasta_file(args.seqfile2):
            seq2_record = record
            break # Take the first sequence
        if seq2_record is None:
            print(f"Error: No sequences found in {args.seqfile2}")
            return
        seq2_str = str(seq2_record.seq).upper() # Ensure uppercase for matching
    except FileNotFoundError:
        print(f"Error: File not found {args.seqfile2}")
        return
    except Exception as e:
        print(f"Error reading {args.seqfile2}: {e}")
        return
        
    if not seq1_str or not seq2_str:
        print("Error: One or both sequences are empty.")
        return

    # Calculate matches and write to dot file
    try:
        with open(args.dotfile, "w") as dot_fh:
            for k in range(0, len(seq1_str) - args.wordlen + 1, args.step):
                word = seq1_str[k : k + args.wordlen]

                if 'N' in word.upper(): # Check for N, case-insensitive
                    continue

                # Forward strand matches
                positions = find_all_occurrences(word, seq2_str)
                for pos in positions:
                    dot_fh.write(f"{k}\t{pos}\n")

                # Reverse complement matches
                rc_word = str(Seq(word).reverse_complement())
                rc_positions = find_all_occurrences(rc_word, seq2_str)
                for pos in rc_positions:
                    dot_fh.write(f"{k}\t{pos}\n") # Perl script stores positive y for revcomp
    except IOError:
        print(f"Error: Could not write to dotfile {args.dotfile}")
        return
    except Exception as e:
        print(f"An error occurred during match calculation: {e}")
        return

    # Generate plot
    x_coords = []
    y_coords = []
    try:
        with open(args.dotfile, "r") as dot_fh:
            for line in dot_fh:
                parts = line.strip().split("\t")
                if len(parts) == 2:
                    try:
                        x_coords.append(int(parts[0]))
                        y_coords.append(int(parts[1]))
                    except ValueError:
                        print(f"Warning: Skipping invalid line in dotfile: {line.strip()}")
                        continue
    except FileNotFoundError:
        print(f"Error: Dotfile {args.dotfile} not found for plotting.")
        return
    except Exception as e:
        print(f"Error reading dotfile for plotting: {e}")
        return

    if not x_coords: # No matches found
        print("No matches found to plot.")
        # Create an empty plot as per user expectation for an output file
        # Or, we could skip creating a file, but the Perl script likely creates one.

    fig_width_px = 700
    fig_height_px = 730
    dpi = 100 # Standard DPI for Matplotlib, can be adjusted

    fig_width_in = fig_width_px / dpi
    fig_height_in = fig_height_px / dpi
    
    fig, ax = plt.subplots(figsize=(fig_width_in, fig_height_in), dpi=dpi)

    ax.scatter(x_coords, y_coords, s=1, color='black', marker='.')
    
    plot_title = args.title if args.title else f"Dot Plot: {seq1_record.id} vs {seq2_record.id}"
    ax.set_title(f"{plot_title} -w={args.wordlen} -s={args.step}")

    # Set axis limits
    ax.set_xlim(0, len(seq1_str))
    ax.set_ylim(0, len(seq2_str))
    # The Perl script's GD drawing implies origin at top-left for sequence coordinates.
    # Matplotlib's default is bottom-left. To match the visual, if y-coords are plotted directly,
    # then larger y values are higher on the plot.
    # The Perl script maps y to $height - ($y0 + $del*$x[1]). If we plot Y as is and don't invert,
    # (0,0) seq coord is at bottom-left. If we invert, (0,0) is top-left.
    # Let's try without inverting first, as it's more standard for scatter plots.
    # ax.invert_yaxis() # Use if (0,0) should be top-left for seq coordinates

    # Draw border lines
    # ax.axvline(x=0, color='gray', linewidth=0.5)
    # ax.axvline(x=len(seq1_str), color='gray', linewidth=0.5)
    # ax.axhline(y=0, color='gray', linewidth=0.5)
    # ax.axhline(y=len(seq2_str), color='gray', linewidth=0.5)
    # The Perl script draws a box using $gd->rectangle.
    # A neater way in Matplotlib is to set spines.
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_visible(True)
    # Or, if we want to mimic the exact box:
    # rect = plt.Rectangle((0,0), len(seq1_str), len(seq2_str), fill=False, color="gray", linewidth=0.5)
    # ax.add_patch(rect)


    # Custom axis labels and remove ticks (similar to Perl script)
    ax.set_xticks([])
    ax.set_yticks([])

    # X-axis label (bottom)
    ax.text(len(seq1_str) / 2, -0.03 * len(seq2_str), f"{len(seq1_str)} bp", 
            ha='center', va='top', transform=ax.transData)
    
    # Y-axis label (left)
    ax.text(-0.03 * len(seq1_str), len(seq2_str) / 2, f"{len(seq2_str)} bp", 
            ha='right', va='center', rotation='vertical', transform=ax.transData)

    # To better utilize plot area, adjust subplot parameters
    # The Perl script has margin x0=60, y0=30, y1=60 for labels and title.
    # $del = ($width-2*$x0)/$max_len; effectively determines plot area.
    # Matplotlib's tight_layout or constrained_layout might handle this.
    plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.1) # Manual adjustment
    # plt.tight_layout(pad=3.0) # pad is in inches

    try:
        plt.savefig(args.outfile)
        print(f"Plot saved to {args.outfile}")
    except Exception as e:
        print(f"Error saving plot: {e}")

if __name__ == "__main__":
    main()
