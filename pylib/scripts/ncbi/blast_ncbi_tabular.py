#!/usr/bin/env python3

"""
Submits a sequence to NCBI BLAST, polls for results, and prints them in a tabular format.
"""

import argparse
import requests
import time
import re
import sys
import os

# Add the parent directory of 'ncbi' to sys.path to allow importing utility_functions
# This assumes the script is in 'pylib/scripts/ncbi/' and 'utility_functions.py' is in 'pylib/utils/'
current_dir = os.path.dirname(os.path.abspath(__file__))
pylib_scripts_dir = os.path.dirname(current_dir)
pylib_dir = os.path.dirname(pylib_scripts_dir)
if pylib_dir not in sys.path:
    sys.path.insert(0, pylib_dir)

try:
    from utils import utility_functions
except ImportError:
    # Fallback if the above structure isn't perfect, try direct relative if utils is in scripts
    try:
        # This path might be needed if utility_functions is directly in pylib/utils
        # and the script is run from a different working directory.
        # For a package structure, direct imports like 'from ..utils import utility_functions'
        # would be preferred if the package is installed or PYTHONPATH is set.
        # However, for script execution as seen in tests, this dynamic path adjustment is common.
        # Assuming utility_functions.py is in pylib/utils
        utils_path = os.path.join(pylib_dir, 'utils')
        if utils_path not in sys.path:
            sys.path.insert(0, utils_path) # Add utils directory to path
        import utility_functions
    except ImportError:
        print("Error: Unable to import utility_functions. Ensure it's in the correct utils directory relative to scripts.", file=sys.stderr)
        sys.exit(1)


# Default NCBI BLAST URL
DEFAULT_NCBI_URL = "https://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
USER_AGENT = utility_functions.get_user_agent()


def submit_blast_query(sequence, database, program, filter_opt, hitlist_size, evalue, ncbi_url):
    """
    Submits the BLAST query to NCBI.
    Returns the Request ID (RID) if successful, otherwise None.
    """
    params = {
        'CMD': 'Put',
        'QUERY': sequence,
        'DATABASE': database,
        'PROGRAM': program,
        'FILTER': filter_opt,
        'HITLIST_SIZE': hitlist_size,
        'EXPECT': evalue,
        'FORMAT_OBJECT': 'Alignment', # Request alignment object
        'FORMAT_TYPE': 'Text' # To get RID in plain text
    }
    try:
        print(f"Submitting BLAST query to {ncbi_url}...", file=sys.stderr)
        print(f"Params: {params}", file=sys.stderr)
        response = requests.post(ncbi_url, data=params, headers={'User-Agent': USER_AGENT})
        response.raise_for_status() # Raise HTTPError for bad responses (4XX or 5XX)
        
        content = response.text
        print(f"NCBI Response Content for RID extraction:\n{content[:1000]}...", file=sys.stderr) # Log part of response

        # Try to extract RID using regex, similar to blast_ncbi_seq_def.py
        # Look for name="RID" value="YOUR_RID_HERE"
        rid_match_html = re.search(r'name="RID" value="(\S+)"', content)
        if rid_match_html:
            rid = rid_match_html.group(1)
            print(f"Extracted RID (HTML): {rid}", file=sys.stderr)
            return rid

        # Look for lines like RID = YOUR_RID_HERE
        rid_match_text = re.search(r'^\s*RID\s*=\s*(\S+)', content, re.MULTILINE)
        if rid_match_text:
            rid = rid_match_text.group(1)
            print(f"Extracted RID (text): {rid}", file=sys.stderr)
            return rid
        
        # Fallback: check for QBlastInfo_RID
        qblast_rid_match = re.search(r'QBlastInfoBegin\s*^\s*RID\s*=\s*(\S+)\s*QBlastInfoEnd', content, re.MULTILINE | re.DOTALL)
        if qblast_rid_match:
            rid = qblast_rid_match.group(1)
            print(f"Extracted RID (QBlastInfo): {rid}", file=sys.stderr)
            return rid

        print("Error: Could not find RID in NCBI response.", file=sys.stderr)
        print(f"Full response text was:\n{content}", file=sys.stderr)
        return None

    except requests.exceptions.RequestException as e:
        print(f"Error submitting BLAST request: {e}", file=sys.stderr)
        return None


def poll_blast_results(rid, poll_interval, max_poll_attempts, ncbi_url):
    """
    Polls NCBI for BLAST results using the RID.
    Returns the raw tabular BLAST report text if successful, otherwise None.
    """
    print(f"Polling for results with RID: {rid}", file=sys.stderr)
    attempts = 0
    while attempts < max_poll_attempts:
        attempts += 1
        print(f"Polling attempt {attempts}/{max_poll_attempts}...", file=sys.stderr)
        time.sleep(poll_interval)

        params = {
            'CMD': 'Get',
            'RID': rid,
            'ALIGNMENT_VIEW': 'Tabular',
            'FORMAT_TYPE': 'Text' # Ensure we get plain text tabular
        }
        try:
            response = requests.post(ncbi_url, data=params, headers={'User-Agent': USER_AGENT}) # Changed to POST
            response.raise_for_status()
            content = response.text

            # Check status
            if "Status=WAITING" in content:
                print("Status: WAITING", file=sys.stderr)
                status_match = re.search(r"Status=WAITING", content)
                if status_match: # More robust check
                    delay_match = re.search(r"کند.(\d+)", content) # Example from Perl: "کند.(10)" -> wait 10s
                    if delay_match:
                        try:
                            specific_delay = int(delay_match.group(1))
                            print(f"NCBI suggests specific delay: {specific_delay}s. Overriding poll_interval.", file=sys.stderr)
                            time.sleep(specific_delay)
                        except ValueError:
                            pass # Ignore if not an int
                continue 
            elif "Status=FAILED" in content:
                print("Error: BLAST query failed.", file=sys.stderr)
                print(content, file=sys.stderr)
                return None
            elif "Status=UNKNOWN" in content:
                print("Error: BLAST query RID expired or unknown.", file=sys.stderr)
                print(content, file=sys.stderr)
                return None
            elif "Status=READY" in content:
                print("Status: READY", file=sys.stderr)
                # Further check if results are actually there
                if "No hits found" in content:
                    print("No hits found for the query.", file=sys.stdout)
                    return "" # Return empty string to signify no hits but successful run
                # Heuristic: If it says READY and doesn't say "No hits found", assume results are in the content
                # The actual tabular data might be preceded by other info.
                return content
            else:
                # If status is not clear, but doesn't seem to be an error,
                # check for tabular data markers.
                if content.strip().startswith("# BLAST") and "# Fields:" in content:
                    print("Status: Assumed READY (tabular data found without explicit Status=READY).", file=sys.stderr)
                    return content
                elif "No hits found" in content: # Should be caught by READY block too, but good fallback
                    print("No hits found for the query (detected without explicit READY).", file=sys.stdout)
                    return ""

                print(f"Warning: Unknown status or unexpected response format during polling (attempt {attempts}):", file=sys.stderr)
                print(content[:1000] + "...", file=sys.stderr) # Print a snippet
                # Continue polling for a few more times in case it's a transient issue
                if attempts < max_poll_attempts / 2 : # Only continue if less than half attempts used
                    continue
                else:
                    print("Error: Too many attempts with unclear status. Aborting.", file=sys.stderr)
                    return None


        except requests.exceptions.RequestException as e:
            print(f"Error polling for results: {e}", file=sys.stderr)
            # Don't immediately exit, could be a transient network issue
            if attempts >= max_poll_attempts:
                return None # Exit if max attempts reached
            print("Retrying...", file=sys.stderr)


    print("Error: Timed out waiting for BLAST results.", file=sys.stderr)
    return None


def parse_subject_id(raw_subject_id):
    """
    Parses the raw subject ID string into a more readable format.
    Example: 'gb|AEG74000.1|emb|CAM39480.1|' -> 'gb:AEG74000.1 emb:CAM39480.1'
    Handles cases with fewer parts as well.
    """
    if not raw_subject_id:
        return ""
    
    # Remove leading/trailing pipes and then split
    parts = raw_subject_id.strip('|').split('|')
    
    processed_parts = []
    for i in range(0, len(parts) -1, 2): # Iterate in steps of 2
        try:
            processed_parts.append(f"{parts[i]}:{parts[i+1]}")
        except IndexError:
            # If there's an odd number of parts, the last one might be a standalone identifier
            processed_parts.append(parts[i]) 
            
    return " ".join(processed_parts)


def parse_and_print_tabular_results(report_text):
    """
    Parses the raw tabular BLAST report and prints selected fields.
    """
    if not report_text.strip(): # Handles case where poll_blast_results returns "" for "No hits found"
        print("No results to parse.", file=sys.stderr)
        return

    lines = report_text.splitlines()
    
    header_line = None
    header_index = -1
    for i, line in enumerate(lines):
        if line.startswith("# Fields:"):
            header_line = line
            header_index = i
            break
            
    if not header_line:
        # Fallback: if # Fields: is missing, but we have lines that look like data
        # and the report_text itself contains typical BLAST headers like "# BLASTP"
        # This is a less ideal scenario.
        if any(line.startswith("# BLAST") for line in lines) and \
           any(not line.startswith("#") and len(line.split('\t')) > 10 for line in lines):
            print("Warning: '# Fields:' line not found. Attempting to parse data assuming standard columns.", file=sys.stderr)
            # Default headers for common BLAST tabular output (12 columns)
            # query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
            headers = ["query id", "subject id", "% identity", "alignment length", "mismatches", "gap opens", 
                       "q. start", "q. end", "s. start", "s. end", "evalue", "bit score"]
        else:
            print("Error: Could not find '# Fields:' header line in BLAST report.", file=sys.stderr)
            print("--- Report Snippet ---", file=sys.stderr)
            for i, l in enumerate(lines[:20]): # Print first 20 lines for debugging
                print(f"Line {i}: {l}", file=sys.stderr)
            print("--- End Snippet ---", file=sys.stderr)
            return
    else:
        # Clean up header: remove '# Fields: ' and split by ', '
        headers = [h.strip() for h in header_line.replace("# Fields: ", "").split(', ')]

    # Determine column indices (0-based)
    # subject_id_col = 1 (subject acc.ver)
    # identity_col = 2 (% identity)
    # length_col = 3 (alignment length)
    # evalue_col = 10 (evalue)
    # These are standard for default tabular output.
    # If headers were dynamically parsed, one could search for them by name.
    # For simplicity, using fixed indices as per typical BLAST output.

    print(f"{'#':<3} {'Subject ID (Processed)':<30} {'Identity':<8} {'Length':<7} {'E-value':<10}")
    
    count = 0
    for i, line in enumerate(lines):
        if i <= header_index: # Skip lines at or before the header line
            continue
        if line.startswith("#"): # Skip any other comment lines after headers
            continue
        if not line.strip(): # Skip empty lines
            continue

        fields = line.split('\t')
        if len(fields) < 11: # Need at least 11 columns for e-value
            print(f"Warning: Skipping malformed line (not enough columns): {line}", file=sys.stderr)
            continue
        
        count += 1
        raw_subject_id = fields[1]
        processed_subject_id = parse_subject_id(raw_subject_id)
        identity = fields[2]
        align_length = fields[3]
        evalue = fields[10]

        # Fixed-width printing
        # Using f-strings for alignment:
        # {value:<width} for left-align
        # {value:>width} for right-align
        # {value:^width} for center-align
        print(f"{count:<3} {processed_subject_id:<30.30} {identity:>8} {align_length:>7} {evalue:>10}")


def main():
    parser = argparse.ArgumentParser(description="Submit a sequence to NCBI BLAST and get tabular results.")
    parser.add_argument("--query", required=True, help="Input sequence (string).")
    parser.add_argument("--database", default="nr", help="BLAST database (default: nr).")
    parser.add_argument("--program", default="blastp", help="BLAST program (default: blastp).")
    parser.add_argument("--filter", default="L", help="Filter option (default: L for low complexity). Use 'F' for no filter, or 'm S' for mask by sequence ID.")
    parser.add_argument("--hitlist_size", type=int, default=20, help="Number of hits to return (default: 20).")
    parser.add_argument("--evalue", type=float, default=0.01, help="E-value threshold (default: 0.01).")
    parser.add_argument("--poll_interval", type=int, default=10, help="Seconds between polls (default: 10).")
    parser.add_argument("--max_poll_attempts", type=int, default=120, help="Max poll attempts before timeout (default: 120, e.g., 120*10s = 20 minutes).")
    parser.add_argument("--ncbi_url", default=DEFAULT_NCBI_URL, help=f"NCBI BLAST URL (default: {DEFAULT_NCBI_URL}).")
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
        
    args = parser.parse_args()

    # Initial sleep as per original Perl script's implicit behavior (e.g. before first poll)
    # However, the main polling loop already has a poll_interval sleep *before* the first check.
    # The original script had a fixed 20s initial sleep in blast_ncbi_seq_def.pl
    # Adding a small initial delay before submission, though NCBI might not require this.
    # time.sleep(2) # Small courtesy delay

    rid = submit_blast_query(
        args.query, args.database, args.program, args.filter,
        args.hitlist_size, args.evalue, args.ncbi_url
    )

    if not rid:
        print("Exiting due to error in BLAST submission.", file=sys.stderr)
        sys.exit(1)
    
    # A short delay after RID retrieval before starting to poll, similar to blast_ncbi_seq_def.py's initial 20s.
    # This can prevent hitting the server too quickly.
    print(f"Obtained RID: {rid}. Waiting {args.poll_interval}s before first poll status check.", file=sys.stderr)
    # The poll_blast_results function itself will sleep for poll_interval *before* the first check.

    results_text = poll_blast_results(rid, args.poll_interval, args.max_poll_attempts, args.ncbi_url)

    if results_text is None: # Indicates an error or timeout during polling
        print("Exiting due to error or timeout during polling.", file=sys.stderr)
        sys.exit(1)
    elif not results_text.strip() and "No hits found" not in results_text: # Check if it was an empty result from "No hits found"
        # This condition might be redundant if "No hits found" cases return "" and are handled by parse_and_print
        print("No results returned from BLAST query or results were empty.", file=sys.stderr)
        # sys.exit(0) # Successfully determined no hits, or empty result.
        # Let parse_and_print_tabular_results handle empty string.
    
    parse_and_print_tabular_results(results_text)
    print("\nBLAST processing finished.", file=sys.stderr)

if __name__ == "__main__":
    main()
