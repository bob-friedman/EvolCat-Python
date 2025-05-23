import argparse
import time
import re
import requests

MAX_ATTEMPTS = 120 # 20 minutes (120 * 10 seconds)

def main():
    parser = argparse.ArgumentParser(description="BLAST NCBI sequence defragmenter")
    parser.add_argument("fasta_file", help="Path to the input FASTA file")
    parser.add_argument("e_value", help="Expect value threshold for BLAST")
    args = parser.parse_args()

    print(f"Input FASTA file: {args.fasta_file}")
    print(f"E-value: {args.e_value}")

    time.sleep(20) # Initial delay

    # Read FASTA file
    sequences = parse_fasta(args.fasta_file)

    for seq_id, sequence in sequences.items():
        print(f"Processing sequence: {seq_id}")
        rid = submit_blast_request(sequence, args.e_value)
        if rid == "UNKNOWN":
            print("bad input")
            exit()
        
        poll_for_results(rid)
    exit() # Ensure main always exits, facilitating SystemExit checks in tests.


def parse_fasta(fasta_file):
    sequences = {}
    current_seq_id = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if line.startswith('>'):
                if current_seq_id:
                    sequences[current_seq_id] = "".join(current_seq)
                current_seq_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line.replace(" ", ""))
        if current_seq_id: # Add the last sequence
            sequences[current_seq_id] = "".join(current_seq)
    return sequences

def submit_blast_request(sequence, e_value):
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
    }
    params = {
        'CMD': 'Put',
        'QUERY': sequence,
        'DATABASE': 'nr',
        'PROGRAM': 'blastp',
        'EXPECT': e_value,
        'HITLIST_SIZE': '1000'
    }
    try:
        response = requests.post("https://www.ncbi.nlm.nih.gov/blast/Blast.cgi", data=params, headers=headers)
        response.raise_for_status() # Raise HTTPError for bad responses (4XX or 5XX)
        
        # Extract RID
        rid_match = re.search(r"RID = (\S+)", response.text)
        if rid_match:
            return rid_match.group(1)
        else:
            # Check for "UNKNOWN" RID
            if "RID = UNKNOWN" in response.text:
                return "UNKNOWN"
            print("Could not find RID in response:")
            print(response.text)
            exit() # Exit if RID not found and not UNKNOWN
            
    except requests.exceptions.RequestException as e:
        print(f"Error submitting BLAST request: {e}")
        exit()


def poll_for_results(rid):
    attempts = 0

    while attempts < MAX_ATTEMPTS:
        attempts += 1
        time.sleep(10) # Delay between polls

        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
        }
        params = {
            'CMD': 'Get',
            'RID': rid,
            'FORMAT_TYPE': 'Text'
        }
        try:
            response = requests.post("https://www.ncbi.nlm.nih.gov/blast/Blast.cgi", data=params, headers=headers)
            response.raise_for_status()

            # Extract Status
            status_match = re.search(r"Status=(\S+)", response.text)
            if status_match:
                status = status_match.group(1)
                print(f"{attempts}\t{status}")

                if status != "WAITING":
                    print(response.text) # Print BLAST report
                    return
            else:
                print("Could not find Status in response:")
                print(response.text)
                # Decide if to continue or exit. For now, let's assume it's an error and exit.
                # If the HTML structure changes, this part might need adjustment.
                # Consider if "READY" without "Status=" line means success.
                # For now, strict check for "Status="
                if "READY" in response.text and "Status=" not in response.text : # A more lenient check if status line is missing but results seem ready
                    print("Assuming READY based on content, though Status line missing.")
                    print(response.text)
                    return
                exit()


        except requests.exceptions.RequestException as e:
            print(f"Error polling for results: {e}")
            # Decide if to continue or exit. For now, let's exit on poll error.
            exit()

    print("timed out at 20 minutes")
    exit()

if __name__ == "__main__":
    main()
