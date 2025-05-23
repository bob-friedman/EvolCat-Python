import argparse
import requests
import re
import time

SEARCH_URL = "https://www.ncbi.nlm.nih.gov/entrez/query.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
USER_AGENT = 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'

def main():
    parser = argparse.ArgumentParser(description="Query NCBI Nucleotide database and retrieve GenBank records.")
    parser.add_argument("search_term", help="Term to search for in NCBI's Nucleotide database.")
    args = parser.parse_args()

    print(f"Search Term: {args.search_term}")

    try:
        # Step 1: Search for the term and get UIDs
        uids = search_ncbi(args.search_term)

        if not uids:
            print(f"No UIDs found for search term: {args.search_term}")
        else:
            print(f"Found UIDs: {', '.join(uids)}")
            # Step 2: Fetch and print GenBank records for each UID
            for i, uid in enumerate(uids):
                if i > 0: # Add delay before the second request and onwards
                    time.sleep(1) # 1-second delay to respect NCBI rate limits
                print(f"\nFetching GenBank record for UID: {uid}")
                genbank_record = fetch_genbank_record(uid)
                if genbank_record:
                    # Print lines starting from "LOCUS"
                    print_record_from_locus(genbank_record)
                else:
                    print(f"Could not retrieve GenBank record for UID: {uid}")

    except requests.exceptions.RequestException as e:
        print(f"An HTTP error occurred: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
    finally:
        # Step 3: Final delay
        time.sleep(5)
        print("\nScript finished.")

def search_ncbi(term):
    """
    Searches NCBI Nucleotide database for a term and extracts UIDs.
    """
    headers = {'User-Agent': USER_AGENT}
    params = {
        'cmd': 'search',
        'term': term,
        'db': 'Nucleotide'
    }
    print(f"Submitting search request to NCBI for term: {term}")
    response = requests.post(SEARCH_URL, data=params, headers=headers)
    response.raise_for_status() # Raise an exception for HTTP errors

    # --- Debug: Print response text ---
    # print("--- Full Search Response Text ---")
    # print(response.text)
    # print("--- End of Full Search Response Text ---")
    # --- End Debug ---

    # Extract UIDs using regex, looking for <!-- docsumok ... -->
    # Example: <!-- docsumok 12345 67890 -->
    # The Perl script uses: /<!-- docsumok (.*) -->/
    # This regex captures everything between "docsumok " and " -->"
    # Then splits the captured string by space.
    
    uids = []
    # New regex to find UIDs in input tags like:
    # <input name="EntrezSystem2.PEntrez.Nuccore.Sequence_ResultsPanel.Sequence_RVDocSum.uid" sid="1" type="checkbox" id="UidCheckBox2194972897" value="2194972897" acc="NC_060941" />
    # Capture the value attribute.
    matches = re.findall(r'<input\s+name="[^"]*RVDocSum\.uid"[^>]*value="(\d+)"', response.text)
    uids.extend(matches) # re.findall will return a list of captured strings directly.
        
    return uids

def fetch_genbank_record(uid):
    """
    Fetches a GenBank record for a given UID using EFetch.
    """
    headers = {'User-Agent': USER_AGENT}
    params = {
        'db': 'nuccore',
        'id': uid,
        'rettype': 'gb', # GenBank format
        'retmode': 'text'
    }
    print(f"Fetching GenBank record from NCBI EFetch for UID: {uid}")
    response = requests.get(EFETCH_URL, params=params, headers=headers) # GET request for EFetch
    response.raise_for_status()
    return response.text.strip() # Should be plain text GenBank


def print_record_from_locus(record_text):
    """
    Prints a record starting from the line containing "LOCUS ".
    """
    lines = record_text.splitlines()
    start_printing = False
    for line in lines:
        if line.startswith("LOCUS "): # More precise check
            start_printing = True
        if start_printing:
            print(line)
    if not start_printing:
        print("LOCUS line not found in the record.")

if __name__ == "__main__":
    main()
