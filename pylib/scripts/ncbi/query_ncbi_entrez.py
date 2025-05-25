#!/usr/bin/env python3

"""
Performs a search on NCBI Entrez and fetches the corresponding records.
Uses NCBI E-utilities (ESearch and EFetch).
"""

import argparse
import requests
import xml.etree.ElementTree as ET
import time
import sys
import os

# Add the parent directory of 'ncbi' to sys.path to allow importing utility_functions
# This assumes the script is in 'pylib/scripts/ncbi/'
current_dir = os.path.dirname(os.path.abspath(__file__))
pylib_scripts_dir = os.path.dirname(current_dir)
pylib_dir = os.path.dirname(pylib_scripts_dir)
if pylib_dir not in sys.path:
    sys.path.insert(0, pylib_dir)

# --- Global Request Settings ---
# These can be modified by command-line arguments
BASE_EUTILS_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
USER_AGENT = 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
REQUEST_COUNT = 0 # To track requests for delays

def get_request_delay(api_key_provided):
    """Returns the appropriate delay based on API key presence."""
    if api_key_provided:
        return 0.11  # 0.1 seconds + small buffer = ~10 requests/sec
    else:
        return 0.34  # 0.333 seconds + small buffer = ~3 requests/sec

def _make_request_with_retries(url, params, method='GET', data=None, max_retries=3, api_key_provided=False):
    """
    Makes an HTTP request with retries and appropriate delays.
    Returns the requests.Response object or None on failure.
    """
    global REQUEST_COUNT
    
    # Apply delay *before* making the request (except for the very first one)
    # NCBI guidelines: delay between requests.
    # This simple model applies it before every request.
    if REQUEST_COUNT > 0: # Don't sleep before the very first request
        delay = get_request_delay(api_key_provided)
        # print(f"Sleeping for {delay}s before next request...", file=sys.stderr)
        time.sleep(delay)

    for attempt in range(max_retries):
        try:
            REQUEST_COUNT += 1 # Increment after successful or failed attempt that consumes quota
            print(f"Making {method} request to {url} with params {params} (Attempt {attempt + 1})", file=sys.stderr)
            if method == 'GET':
                response = requests.get(url, params=params, headers={'User-Agent': USER_AGENT})
            elif method == 'POST':
                response = requests.post(url, data=params, headers={'User-Agent': USER_AGENT}) # params are in data for POST
            else:
                raise ValueError(f"Unsupported HTTP method: {method}")
            
            response.raise_for_status() # Raise HTTPError for bad responses (4XX or 5XX)
            return response
        except requests.exceptions.RequestException as e:
            print(f"Request failed (Attempt {attempt + 1}/{max_retries}): {e}", file=sys.stderr)
            if attempt < max_retries - 1:
                time.sleep(1 + attempt) # Exponential backoff (simple version)
            else:
                print("Max retries reached. Request failed.", file=sys.stderr)
                return None
    return None


def parse_eutils_error(xml_content):
    """
    Checks for <ePostResult><Error> or <eSearchResult><ErrorList> in NCBI XML response.
    Returns error message string if found, else None.
    """
    try:
        root = ET.fromstring(xml_content)
        # Common error paths in ESearch/EFetch XML
        error_message = None
        if root.find('ERROR') is not None: # Top-level ERROR
            error_message = root.find('ERROR').text
        elif root.find('ErrorList/PhraseNotFound') is not None: # eSearchResult
             error_message = "PhraseNotFound: " + root.find('ErrorList/PhraseNotFound').text
        elif root.find('ErrorList/FieldNotFound') is not None: # eSearchResult
             error_message = "FieldNotFound: " + root.find('ErrorList/FieldNotFound').text
        # For efetch, errors might be in <eFetchResult><ERROR>
        elif root.tag == 'eFetchResult' and root.find('ERROR') is not None:
            error_message = root.find('ERROR').text
        # For epost, errors might be in <ePostResult><Error>
        elif root.tag == 'ePostResult' and root.find('Error') is not None: # Note: case sensitive 'Error'
            error_message = root.find('Error').text

        if error_message:
            return f"NCBI E-utility Error: {error_message.strip()}"
            
    except ET.ParseError:
        # If it's not XML, it might be a plain text error from NCBI or a proxy
        if "error" in xml_content.lower() or "fail" in xml_content.lower():
            return f"Potential non-XML error from NCBI: {xml_content[:200]}" # Show snippet
    return None


def search_entrez(term, search_db, base_url, email, api_key, max_retries):
    """
    Performs an ESearch query on NCBI Entrez.
    Returns (uid_list, webenv, query_key) or (None, None, None) on error.
    """
    esearch_url = base_url + "esearch.fcgi"
    params = {
        'db': search_db,
        'term': term,
        'retmode': 'xml',
        'usehistory': 'y'
    }
    if email:
        params['email'] = email
    if api_key:
        params['api_key'] = api_key

    print(f"Searching Entrez ({search_db}) for term: '{term}'...", file=sys.stderr)
    response = _make_request_with_retries(esearch_url, params, method='GET', max_retries=max_retries, api_key_provided=bool(api_key))

    if not response:
        return None, None, None

    xml_content = response.text
    
    # Check for NCBI specific errors first
    ncbi_error = parse_eutils_error(xml_content)
    if ncbi_error:
        print(ncbi_error, file=sys.stderr)
        return None, None, None

    try:
        root = ET.fromstring(xml_content)
        uid_list = [uid_element.text for uid_element in root.findall('.//IdList/Id')]
        
        webenv_element = root.find('WebEnv')
        webenv = webenv_element.text if webenv_element is not None else None
        
        query_key_element = root.find('QueryKey')
        query_key = query_key_element.text if query_key_element is not None else None

        if not uid_list:
            print("No UIDs found for the search term.", file=sys.stderr)
            # Still return WebEnv and QueryKey as they might be present even with 0 results
            return [], webenv, query_key 
            
        print(f"Found {len(uid_list)} UIDs. WebEnv: {webenv}, QueryKey: {query_key}", file=sys.stderr)
        return uid_list, webenv, query_key

    except ET.ParseError as e:
        print(f"Error parsing ESearch XML response: {e}", file=sys.stderr)
        print(f"XML content snippet: {xml_content[:500]}", file=sys.stderr)
        return None, None, None


def fetch_entrez_records(fetch_db, rettype, retmode, base_url, email, api_key, max_retries,
                         uid_list=None, webenv=None, query_key=None):
    """
    Fetches records from NCBI Entrez using EFetch.
    Prints records directly to stdout.
    """
    efetch_url = base_url + "efetch.fcgi"
    params = {
        'db': fetch_db,
        'rettype': rettype,
        'retmode': retmode
    }
    if email:
        params['email'] = email
    if api_key:
        params['api_key'] = api_key

    method = 'GET' # Default method

    if webenv and query_key:
        params['WebEnv'] = webenv
        params['query_key'] = query_key
        print(f"Fetching records using WebEnv and QueryKey from {fetch_db}...", file=sys.stderr)
    elif uid_list:
        id_string = ','.join(uid_list)
        if len(uid_list) > 200 and len(id_string) > 1000 : # Heuristic for when to use POST
             print(f"Fetching {len(uid_list)} records using POST from {fetch_db} (UID list potentially long)...", file=sys.stderr)
             params['id'] = id_string
             method = 'POST'
        else:
            params['id'] = id_string
            print(f"Fetching {len(uid_list)} records using GET from {fetch_db}...", file=sys.stderr)
    else:
        print("Error: Neither UID list nor WebEnv/QueryKey provided for EFetch.", file=sys.stderr)
        return False

    response = _make_request_with_retries(efetch_url, params, method=method, max_retries=max_retries, api_key_provided=bool(api_key))

    if not response:
        return False

    content = response.text
    
    # Check for NCBI specific errors in XML/text (EFetch can return non-XML errors too)
    # For XML rettypes, parse_eutils_error works. For text (like gb, fasta), need careful checks.
    if retmode == 'xml' or rettype == 'xml': # If expected output is XML
        ncbi_error = parse_eutils_error(content)
        if ncbi_error:
            print(ncbi_error, file=sys.stderr)
            return False
    else: # For text formats, check for common error phrases if content is small
        if len(content) < 1000 and ("error" in content.lower() or "fail" in content.lower() or "not found" in content.lower()):
            # Avoid printing whole huge GenBank/FASTA file if it's actually data
            possible_error = parse_eutils_error(content) # Try parsing even if not XML
            if possible_error:
                 print(possible_error, file=sys.stderr)
                 return False
            elif "phrase not found" in content.lower() or "id list is empty" in content.lower():
                 print(f"NCBI E-utility Message: {content.strip()}", file=sys.stderr)
                 return False # Treat as error for fetching

    # Print results to stdout
    # Ensure final newline if not present, helps when redirecting output
    if content and not content.endswith('\n'):
        sys.stdout.write(content + '\n')
    else:
        sys.stdout.write(content)
    
    sys.stdout.flush()
    return True


def main():
    global BASE_EUTILS_URL # Allow modification by args

    parser = argparse.ArgumentParser(
        description="Search NCBI Entrez using ESearch and fetch records using EFetch.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--term", required=True, help="Search term (e.g., 'insulin human').")
    parser.add_argument("--search_db", default='protein', help="Database for ESearch (e.g., 'protein', 'nuccore').")
    parser.add_argument("--fetch_db", default='protein', help="Database for EFetch (default will match search_db if not specified, but explicit is better).")
    parser.add_argument("--rettype", default='gb', help="Retrieval type for EFetch (e.g., 'gb', 'fasta', 'xml').")
    parser.add_argument("--retmode", default='text', help="Retrieval mode for EFetch (e.g., 'text', 'xml').")
    parser.add_argument("--email", help="User email for NCBI E-utils (recommended).")
    parser.add_argument("--api_key", help="NCBI API key (recommended for higher request rates).")
    parser.add_argument("--max_retries", type=int, default=3, help="Max retries for HTTP requests.")
    parser.add_argument("--base_eutils_url", default=BASE_EUTILS_URL, help="Base URL for NCBI E-utilities.")
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
        
    args = parser.parse_args()

    # If fetch_db is not explicitly set by user, make it same as search_db for consistency
    # The previous default for fetch_db was 'nuccore', which might be confusing if search_db is 'protein'.
    # This ensures if user only gives --search_db protein, fetch_db also becomes protein by default.
    if args.fetch_db == 'protein' and parser.get_default('fetch_db') == 'protein' and args.search_db != 'protein':
        # This condition means fetch_db is at its default 'protein', but search_db is something else.
        # This is fine, user might want to search protein and fetch related nucleotides, or vice-versa.
        # However, if the user *didn't* specify fetch_db, it should probably match search_db.
        # Let's refine: if 'fetch_db' argument was NOT provided by user, then sync it.
        # Checking if 'fetch_db' is in sys.argv is a bit hacky. A better way is to see if it's different from default.
        # The ArgumentDefaultsHelpFormatter makes this tricky.
        # For now, let's assume if search_db is 'nuccore' and fetch_db default 'protein' is kept, it's intentional.
        # A common scenario: search 'nuccore', fetch 'nuccore'. search 'protein', fetch 'protein'.
        # If user only specifies --search_db=nuccore, and not --fetch_db, it's better if fetch_db defaults to nuccore.
        # The 'default' in add_argument for fetch_db should ideally be dynamic or user warned.
        # Simplest: if user provides search_db but not fetch_db, make fetch_db = search_db.
        # This requires checking if fetch_db was explicitly set.
        # A common way:
        if any(arg.startswith('--fetch_db') for arg in sys.argv):
            pass # User explicitly set fetch_db
        else:
            args.fetch_db = args.search_db
            print(f"Defaulting fetch_db to match search_db: '{args.fetch_db}'", file=sys.stderr)


    BASE_EUTILS_URL = args.base_eutils_url # Update global from arg

    uid_list, webenv, query_key = search_entrez(
        args.term, args.search_db, args.base_eutils_url,
        args.email, args.api_key, args.max_retries
    )

    if uid_list is None: # Indicates an error during search
        print("Exiting due to error in ESearch.", file=sys.stderr)
        sys.exit(1)
    
    if not uid_list: # No UIDs found, but search itself was successful
        print(f"No records found for term '{args.term}' in database '{args.search_db}'.", file=sys.stdout)
        sys.exit(0)

    print(f"Fetching records for {len(uid_list)} UIDs: {', '.join(uid_list[:10])}{'...' if len(uid_list) > 10 else ''}", file=sys.stderr)
    
    success = fetch_entrez_records(
        args.fetch_db, args.rettype, args.retmode, args.base_eutils_url,
        args.email, args.api_key, args.max_retries,
        uid_list=uid_list, webenv=webenv, query_key=query_key # Pass all to fetch
    )

    if success:
        print(f"Successfully fetched records.", file=sys.stderr)
    else:
        print(f"Failed to fetch all records.", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
