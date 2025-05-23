import math
import re

def format_blast_coordinates(hsp):
    """Converts Biopython HSP 0-based, exclusive-end coords to 1-based, inclusive-end."""
    q_start = hsp.query_start + 1
    q_end = hsp.query_end
    s_start = hsp.hit_start + 1
    s_end = hsp.hit_end
    return q_start, q_end, s_start, s_end

def get_hit_strand_str(hsp):
    """Returns 'Plus' or 'Minus' for the hit strand based on hsp.hit_strand."""
    return "Minus" if hsp.hit_strand == -1 else "Plus"

def format_evalue(evalue):
    """Formats an E-value into mantissa and exponent strings."""
    if evalue == 0.0:
        return "0.00", "0" # Ensure "0.00" for 0.0 evalue as per common BLAST output
    # Simplified approach based on previous worker implementation
    s = f"{evalue:.2e}" 
    parts = s.split('e')
    mantissa = parts[0]
    exponent = parts[1] if len(parts) > 1 else "0" # Should always have exponent part from .2e
    
    # Clean up exponent: remove leading '+' and then leading '0's unless it's just '0'
    if exponent.startswith('+'):
        exponent = exponent[1:]
    
    # Remove leading zeros from exponent, but keep '0' if that's what it is
    # e.g., -05 -> -5, +05 -> 5, 00 -> 0
    if len(exponent) > 1 and exponent.startswith('0'): # e.g. "05"
        exponent = exponent.lstrip('0')
        if not exponent: exponent = "0" 
    elif len(exponent) > 2 and exponent.startswith('-0'): # e.g. "-05"
        exponent = "-" + exponent[2:].lstrip('0')
        if exponent == "-": exponent = "-0" # Should not happen if original was e.g. -00

    if not exponent: # Handles case where exponent was "00" and became empty
        exponent = "0"
        
    return mantissa, exponent

def calculate_percent_metric(count, span, decimals=2):
    """Calculates percentage, typically for identity or positives."""
    if span == 0:
        return f"{0.0:.{decimals}f}"
    percent = (count / span) * 100.0
    return f"{percent:.{decimals}f}"

def parse_ncbi_header(header_string):
    """
    Attempts to parse NCBI-style headers (e.g., from hit.description or hit.id).
    Example: "gi|12345|gb|U00001.1|LOCUS1 Human herpesvirus 1"
    Returns: (accession, locus, title_remainder)
    More complex parsing might be needed for all NCBI header variants.
    This is a basic attempt.
    """
    accession = "N/A"
    locus = header_string # Default locus to the whole header if no parse
    title_remainder = header_string # Default title to whole header

    # Regex for typical NCBI headers like:
    # gi|12345|gb|U00001.1|LOCUS_NAME optional description
    # sp|P12345|PROT_NAME optional description
    # tr|Q12345|PROT_NAME optional description
    # Also handles cases like >P0DTC2 ORF1ab polyprotein
    # Or >some_id description here
    
    # Try a common NCBI pattern first
    ncbi_match = re.match(r"(?:gi\|\d+\|)?(?:gb|emb|dbj|ref|sp|tr|pir|prf|pdb)\|([\w\.]+)(?:\|([\w\.-]+))?(?:\s+(.*))?", header_string)
    # Group 1: Accession.Version (e.g., U00001.1 or P12345)
    # Group 2: Locus/Protein Name (optional) (e.g., LOCUS_NAME or PROT_NAME)
    # Group 3: Rest of description (optional)
    
    if ncbi_match:
        acc_ver = ncbi_match.group(1)
        accession = acc_ver.split('.')[0] # Remove version if present
        
        parsed_locus = ncbi_match.group(2)
        parsed_desc = ncbi_match.group(3)

        if parsed_locus:
            locus = parsed_locus
            title_remainder = parsed_desc if parsed_desc else parsed_locus # If no further desc, locus is title
        else: # No explicit locus/protein name after accession
            locus = acc_ver # Use accession.version as locus
            title_remainder = parsed_desc if parsed_desc else acc_ver
        
        # If title_remainder ended up being just the locus again, and there was no specific description part
        if title_remainder == locus and not parsed_desc :
            # This means the original header was likely just something like >sp|P12345|PROT_NAME
            # and title_remainder became PROT_NAME. A more general title is often missing.
             pass


    else: # Not a clear pipe-based NCBI header, try space separation
        parts = header_string.split(' ', 1)
        locus = parts[0] # First word is locus/ID
        accession = parts[0] # Assume first word might also be an accession
        if len(parts) > 1:
            title_remainder = parts[1]
        else:
            title_remainder = locus # If no space, locus is the whole title

    return accession, locus, title_remainder
