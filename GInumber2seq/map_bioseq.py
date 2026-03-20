import pandas as pd
from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML
import time
import argparse
import sys

def parse_fasta(fasta_data):
    """Parses FASTA string into (header, accession, sequence)."""
    lines = fasta_data.splitlines()
    if not lines: return None
    header = lines[0]
    sequence = "".join(lines[1:])
    
    accession = "N/A"
    if "|" in header:
        parts = header.split("|")
        for p in parts:
            if "_" in p and len(p) > 4: 
                accession = p
                break
    else:
        accession = header.split()[0].replace(">", "")
    return header, accession, sequence

def fetch_single_gi(gi, email, db="protein", retries=3):
    """Standard NCBI E-fetch with retry logic."""
    Entrez.email = email
    attempt = 0
    while attempt < retries:
        try:
            handle = Entrez.efetch(db=db, id=str(gi), rettype="fasta", retmode="text")
            fasta_data = handle.read().strip()
            handle.close()
            if fasta_data:
                return parse_fasta(fasta_data)
        except Exception:
            attempt += 1
            time.sleep(1)
    return None

def blast_backup(fragment_seq, email):
    """BLAST fragment and verify 100% identity match in the subject."""
    Entrez.email = email
    try:
        result_handle = NCBIWWW.qblast("blastp", "nr", fragment_seq, hitlist_size=3)
        blast_record = NCBIXML.read(result_handle)
        result_handle.close()

        for alignment in blast_record.alignments:
            full_record = fetch_single_gi(alignment.accession, email)
            if full_record:
                _, _, full_seq = full_record
                if fragment_seq.upper() in full_seq.upper():
                    return full_record
    except Exception as e:
        print(f"      BLAST failed: {e}")
    return None

def main():
    parser = argparse.ArgumentParser(description="Map GI numbers to full protein sequences using NCBI Entrez and optional BLAST backup.")
    parser.add_argument("-i", "--input", required=True, help="Path to input CSV file")
    parser.add_argument("-o", "--output", required=True, help="Path to save output CSV file")
    parser.add_argument("-e", "--email", required=True, help="Your email (required by NCBI)")
    parser.add_argument("--gi_col", default="gi", help="Name of the column containing GI/Accession numbers")
    parser.add_argument("--seq_col", help="Optional: Name of column containing fragments for BLAST backup")
    
    args = parser.parse_args()

    print(f"Reading {args.input}...")
    try:
        df = pd.read_csv(args.input)
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)

    results = []
    total = len(df)

    for idx, row in df.iterrows():
        gi = row.get(args.gi_col)
        fragment = row.get(args.seq_col) if args.seq_col else None
        
        print(f"[{idx+1}/{total}] Querying ID: {gi}")
        
        # Step 1: GI Lookup
        res = None
        if pd.notna(gi):
            search_id = str(int(gi)) if isinstance(gi, (float, int)) else str(gi)
            res = fetch_single_gi(search_id, args.email)

        # Step 2: BLAST Backup (Optional)
        if res is None and fragment and pd.notna(fragment):
            print(f"    ID {gi} not found. Running BLAST backup...")
            res = blast_backup(str(fragment), args.email)
            if res: print(f"    Found via BLAST: {res[1]}")

        results.append(res if res else (None, None, None))
        time.sleep(0.4) # NCBI rate limit safety

    df['fasta_header'], df['accession'], df['full_sequence'] = zip(*results)
    df.to_csv(args.output, index=False)
    print(f"Success! Saved to {args.output}")

if __name__ == "__main__":
    main()
