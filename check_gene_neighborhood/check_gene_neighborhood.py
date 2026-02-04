"""
Genomic Distance Calculator for UniProt Protein Pairs
=====================================================

Author: Gesa Freimann
Date: 2026-02-04
Version: 2.0

Description
-----------
This script computes genomic proximity metrics between pairs of 
proteins identified by UniProt accessions or Entry Names. It is designed 
specifically for prokaryotic comparative genomics and operon analysis.

Workflow:
1. NCBI Data Acquisition: Downloads the RefSeq genomic annotation (GFF) for 
    a specified genome assembly.
2. GFF Indexing: Parses gene coordinates and builds a positional 
    index using multiple identifiers (NCBI GeneID, Locus Tag, Old Locus Tag, 
    and Gene Name) to maximize cross-database compatibility.
3. ID Mapping: Maps UniProt queries to genomic identifiers using a 
    two-step strategy:
   - Primary: Bulk mapping of UniProt IDs to NCBI GeneIDs via REST API.
   - Fallback: Direct UniProt metadata search for Ordered Locus Names (OLN) 
     or Primary Gene Names for any IDs that failed standard mapping.
4. Proximity Metrics Calculation: For each protein pair, the script computes:
   - Intergenic Base-pair Gap: Precise physical distance between gene boundaries.
   - Intervening Gene Count: Total number of genes located between the pair.
   - Strand-Specific Distance: Count of genes on the same strand (useful for 
     detecting potential operon co-membership).
5. Reporting: Generates a CSV file containing all genomic 
    coordinates and distances. Critically, every input pair is preserved in 
    the output; pairs that could not be mapped are assigned 'NaN' values 
    and a specific failure status for downstream filtering.


"""

# ============================================================================
# =========================== USER CONFIGURATION ==============================
# ============================================================================

# RefSeq genome accession (GCF_*)
REF_SEQ_ACC = "GCF_000091545.1" # thermus thermophilus H8

# Directory for downloaded GFF and temporary files
GFF_DIR = "."

# Input CSV containing UniProt ID pairs (two columns)
INPUT_CSV = "uids_pairs_test.csv"

# Output CSV (auto-derived from input name)
OUTPUT_CSV = f"{INPUT_CSV.replace('.csv', '')}_gene_distances.csv"

# Verbosity toggle
VERBOSE = True


# ============================================================================
# =============================== IMPORTS =====================================
# ============================================================================

import time
import requests
import os
import shutil
import zipfile
import csv
import re


# ============================================================================
# =========================== NCBI DATA ACCESS ================================
# ============================================================================

def download_ncbi_gff(ref_seq_acc, output_dir, verbose=True):
    """
    Download and extract the genomic GFF annotation for a RefSeq genome.

    Parameters
    ----------
    ref_seq_acc : str
        RefSeq genome accession (e.g. GCF_000091545.1)
    output_dir : str
        Directory to store the extracted GFF
    verbose : bool
        Print progress messages

    Returns
    -------
    str or None
        Path to extracted GFF file, or None on failure
    """
    os.makedirs(output_dir, exist_ok=True)

    url = (
        "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/"
        f"accession/{ref_seq_acc}/download"
        "?include_annotation_type=GENOME_GFF"
    )
    headers = {"Accept": "application/zip"}
    zip_path = os.path.join(output_dir, f"temp_{ref_seq_acc}.zip")

    if verbose:
        print(f"Downloading GFF for {ref_seq_acc}...")

    response = requests.get(url, headers=headers, timeout=600)
    if response.status_code != 200:
        print(f"Download failed (HTTP {response.status_code})")
        return None

    with open(zip_path, "wb") as f:
        f.write(response.content)

    with zipfile.ZipFile(zip_path, "r") as zip_ref:
        target_file = f"ncbi_dataset/data/{ref_seq_acc}/genomic.gff"
        zip_ref.extract(target_file, path=output_dir)

    final_gff = os.path.join(output_dir, f"{ref_seq_acc}.gff")
    source_path = os.path.join(output_dir, target_file)

    shutil.move(source_path, final_gff)
    shutil.rmtree(os.path.join(output_dir, "ncbi_dataset"))
    os.remove(zip_path)

    if verbose:
        print(f"GFF saved to: {final_gff}")

    return final_gff


# ============================================================================
# ============================ INPUT HANDLING =================================
# ============================================================================

def load_protein_pairs(input_file):
    """
    Load UniProt accession pairs from a CSV file.

    Parameters
    ----------
    input_file : str
        CSV file with two UniProt accessions per row

    Returns
    -------
    tuple
        (list of (u1, u2) pairs, set of all unique UniProt IDs)
    """
    pairs = []
    unique_ids = set()

    with open(input_file, "r") as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) < 2:
                continue
            u1, u2 = row[0].strip(), row[1].strip()
            pairs.append((u1, u2))
            unique_ids.update([u1, u2])

    return pairs, unique_ids


# ============================================================================
# ============================ UNIPROT MAPPING ================================
# ============================================================================

def chunk_list(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i : i + n]

import re

def fetch_locus_tags(unmapped_ids):
    """
    Fallback: Fetch Ordered Locus Names or Gene Names from UniProt.
    Supports both UniProt Accessions and Entry Names (IDs).
    """
    locus_mappings = {}
    if not unmapped_ids:
        return locus_mappings

    print(f"  Attempting metadata lookup for {len(unmapped_ids)} unmapped IDs...")
    
    # 1. Build query: Check if ID looks like an Accession or an Entry Name (ID)
    # Accessions: [OPQ][0-9][A-Z0-9]{3}[0-9] or [A-N,R-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}
    # Entry Names: usually Word_Word (e.g., RECA_THETH)
    query_parts = []
    for uid in unmapped_ids:
        if re.match(r"^[A-Z][0-9][A-Z0-9]{3}[0-9](-[0-9]+)?$", uid) or re.match(r"^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$", uid):
            query_parts.append(f"accession:{uid}")
        else:
            query_parts.append(f"id:{uid}")

    query = " OR ".join(query_parts)
    url = f"https://rest.uniprot.org/uniprotkb/search?query={query}&fields=id,accession,gene_oln,gene_primary&format=json"

    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        data = response.json()

        for entry in data.get("results", []):
            acc = entry.get("primaryAccession")
            entry_name = entry.get("uniProtkbId")
            
            # Determine which input ID matches this result
            # (Mapping back to the original input used in the query)
            match_key = None
            if acc in unmapped_ids:
                match_key = acc
            elif entry_name in unmapped_ids:
                match_key = entry_name

            if not match_key:
                continue

            genes = entry.get("genes", [])
            target_value = None

            for gene in genes:
                # Priority 1: Ordered Locus Name (OLN)
                olns = gene.get("orderedLocusNames", [])
                if olns:
                    target_value = olns[0].get("value")
                    break
                
                # Priority 2: Primary Gene Name
                primary = gene.get("geneName")
                if primary:
                    target_value = primary.get("value")
                    break

            if target_value:
                locus_mappings[match_key] = target_value

    except Exception as e:
        print(f"    Metadata lookup failed: {e}")

    if VERBOSE:
        print(f"    Successfully found alternative identifiers for {len(locus_mappings)} entries.")

    print(locus_mappings)

    return locus_mappings

def map_uniprot_to_ncbi(uniprot_ids, chunk_size=25):
    """
    Map UniProt accessions to NCBI GeneIDs, with fallback to Locus Tags.
    """
    all_mappings = {}
    
    uniprot_list = list(uniprot_ids)
    total_chunks = (len(uniprot_list) + chunk_size - 1) // chunk_size

    print(f"Mapping {len(uniprot_list)} IDs in {total_chunks} chunks...")

    for i, chunk in enumerate(chunk_list(uniprot_list, chunk_size), 1):
        print(f"  Processing chunk {i}/{total_chunks} ({len(chunk)} IDs)...")
        
        run_url = "https://rest.uniprot.org/idmapping/run"
        params = {
            "from": "UniProtKB_AC-ID",
            "to": "GeneID",
            "ids": ",".join(chunk)
        }

        try:
            # Submit Job
            response = requests.post(run_url, data=params, timeout=30)
            response.raise_for_status()
            job_id = response.json().get("jobId")

            # Poll Status
            status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
            while True:
                status_resp = requests.get(status_url, timeout=30)
                status_resp.raise_for_status()
                status = status_resp.json()
                
                if status.get("jobStatus") == "FINISHED" or "results" in status:
                    break
                if status.get("jobStatus") == "FAILED":
                    print(f"    Warning: Chunk {i} failed.")
                    break
                time.sleep(2)

            # Get Results
            results_url = f"https://rest.uniprot.org/idmapping/results/{job_id}"
            results_resp = requests.get(results_url, timeout=30)
            results_resp.raise_for_status()
            
            chunk_results = results_resp.json().get("results", [])
            for item in chunk_results:
                all_mappings[item["from"]] = item["to"]

        except Exception as e:
            print(f"    Error processing chunk {i}: {e}")
            continue # Move to next chunk even if one fails

    unmapped = set(uniprot_ids) - set(all_mappings.keys())
    
    if unmapped:
        locus_tags = fetch_locus_tags(list(unmapped))
        all_mappings.update(locus_tags) # Add locus tags to the main mapping dict

    # Final reporting
    still_unmapped = set(uniprot_ids) - set(all_mappings.keys())
    if still_unmapped:
        print(f"WARNING: {len(still_unmapped)} IDs total did not map to GeneID or Locus Tag.")

    return all_mappings

# ============================================================================
# ============================= GFF PARSING ==================================
# ============================================================================

def parse_gff_neighborhood(gff_path):
    """
    Parse gene features from a GFF file and build a positional index.

    Parameters
    ----------
    gff_path : str
        Path to RefSeq GFF file

    Returns
    -------
    tuple
        (sorted list of gene dicts, index mapping GeneID/locus_tag -> list index)
    """
    genes = []

    with open(gff_path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.rstrip().split("\t")
            if parts[2] != "gene":
                continue

            attributes = dict(
                item.split("=", 1)
                for item in parts[8].split(";")
                if "=" in item
            )

            genes.append({
                "seqid": parts[0],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "strand": parts[6],
                "locus_tag": attributes.get("locus_tag", ""),
                "old_locus_tag": attributes.get("old_locus_tag", ""),
                "gene_name": attributes.get("gene", ""),
                "dbxref": attributes.get("Dbxref", "")
            })

    genes.sort(key=lambda x: (x["seqid"], x["start"]))

    # Build the multi-key index
    index = {}
    for i, g in enumerate(genes):
        # Index by locus_tag
        if g["locus_tag"]:
            index[g["locus_tag"]] = i
        # Index by old locus_tag
        if g["old_locus_tag"]:
            index[g["old_locus_tag"]] = i
        # Index by Gene Name (e.g., 'recA')
        if g["gene_name"]:
            index[g["gene_name"]] = i
        # Index by NCBI GeneID
        for ref in g["dbxref"].split(","):
            if ref.startswith("GeneID:"):
                index[ref.split(":")[1]] = i

    return genes, index


# ============================================================================
# ========================== DISTANCE CALCULATION =============================
# ============================================================================

def calculate_distances(input_csv, output_csv, genes, gff_index, uniprot_to_gene):
    """
    Calculate distances for ALL pairs. Unmapped pairs will have NaN/None values.
    """
    results = []

    with open(input_csv, "r") as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) < 2:
                continue

            u1, u2 = row[0].strip(), row[1].strip()
            id_val1 = uniprot_to_gene.get(u1)
            id_val2 = uniprot_to_gene.get(u2)

            # Initialize all values as NaN or None
            # same_contig and same_strand are now NaN by default instead of False
            res = {
                "uniprot_id1": u1, 
                "mapped_id1": id_val1,
                "uniprot_id2": u2, 
                "mapped_id2": id_val2,
                "contig1": None, 
                "strand1": None,
                "contig2": None, 
                "strand2": None,
                "same_contig": float("nan"), 
                "same_strand": float("nan"),
                "gene_dist_all": float("nan"),
                "gene_dist_same_strand": float("nan"),
                "base_gap": float("nan"),
                "status": "fail"
            }

            # Attempt to find indices in GFF
            idx1 = gff_index.get(str(id_val1)) if id_val1 else None
            idx2 = gff_index.get(str(id_val2)) if id_val2 else None

            if idx1 is not None and idx2 is not None:
                g1, g2 = genes[idx1], genes[idx2]
                
                same_contig = g1["seqid"] == g2["seqid"]
                same_strand = g1["strand"] == g2["strand"]

                res.update({
                    "contig1": g1["seqid"], "strand1": g1["strand"],
                    "contig2": g2["seqid"], "strand2": g2["strand"],
                    "same_contig": same_contig,
                    "same_strand": same_strand,
                    "status": "success" if same_contig else "different_contigs"
                })

                if same_contig:
                    # Physical base-pair gap
                    if idx1 < idx2:
                        res["base_gap"] = max(0, g2['start'] - g1['end'] - 1)
                    else:
                        res["base_gap"] = max(0, g1['start'] - g2['end'] - 1)
                    
                    # Intervening gene count
                    res["gene_dist_all"] = abs(idx1 - idx2) - 1
                    
                    # Intervening same-strand count
                    if same_strand:
                        start_idx, end_idx = min(idx1, idx2), max(idx1, idx2)
                        intervening = genes[start_idx + 1 : end_idx]
                        res["gene_dist_same_strand"] = sum(1 for gx in intervening if gx['strand'] == g1['strand'])
            else:
                if VERBOSE:
                    reason = "Unmapped ID" if not id_val1 or not id_val2 else "Not in GFF"
                    print(f"Row ({u1}, {u2}): Distances set to NaN ({reason})")
                res["status"] = "mapping_failed"

            results.append(res)

    if results:
        fieldnames = [
            "uniprot_id1", "mapped_id1", "uniprot_id2", "mapped_id2",
            "contig1", "strand1", "contig2", "strand2",
            "same_contig", "same_strand", "gene_dist_all",
            "gene_dist_same_strand", "base_gap", "status"
        ]
        with open(output_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
        print(f"Processed {len(results)} pairs. Results saved to {output_csv}")


# ============================================================================
# =============================== MAIN =======================================
# ============================================================================

def main():
    gff_file = download_ncbi_gff(REF_SEQ_ACC, GFF_DIR, VERBOSE)
    if not gff_file:
        return

    pairs, unique_ids = load_protein_pairs(INPUT_CSV)
    uniprot_to_gene = map_uniprot_to_ncbi(list(unique_ids))

    genes, gff_index = parse_gff_neighborhood(gff_file)

    calculate_distances(
        INPUT_CSV,
        OUTPUT_CSV,
        genes,
        gff_index,
        uniprot_to_gene
    )


if __name__ == "__main__":
    main()
