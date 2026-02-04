"""
Gene Neighborhood Extraction and UniProt Mapping Pipeline
=========================================================

Author: Gesa Freimann
Date: 2026-02-04
Version: 2.0

Description
-----------
This pipeline performs high-resolution genomic context analysis by extracting 
neighboring genes for specific protein queries and cross-referencing them 
across NCBI RefSeq and UniProtKB databases.

Workflow:
1. NCBI Data Acquisition: Downloads the RefSeq genomic annotation (GFF) for 
   the specified assembly.
2. GFF Parsing: Constructs a positional index using multiple keys 
   (NCBI GeneID, Locus Tag, Old Locus Tag, and Gene Name) to maximize 
   mapping compatibility.
3. Query Mapping: Uses a tiered strategy to locate query UniProt IDs 
   on the genome:
   - Primary: Bulk mapping of UniProt Accessions to NCBI GeneIDs.
   - Fallback: Direct UniProt metadata search for Ordered Locus Names (OLN) 
     or Primary Gene Names if standard mapping fails or is unindexed.
4. Neighborhood Extraction: Identifies X upstream and downstream neighbors 
   located on the same strand as the query gene.
5. Neighbor Mapping: Maps discovered neighbors back to 
   UniProtKB via RefSeq Protein IDs, GeneIDs, or direct cross-reference 
   searches.
6. Output: Generates a JSON file containing all metadata. 
   Critically, the output maintains 1:1 parity with the input; queries 
   that could not be mapped are included as empty entries for data integrity.

"""

# ============================================================================
# =========================== USER CONFIGURATION ==============================
# ============================================================================

# RefSeq genome accession (GCF_*)
REF_SEQ_ACC = "GCF_000091545.1" # thermus thermophilus H8

# Input CSV containing query UniProt accessions (one per line)
INPUT_CSV = "uids_test.csv"

# Output JSON file
OUTPUT_JSON = f"{INPUT_CSV.strip('.csv')}_neighbor_data.json"

# Number of upstream/downstream neighbors to collect per gene
X_NEIGHBORS = 3

# Output directory for downloaded and intermediate files
OUTPUT_DIR = "."

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
import json
import re
from collections import defaultdict


# ============================================================================
# =========================== NCBI DATA ACCESS ================================
# ============================================================================

def download_ncbi_gff(ref_seq_acc, output_dir, verbose=True):
    """
    Download and extract the genomic GFF file for a RefSeq genome accession.

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
    str
        Path to extracted GFF file
    """
    os.makedirs(output_dir, exist_ok=True)

    url = (
        f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/"
        f"accession/{ref_seq_acc}/download?include_annotation_type=GENOME_GFF"
    )
    headers = {"Accept": "application/zip"}
    zip_path = os.path.join(output_dir, f"temp_{ref_seq_acc}.zip")

    if verbose:
        print(f"Downloading GFF for {ref_seq_acc}...")

    r = requests.get(url, headers=headers, timeout=600)
    r.raise_for_status()

    with open(zip_path, "wb") as f:
        f.write(r.content)

    with zipfile.ZipFile(zip_path, "r") as zip_ref:
        target = f"ncbi_dataset/data/{ref_seq_acc}/genomic.gff"
        zip_ref.extract(target, path=output_dir)

    final_gff = os.path.join(output_dir, f"{ref_seq_acc}.gff")
    shutil.move(os.path.join(output_dir, target), final_gff)

    shutil.rmtree(os.path.join(output_dir, "ncbi_dataset"))
    os.remove(zip_path)

    return final_gff


# ============================================================================
# ============================ UNIPROT MAPPING ================================
# ============================================================================

def fetch_locus_tags_and_names(unmapped_ids):
    """
    Fallback: Fetch Ordered Locus Names or Gene Names from UniProt.
    Used for IDs that failed bulk mapping.
    """
    metadata_map = {}
    if not unmapped_ids:
        return metadata_map

    if VERBOSE:
        print(f"  Attempting metadata lookup for {len(unmapped_ids)} unmapped IDs...")
    
    # Identify if input is Accession or Entry Name for the query
    query_parts = []
    for uid in unmapped_ids:
        if re.match(r"^[A-Z][0-9][A-Z0-9]{3}[0-9](-[0-9]+)?$", uid) or \
           re.match(r"^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$", uid):
            query_parts.append(f"accession:{uid}")
        else:
            query_parts.append(f"id:{uid}")

    query = " OR ".join(query_parts)
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": query,
        "fields": "id,accession,gene_oln,gene_primary",
        "format": "json"
    }

    try:
        r = requests.get(url, params=params, timeout=30)
        r.raise_for_status()
        data = r.json()

        for entry in data.get("results", []):
            acc = entry.get("primaryAccession")
            entry_name = entry.get("uniProtkbId")
            
            # Map back to whichever input ID was used
            match_key = acc if acc in unmapped_ids else entry_name
            if not match_key: continue

            target_value = None
            for gene in entry.get("genes", []):
                # Priority 1: Ordered Locus Name
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
                metadata_map[match_key] = target_value
    except Exception as e:
        print(f"    Metadata fallback failed: {e}")

    return metadata_map

def uniprot_id_mapping(ids, from_db, to_db):
    """
    Perform bulk UniProt ID mapping using the UniProt REST API.

    Parameters
    ----------
    ids : iterable of str
        Source identifiers
    from_db : str
        Source database name (e.g. RefSeq_Protein, GeneID)
    to_db : str
        Target database name (e.g. UniProtKB)

    Returns
    -------
    dict
        Mapping: source_id -> list of target UniProt accessions
    """
    ids = sorted(set(str(i) for i in ids if i and str(i).lower() != "n/a"))
    if not ids:
        return {}

    print(f"Mapping {len(ids)} IDs from {from_db} to {to_db}...")

    run_url = "https://rest.uniprot.org/idmapping/run"
    params = {"from": from_db, "to": to_db, "ids": ",".join(ids)}

    r = requests.post(run_url, data=params)
    r.raise_for_status()
    job_id = r.json()["jobId"]

    status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
    while True:
        status = requests.get(status_url).json()
        if status.get("jobStatus") == "FINISHED" or "results" in status:
            break
        time.sleep(2)

    results_url = f"https://rest.uniprot.org/idmapping/results/{job_id}"
    results = requests.get(results_url).json()

    mapping = defaultdict(list)
    for item in results.get("results", []):
        src = item["from"]
        tgt = item["to"]
        if isinstance(tgt, dict):
            tgt = tgt.get("primaryAccession")
        mapping[src].append(tgt)

    return dict(mapping)


def uniprot_search_by_refseq(refseq):
    """
    Direct UniProtKB search using a RefSeq protein cross-reference.

    This is used as a fallback when bulk idmapping fails.

    Parameters
    ----------
    refseq : str
        RefSeq protein accession (e.g. WP_012345678.1)

    Returns
    -------
    list of str
        UniProtKB accessions
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"xref:RefSeq-{refseq}",
        "fields": "accession",
        "format": "json",
        "size": 20
    }
    r = requests.get(url, params=params, timeout=30)
    if r.status_code != 200:
        return []

    data = r.json()
    return [e["primaryAccession"] for e in data.get("results", [])]


# ============================================================================
# ============================= GFF PARSING ==================================
# ============================================================================

def parse_gff_metadata(gff_path):
    """
    Enhanced GFF Parser: Indexes GeneID, locus_tag, and gene name.
    """
    genes = []
    cds_info = {}

    with open(gff_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9: continue
            
            attr = dict(item.split("=", 1) for item in parts[8].split(";") if "=" in item)

            if parts[2] == "gene":
                gene_id = "N/A"
                for ref in attr.get("Dbxref", "").split(","):
                    if ref.startswith("GeneID:"):
                        gene_id = ref.split(":")[1]

                genes.append({
                    "seqid": parts[0],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "strand": parts[6],
                    "locus_tag": attr.get("locus_tag", "N/A"),
                    "old_locus_tag": attr.get("old_locus_tag", "N/A"),
                    "gene_name": attr.get("gene") or attr.get("Name", "N/A"),
                    "gene_id": gene_id,
                    "product": "N/A",
                    "protein_id": "N/A",
                })
            elif parts[2] == "CDS" and "locus_tag" in attr:
                cds_info.setdefault(attr["locus_tag"], {})
                if "product" in attr: cds_info[attr["locus_tag"]]["product"] = attr["product"]
                if "protein_id" in attr: cds_info[attr["locus_tag"]]["protein_id"] = attr["protein_id"]

    for g in genes:
        info = cds_info.get(g["locus_tag"], {})
        g["product"] = info.get("product", "N/A")
        g["protein_id"] = info.get("protein_id", "N/A")

    genes.sort(key=lambda x: (x["seqid"], x["start"]))
    
    # Robust Multi-Key Index
    index = {}
    for i, g in enumerate(genes):
        if g["gene_id"] != "N/A": index[str(g["gene_id"])] = i
        if g["locus_tag"] != "N/A": index[g["locus_tag"]] = i
        if g["old_locus_tag"] != "N/A": index[g["old_locus_tag"]] = i
        if g["gene_name"] != "N/A": index[g["gene_name"]] = i

    return genes, index


# ============================================================================
# =============================== MAIN =======================================
# ============================================================================

def main():
    # 1. Load input IDs
    if not os.path.exists(INPUT_CSV):
        print(f"Error: Input file {INPUT_CSV} not found.")
        return

    with open(INPUT_CSV) as f:
        uids = [row[0].strip() for row in csv.reader(f) if row]

    # Initialize EVERY input ID with an empty list to ensure they appear in JSON
    final_output = {uid: [] for uid in uids}

    # 2. Download and Parse GFF
    try:
        gff_file = download_ncbi_gff(REF_SEQ_ACC, OUTPUT_DIR, VERBOSE)
    except Exception as e:
        print(f"Failed to download GFF: {e}")
        return
        
    genes, gff_index = parse_gff_metadata(gff_file)

    # 3. Phase 1: Bulk Mapping (UniProt -> GeneID)
    query_to_gene = uniprot_id_mapping(uids, "UniProtKB_AC-ID", "GeneID")

    # 4. Phase 2: Metadata Fallback for IDs missing from GFF
    missing_from_gff = []
    for uid in uids:
        mapped_gids = query_to_gene.get(uid, [])
        # Check if mapping failed OR if none of the mapped GeneIDs are in our GFF index
        if not mapped_gids or not any(str(gid) in gff_index for gid in mapped_gids):
            missing_from_gff.append(uid)
    
    if missing_from_gff:
        metadata_backups = fetch_locus_tags_and_names(missing_from_gff)
        for uid, backup_val in metadata_backups.items():
            # Add the backup identifier (OLN or Gene Name) to the searchable list
            query_to_gene.setdefault(uid, []).append(backup_val)

    # 5. Neighborhood Collection
    neighborhood = {}
    all_proteins = set()
    all_geneids = set()

    for uid in uids:
        # Get the first identifier that successfully hits the GFF index
        gene_id = next((gid for gid in query_to_gene.get(uid, []) if str(gid) in gff_index), None)
        
        if not gene_id:
            if VERBOSE:
                print(f"Warning: No valid genomic location found for {uid}")
            continue

        center = gff_index[str(gene_id)]
        target = genes[center]
        neighbors = []

        # Helper to collect neighbors based on strand and distance
        def collect(step):
            found = 0
            idx = center + step
            while 0 <= idx < len(genes) and found < X_NEIGHBORS:
                g = genes[idx]
                if g["seqid"] != target["seqid"]:
                    break
                # Only collect neighbors on the same strand
                if g["strand"] == target["strand"]:
                    found += 1
                    n = g.copy()
                    n["relative_position"] = found * (1 if step > 0 else -1)
                    neighbors.append(n)
                    if n["protein_id"] != "N/A":
                        all_proteins.add(n["protein_id"])
                    if n["gene_id"] != "N/A":
                        all_geneids.add(n["gene_id"])
                idx += step

        collect(-1) # Upstream
        collect(1)  # Downstream
        neighborhood[uid] = neighbors

    # 6. Map the Neighboring Genes back to UniProt
    protein_to_uniprot = uniprot_id_mapping(list(all_proteins), "RefSeq_Protein", "UniProtKB")
    gene_to_uniprot = uniprot_id_mapping(list(all_geneids), "GeneID", "UniProtKB")
    refseq_search_cache = {}

    # 7. Compile final data
    for uid, neighbors in neighborhood.items():
        processed_neighbors = []
        for n in neighbors:
            # Tiered mapping for the neighbor
            if n["protein_id"] in protein_to_uniprot:
                n["uniprot_ids"] = protein_to_uniprot[n["protein_id"]]
                n["uniprot_mapping_source"] = "RefSeq_Protein"
            elif n["gene_id"] in gene_to_uniprot:
                n["uniprot_ids"] = gene_to_uniprot[n["gene_id"]]
                n["uniprot_mapping_source"] = "GeneID"
            elif n["protein_id"] != "N/A":
                if n["protein_id"] not in refseq_search_cache:
                    refseq_search_cache[n["protein_id"]] = uniprot_search_by_refseq(n["protein_id"])
                n["uniprot_ids"] = refseq_search_cache[n["protein_id"]]
                n["uniprot_mapping_source"] = "UniProt_search_RefSeq" if n["uniprot_ids"] else "None"
            else:
                n["uniprot_ids"] = []
                n["uniprot_mapping_source"] = "None"
            
            processed_neighbors.append(n)
        
        final_output[uid] = processed_neighbors

    # 8. Save results
    with open(OUTPUT_JSON, "w") as f:
        json.dump(final_output, f, indent=4)

    if VERBOSE:
        success_count = sum(1 for v in final_output.values() if v)
        print(f"Analysis complete. Found neighborhoods for {success_count}/{len(uids)} IDs.")
        print(f"JSON saved to {OUTPUT_JSON}")

if __name__ == "__main__":
    main()


if __name__ == "__main__":
    main()
