"""
Gene Neighborhood Extraction and UniProt Mapping Pipeline
=========================================================

Author: Gesa Freimann
Date: 2026-02-02

Description
-----------
This script:
1. Downloads a RefSeq genome annotation (GFF) from NCBI
2. Parses gene coordinates, products, GeneIDs, and RefSeq protein IDs
3. Identifies genomic neighbors for a list of query UniProt IDs
4. Maps neighboring genes to UniProt accessions using a tiered strategy:
   - RefSeq protein → UniProtKB (bulk mapping)
   - GeneID → UniProtKB (bulk mapping)
   - RefSeq protein → UniProtKB (direct UniProt search fallback)
5. Outputs a structured JSON file

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
    Parse a RefSeq GFF file to extract gene-level metadata.

    Extracted fields:
    - GeneID
    - locus_tag
    - product description
    - RefSeq protein accession

    Returns
    -------
    genes : list of dict
        Sorted list of gene records
    index : dict
        Mapping GeneID -> index in genes list
    """
    genes = []
    cds_info = {}

    with open(gff_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.rstrip().split("\t")
            if len(parts) < 9:
                continue

            attr = dict(
                item.split("=", 1)
                for item in parts[8].split(";")
                if "=" in item
            )

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
                    "gene_id": gene_id,
                    "product": "N/A",
                    "protein_id": "N/A",
                })

            elif parts[2] == "CDS" and "locus_tag" in attr:
                cds_info.setdefault(attr["locus_tag"], {})
                if "product" in attr:
                    cds_info[attr["locus_tag"]]["product"] = attr["product"]
                if "protein_id" in attr:
                    cds_info[attr["locus_tag"]]["protein_id"] = attr["protein_id"]

    for g in genes:
        info = cds_info.get(g["locus_tag"], {})
        g["product"] = info.get("product", "N/A")
        g["protein_id"] = info.get("protein_id", "N/A")

    genes.sort(key=lambda x: (x["seqid"], x["start"]))
    index = {g["gene_id"]: i for i, g in enumerate(genes) if g["gene_id"] != "N/A"}

    return genes, index


# ============================================================================
# =============================== MAIN =======================================
# ============================================================================

def main():
    with open(INPUT_CSV) as f:
        uids = [row[0].strip() for row in csv.reader(f) if row]

    gff_file = download_ncbi_gff(REF_SEQ_ACC, OUTPUT_DIR, VERBOSE)
    genes, gff_index = parse_gff_metadata(gff_file)

    query_to_gene = uniprot_id_mapping(
        uids, "UniProtKB_AC-ID", "GeneID"
    )

    neighborhood = {}
    all_proteins = set()
    all_geneids = set()

    for uid in uids:
        gene_id = next(iter(query_to_gene.get(uid, [])), None)
        if gene_id not in gff_index:
            continue

        center = gff_index[gene_id]
        target = genes[center]
        neighbors = []

        def collect(step):
            found = 0
            idx = center + step
            while 0 <= idx < len(genes) and found < X_NEIGHBORS:
                g = genes[idx]
                if g["seqid"] != target["seqid"]:
                    break
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

        collect(-1)
        collect(1)
        neighborhood[uid] = neighbors

    protein_to_uniprot = uniprot_id_mapping(
        all_proteins, "RefSeq_Protein", "UniProtKB"
    )
    gene_to_uniprot = uniprot_id_mapping(
        all_geneids, "GeneID", "UniProtKB"
    )

    refseq_search_cache = {}

    final = {}
    for uid, neighbors in neighborhood.items():
        final[uid] = []
        for n in neighbors:
            if n["protein_id"] in protein_to_uniprot:
                n["uniprot_ids"] = protein_to_uniprot[n["protein_id"]]
                n["uniprot_mapping_source"] = "RefSeq_Protein"
            elif n["gene_id"] in gene_to_uniprot:
                n["uniprot_ids"] = gene_to_uniprot[n["gene_id"]]
                n["uniprot_mapping_source"] = "GeneID"
            elif n["protein_id"] != "N/A":
                refseq_search_cache.setdefault(
                    n["protein_id"],
                    uniprot_search_by_refseq(n["protein_id"])
                )
                n["uniprot_ids"] = refseq_search_cache[n["protein_id"]]
                n["uniprot_mapping_source"] = (
                    "UniProt_search_RefSeq"
                    if n["uniprot_ids"] else "None"
                )
            else:
                n["uniprot_ids"] = []
                n["uniprot_mapping_source"] = "None"

            final[uid].append(n)

    with open(OUTPUT_JSON, "w") as f:
        json.dump(final, f, indent=4)

    print(f"Analysis complete. JSON saved to {OUTPUT_JSON}")


if __name__ == "__main__":
    main()
