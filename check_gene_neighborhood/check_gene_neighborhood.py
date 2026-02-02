"""
Genomic Distance Calculator for UniProt Protein Pairs
=====================================================

Author: Gesa Freimann
Date: 2026-02-02

Description
-----------
This script computes genomic proximity metrics between pairs of proteins
identified by UniProt accessions. The workflow is:

1. Download a RefSeq genome annotation (GFF) from NCBI
2. Parse gene coordinates and identifiers from the GFF
3. Map UniProt accessions to NCBI GeneIDs via UniProt REST API
4. For each UniProt pair, compute:
   - Intergenic base-pair gap
   - Number of intervening genes (any strand)
   - Number of intervening genes on the same strand
5. Output results to a CSV file

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

def map_uniprot_to_ncbi(uniprot_ids):
    """
    Map UniProt accessions to NCBI GeneIDs using UniProt ID Mapping API.

    Parameters
    ----------
    uniprot_ids : list of str
        UniProt accessions

    Returns
    -------
    dict
        Mapping: UniProt_ID -> GeneID
    """
    print(f"Mapping {len(uniprot_ids)} UniProt IDs to GeneID...")

    run_url = "https://rest.uniprot.org/idmapping/run"
    params = {
        "from": "UniProtKB_AC-ID",
        "to": "GeneID",
        "ids": ",".join(uniprot_ids)
    }

    response = requests.post(run_url, data=params)
    response.raise_for_status()
    job_id = response.json()["jobId"]

    status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
    while True:
        status = requests.get(status_url).json()
        if "results" in status or status.get("jobStatus") == "FINISHED":
            break
        time.sleep(2)

    results_url = f"https://rest.uniprot.org/idmapping/results/{job_id}"
    data = requests.get(results_url).json()

    return {item["from"]: item["to"] for item in data.get("results", [])}


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
                "dbxref": attributes.get("Dbxref", "")
            })

    genes.sort(key=lambda x: (x["seqid"], x["start"]))

    index = {}
    for i, g in enumerate(genes):
        if g["locus_tag"]:
            index[g["locus_tag"]] = i
        for ref in g["dbxref"].split(","):
            if ref.startswith("GeneID:"):
                index[ref.split(":")[1]] = i

    return genes, index


# ============================================================================
# ========================== DISTANCE CALCULATION =============================
# ============================================================================

def calculate_distances(input_csv, output_csv, genes, gff_index, uniprot_to_gene):
    """
    Calculate genomic distance metrics for UniProt protein pairs.

    Metrics include:
    - Intergenic base-pair gap
    - Number of intervening genes (any strand)
    - Number of intervening genes on same strand

    Cross-contig comparisons return NaN distances.
    """
    results = []

    with open(input_csv, "r") as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) < 2:
                continue

            u1, u2 = row[0].strip(), row[1].strip()
            gid1 = uniprot_to_gene.get(u1)
            gid2 = uniprot_to_gene.get(u2)

            if gid1 not in gff_index or gid2 not in gff_index:
                print(f"Warning: mapping failed for {u1} or {u2}")
                continue

            idx1, idx2 = gff_index[gid1], gff_index[gid2]
            g1, g2 = genes[idx1], genes[idx2]

            same_contig = g1["seqid"] == g2["seqid"]
            same_strand = g1["strand"] == g2["strand"]

            gene_dist_all = float("nan")
            gene_dist_same_strand = float("nan")
            base_gap = float("nan")

            if same_contig:
                # Physical distance: Intergenic gap (bases between features)
                if idx1 < idx2:
                    base_gap = max(0, g2['start'] - g1['end'] - 1)
                else:
                    base_gap = max(0, g1['start'] - g2['end'] - 1)
                
                # Neighborhood distance: count of any genes between the pair
                gene_dist_all = abs(idx1 - idx2) - 1
                
                # Functional distance: count of genes on the same strand (operon check)
                if same_strand:
                    start_idx, end_idx = min(idx1, idx2), max(idx1, idx2)
                    intervening = genes[start_idx + 1 : end_idx]
                    gene_dist_same_strand = sum(1 for gx in intervening if gx['strand'] == g1['strand'])

            results.append({
                # Protein / gene identifiers
                "uniprot_id1": u1,
                "gene_id1": gid1,
                "uniprot_id2": u2,
                "gene_id2": gid2,

                # Gene 1 coordinates
                "contig1": g1["seqid"],
                "start1": g1["start"],
                "end1": g1["end"],
                "size1": g1["end"] - g1["start"] + 1,
                "strand1": g1["strand"],

                # Gene 2 coordinates
                "contig2": g2["seqid"],
                "start2": g2["start"],
                "end2": g2["end"],
                "size2": g2["end"] - g2["start"] + 1,
                "strand2": g2["strand"],

                # Structural relationships
                "same_contig": same_contig,
                "same_strand": same_strand,

                # Distance metrics
                "gene_dist_all": gene_dist_all,
                "gene_dist_same_strand": gene_dist_same_strand,
                "base_gap": base_gap
            })


    if results:
        with open(output_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)

        print(f"Results saved to {output_csv}")


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
