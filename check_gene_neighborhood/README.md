# Genomic Neighborhood and Distance Analysis for UniProt Proteins

This repository contains two complementary Python scripts for analyzing **genomic context and proximity** of proteins identified by **UniProt IDs**, using **RefSeq genome annotations** and the **UniProt REST API**.

Together, these scripts enable:

* Extraction of genomic neighborhoods around query proteins
* Mapping between UniProt, GeneID, and RefSeq identifiers
* Calculation of basepair and gene-count distances between protein pairs

---

## Scripts Overview

### 1. Gene Neighborhood Extraction and UniProt Mapping

**Purpose:**
Identify upstream and downstream genomic neighbors of query proteins and map all neighboring genes back to UniProt.

**Key steps:**

1. Download a RefSeq genome annotation (GFF) from NCBI
2. Parse gene coordinates, GeneIDs, locus tags, products, and RefSeq protein IDs
3. Locate genomic neighbors for each query UniProt ID
4. Map neighbors to UniProt using a tiered strategy:

   * RefSeq protein → UniProtKB (bulk mapping)
   * GeneID → UniProtKB (bulk mapping)
   * RefSeq protein → UniProtKB (direct UniProt search fallback)
5. Export results as a structured JSON file

**Output:**
A JSON file containing, for each query UniProt ID, a list of neighboring genes with:

* genomic coordinates
* strand
* relative position
* RefSeq protein IDs
* mapped UniProt accessions (with mapping source)

---

### 2. Genomic Distance Calculator for UniProt Protein Pairs

**Purpose:**
Compute genomic proximity between **pairs of proteins** defined by UniProt IDs.

**Key steps:**

1. Download a RefSeq genome annotation (GFF) from NCBI
2. Parse gene coordinates and identifiers
3. Map UniProt accessions to NCBI GeneIDs
4. For each UniProt pair, calculate:

   * Intergenic base-pair gap
   * Number of intervening genes (any strand)
   * Number of intervening genes on the same strand
5. Export results to a CSV file

**Output:**
A CSV file with one row per UniProt pair, including:

* full gene coordinates for both proteins
* contig and strand information
* structural relationship flags
* distance metrics suitable for downstream statistical analysis

---

## Requirements

* Python 3.12
* Internet connection (for NCBI and UniProt APIs)
* Python packages:

  * `requests`

All other dependencies are from the Python standard library.

---

## Usage

Each script is configured via a **USER CONFIGURATION** section at the top of the file.

### 1. Configure script variables

Open the script you want to run and modify the variables under:

```python
# =========================== USER CONFIGURATION ==============================
```

Typical parameters include:

* `REF_SEQ_ACC` – RefSeq genome accession (e.g. `GCF_000091545.1`)
* `INPUT_CSV` – input file with UniProt IDs (single IDs or pairs)
* `OUTPUT_JSON` / `OUTPUT_CSV` – output file name
* `X_NEIGHBORS` – number of upstream/downstream neighbors (neighborhood script)
* `VERBOSE` – toggle status messages

---

### 2. Run the scripts

From the command line:

**Gene neighborhood extraction**

```bash
python3 get_nearest_genes.py
```

**Genomic distance calculation**

```bash
python3 check_gene_neighborhood.py
```

## Author

**Gesa Freimann**
2026-02-02
