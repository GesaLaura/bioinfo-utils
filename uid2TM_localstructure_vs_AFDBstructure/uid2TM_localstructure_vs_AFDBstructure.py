import os
import pandas as pd
import requests
import subprocess
import numpy as np
from Bio.PDB import PDBParser

# --- CONFIGURATION ---

# 1. Path to your specific LOCAL Reference PDB (The model you are comparing AGAINST)
REFERENCE_PDB_PATH = "/home/gfreimann/Seafile/Freimann_PhD/1_Projects/Af1503_Henri/7TM/structural_similarity/mcx6817393_1_N_terminal_domain.pdb"  # <--- CHANGE THIS

# 2. List of UniProt IDs to fetch and compare
UNIPROT_IDS_TO_COMPARE = [
    "O29843",  
    "O28170",  
    "O28106",  
    "O28155",  
    "O28178",  
    "O28298",  
    "O28711",  
    "O28812",  
    "O28838",  
    "O28879",  
    "O29023", 
    "O29467",  
    "O29726"
]

# 3. Output Directory and Filename
WORKING_DIR = "tmalign_results"
OUTPUT_CSV_NAME = "one_vs_many_results.csv"

# 4. TM-align Command
TMALIGN_CMD = "TMalign"

# ---------------------

def fetch_pdb_url_via_api(uniprot_id):
    """Queries AlphaFold API for the correct PDB URL."""
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        r = requests.get(api_url)
        if r.status_code == 200:
            data = r.json()
            if data and isinstance(data, list):
                return data[0]['pdbUrl']
    except Exception:
        pass
    return None

def fetch_protein_structure(uniprot_id, output_path):
    """Downloads structure from AlphaFold."""
    # 1. Get dynamic URL
    pdb_url = fetch_pdb_url_via_api(uniprot_id)
    if not pdb_url:
        pdb_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"

    # 2. Download
    response = requests.get(pdb_url)
    if response.status_code == 200:
        with open(output_path, 'wb') as f:
            f.write(response.content)
        print(f"Downloaded: {uniprot_id}")
    else:
        raise ValueError(f"Error 404: Could not find PDB for {uniprot_id}")

def calculate_radius_of_gyration(pdb_path):
    """Calculates Rg using BioPython."""
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure('protein', pdb_path)
        atom_coords = np.array([atom.coord for atom in structure.get_atoms() if atom.name == 'CA']) 
        
        if len(atom_coords) == 0:
            atom_coords = np.array([atom.coord for atom in structure.get_atoms()])

        center_of_mass = np.mean(atom_coords, axis=0)
        rg = np.sqrt(np.mean(np.sum((atom_coords - center_of_mass)**2, axis=1)))
        return rg
    except Exception as e:
        print(f"Rg Error: {e}")
        return np.nan

def run_tmalign(mobile_pdb, target_pdb, output_prefix):
    """
    Aligns 'mobile_pdb' ONTO 'target_pdb'.
    Returns: Aligned_Length, RMSD, TM_Score (normalized by mobile length)
    """
    cmd = f"{TMALIGN_CMD} {mobile_pdb} {target_pdb} -o {output_prefix}"
    
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise Exception(f"TM-align failed: {result.stderr}")

    aligned_len = np.nan
    rmsd = np.nan
    tm_score = np.nan
    
    lines = result.stdout.split('\n')
    for line in lines:
        if line.startswith("Aligned length="):
            parts = line.split(',')
            for part in parts:
                if "Aligned length=" in part:
                    aligned_len = int(part.split('=')[1].strip())
                elif "RMSD=" in part:
                    rmsd = float(part.split('=')[1].strip())
        
        # Note: Chain_1 is the first argument (mobile_pdb)
        # We usually care about TM-score normalized by the length of the query (Chain_1)
        elif line.startswith("TM-score=") and "Chain_1" in line:
            tm_score = float(line.split('=')[1].split('(')[0].strip())

    return aligned_len, rmsd, tm_score

# --- MAIN EXECUTION ---

# 1. Setup Directories
if not os.path.exists(WORKING_DIR):
    os.makedirs(WORKING_DIR)

# 2. Check and Process Reference File
if not os.path.exists(REFERENCE_PDB_PATH):
    print(f"CRITICAL ERROR: Local reference file not found at: {REFERENCE_PDB_PATH}")
    exit()

print(f"Calculating Rg for Reference: {os.path.basename(REFERENCE_PDB_PATH)}...")
ref_rg = calculate_radius_of_gyration(REFERENCE_PDB_PATH)
ref_name = os.path.basename(REFERENCE_PDB_PATH)

# 3. Process List
results = []
print(f"Starting comparison of {len(UNIPROT_IDS_TO_COMPARE)} IDs against {ref_name}...")

for uid in UNIPROT_IDS_TO_COMPARE:
    uid = uid.strip() # Clean whitespace
    print(f"Processing {uid}...", end=" ")

    # Define paths
    mobile_pdb_path = os.path.join(WORKING_DIR, f"{uid}.pdb")
    alignment_prefix = os.path.join(WORKING_DIR, f"{uid}_aligned_to_ref")
    
    try:
        # A. Download (if not already there)
        if not os.path.exists(mobile_pdb_path):
            fetch_protein_structure(uid, mobile_pdb_path)
        
        # B. Calculate Rg for Mobile
        mobile_rg = calculate_radius_of_gyration(mobile_pdb_path)
        
        # C. Run TM-align
        # Argument 1: Mobile (The downloaded one)
        # Argument 2: Target (The local reference)
        aln_len, rmsd, tm = run_tmalign(mobile_pdb_path, REFERENCE_PDB_PATH, alignment_prefix)
        
        # D. Store Data
        results.append({
            "Query_ID": uid,
            "Reference_File": ref_name,
            "Rg_Query": mobile_rg,
            "Rg_Reference": ref_rg,
            "RMSD": rmsd,
            "TM-Score": tm,
            "Aligned_Length": aln_len
        })
        print(f"TM-Score: {tm}")

    except Exception as e:
        print(f"Failed: {e}")
        # Log failure
        results.append({
            "Query_ID": uid,
            "Reference_File": ref_name,
            "Rg_Query": np.nan,
            "Rg_Reference": ref_rg,
            "RMSD": np.nan,
            "TM-Score": np.nan,
            "Aligned_Length": np.nan
        })

# 4. Save to CSV
output_path = os.path.join(WORKING_DIR, OUTPUT_CSV_NAME)
df = pd.DataFrame(results)
df.to_csv(output_path, index=False)

print(f"\nProcessing complete! Results saved to: {output_path}")