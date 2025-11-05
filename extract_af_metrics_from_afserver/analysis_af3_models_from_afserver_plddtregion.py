"""
Extract metrics from AF3 outputs from AlphaFold server, 
using CIF files for region-specific pLDDT and JSON files for PAE/ranking scores.

author: Gesa Freimann, 2025
"""

import os
import json
import pandas as pd
import numpy as np
from Bio.PDB import MMCIFParser


def extract_af3_metrics(folder, regions_of_interest=None):
    """
    Extract metrics from AF3 outputs, using CIF files for region-specific pLDDT
    and JSON files for PAE/ranking scores.
    """
    metrics_list = []
    
    # Find all *_summary_confidences_*.json files
    summary_files = [f for f in os.listdir(folder) 
                     if f.endswith("_summary_confidences_0.json") or "_summary_confidences_" in f]
    
    for summary_file in summary_files:
        parts = summary_file.split("_summary_confidences_")
        if len(parts) != 2:
            print(f"Warning: Skipping file with unexpected name format: {summary_file}")
            continue
        
        base_name = parts[0] # The part before _summary...
        index = parts[1].replace(".json", "") # The number (e.g., '0')
        
        # Construct paths for all 3 relevant files
        summary_path = os.path.join(folder, summary_file)
        full_data_file = f"{base_name}_full_data_{index}.json"
        full_data_path = os.path.join(folder, full_data_file)
        cif_file = f"{base_name}_model_{index}.cif" 
        cif_path = os.path.join(folder, cif_file)

        model_metrics = {
            "folder": folder,
            "model_index": index,
            "iptm": None,
            "ptm": None,
            "ranking_score": None,
            "mean_plddt_global": None,
            "mean_pae": None
        }
        
        # Add keys for regional pLDDTs (initialized to nan)
        if regions_of_interest:
            for region in regions_of_interest:
                model_metrics[f"mean_plddt_{region['name']}"] = np.nan

        # 1. Parse summary_confidences_*.json (for ranking scores)
        if os.path.exists(summary_path):
            try:
                with open(summary_path) as f:
                    summary_data = json.load(f)
                model_metrics["iptm"] = summary_data.get("iptm")
                model_metrics["ptm"] = summary_data.get("ptm")
                model_metrics["ranking_score"] = summary_data.get("ranking_score")
            except json.JSONDecodeError:
                print(f"Error: Could not decode JSON from {summary_path}")

        # 2. Parse full_data_*.json (for PAE only)
        if os.path.exists(full_data_path):
            try:
                with open(full_data_path) as f:
                    full_data = json.load(f)
                if "pae" in full_data:
                    pae_array = np.array(full_data["pae"])
                    model_metrics["mean_pae"] = float(np.mean(pae_array))
            except json.JSONDecodeError:
                print(f"Error: Could not decode JSON from {full_data_path}")
        
        # 3. Parse ranked_*.cif (for all pLDDT calculations)
        if os.path.exists(cif_path):
            parser = MMCIFParser(QUIET=True) # QUIET=True suppresses warnings
            try:
                structure = parser.get_structure("model", cif_path)
                model = structure[0] # Get the first model
                
                all_plddts = []
                # Create a dict to hold lists of plddts per region
                region_plddts_lists = {region['name']: [] for region in regions_of_interest}

                for chain in model:
                    chain_id = chain.id
                    for residue in chain:
                        # Skip HETATMs, water, etc. (non-standard residues)
                        if residue.id[0] != ' ':
                            continue
                        
                        res_id = residue.id[1] # The residue number
                        
                        is_in_region = False
                        target_region_name = None
                        
                        # Check if this residue belongs to a region of interest
                        if regions_of_interest:
                            for region in regions_of_interest:
                                if (chain_id == region['chain'] and
                                    res_id >= region['start'] and
                                    res_id <= region['end']):
                                    
                                    is_in_region = True
                                    target_region_name = region['name']
                                    break
                        
                        # Collect pLDDTs from all atoms in this residue
                        for atom in residue:
                            plddt = atom.bfactor # pLDDT is stored in the B-factor column
                            all_plddts.append(plddt)
                            
                            if is_in_region:
                                region_plddts_lists[target_region_name].append(plddt)
                
                # Now calculate the means from the collected lists
                if all_plddts:
                    model_metrics["mean_plddt_global"] = np.mean(all_plddts)
                
                if regions_of_interest:
                    for region_name, plddts_list in region_plddts_lists.items():
                        if plddts_list:
                            model_metrics[f"mean_plddt_{region_name}"] = np.mean(plddts_list)
                        else:
                            # This handles cases where the region was defined but no atoms were found
                            print(f"Warning: Region '{region_name}' had 0 atoms in {cif_path}")
                            model_metrics[f"mean_plddt_{region_name}"] = np.nan

            except Exception as e:
                print(f"Error parsing CIF file {cif_path}: {e}")
        else:
            print(f"Warning: Could not find matching CIF file: {cif_path}. Skipping pLDDT calculations.")
            
        metrics_list.append(model_metrics)
    
    return metrics_list


def main(base_dir=".", REGIONS_OF_INTEREST=None, output_path="af3_metrics_summary_plddt_per_region.csv"):
    
    all_metrics = []
    for root, dirs, files in os.walk(base_dir):
        if any("_summary_confidences_" in f for f in files):
            all_metrics.extend(extract_af3_metrics(root, REGIONS_OF_INTEREST))
    
    if not all_metrics:
        print("No metrics found. Check your base directory and file names.")
        return

    df = pd.DataFrame(all_metrics)

    # Calculate mean of specified regions
    if REGIONS_OF_INTEREST:
        # Get the list of column names for regional pLDDTs
        region_plddt_cols = [f"mean_plddt_{region['name']}" for region in REGIONS_OF_INTEREST]
        
        # Filter for columns that actually exist in the DataFrame
        existing_region_cols = [col for col in region_plddt_cols if col in df.columns]
        
        if existing_region_cols:
            print(f"Calculating 'mean_plddt_regions' from: {existing_region_cols}")
            df['mean_plddt_regions'] = df[existing_region_cols].mean(axis=1, skipna=True)
        else:
            print("Warning: No regional pLDDT columns found; 'mean_plddt_regions' will not be added.")
    
    # Re-order columns
    cols = list(df.columns)
    plddt_cols = sorted([c for c in cols if 'plddt' in c])
    other_cols = [c for c in cols if 'plddt' not in c and c not in ['folder', 'model_index']]
    ordered_cols = ['folder', 'model_index'] + other_cols + plddt_cols
    final_cols = [c for c in ordered_cols if c in df.columns]
    df = df[final_cols]
    
    # save to csv file
    df.to_csv(output_path, index=False)
    
    print(f"Extracted metrics for {len(df)} models. Saved to {output_path}")

if __name__ == "__main__":
    # Define your base directory here (directory containing AF3 model folders)
    BASE_DIR = "/home/gfreimann/Seafile/Freimann_PhD/1_Projects/SolubleAf1503/alphafold3models"
    
    # Define your regions of interest here
    REGIONS_OF_INTEREST = [
        {'name': 'TM1_A', 'chain': 'A', 'start': 7, 'end': 30},
        {'name': 'TM2_A', 'chain': 'A', 'start': 255, 'end': 277},
        {'name': 'TM1_B', 'chain': 'B', 'start': 7, 'end': 30},
        {'name': 'TM2_B', 'chain': 'B', 'start': 255, 'end': 277}
    ]

    # Define your output path here
    output_path = "af3_metrics_summary_plddt_per_region.csv"

    main(BASE_DIR, REGIONS_OF_INTEREST, output_path)