#!/usr/bin/env python3
"""
Ramachandran Plot Generator
Author: Gesa Freimann

Description:
This script reads a CSV or TSV file containing PDB IDs, fetches the corresponding 
PDB files from the RCSB database, calculates the backbone phi and psi dihedral 
angles, and generates a combined Ramachandran plot.
"""

import math
import os
import csv
import argparse
import matplotlib.pyplot as plt
from Bio.PDB import PDBList, PDBParser, PPBuilder

def read_pdb_ids(file_path):
    """
    Reads a CSV or TSV file and extracts PDB IDs from the first column.
    Skips headers automatically by checking if the string is exactly 4 alphanumeric characters.
    """
    pdb_ids = []
    
    if not os.path.isfile(file_path):
        print(f"Error: The file '{file_path}' does not exist.")
        return pdb_ids

    # Determine the delimiter based on the file extension
    delimiter = '\t' if file_path.lower().endswith('.tsv') else ','
    
    with open(file_path, mode='r', encoding='utf-8') as file:
        reader = csv.reader(file, delimiter=delimiter)
        for row in reader:
            if row: # Skip empty rows
                # Extract the first column and clean up whitespace
                potential_id = row[0].strip()
                
                # Simple validation: PDB IDs are exactly 4 alphanumeric characters
                if len(potential_id) == 4 and potential_id.isalnum():
                    pdb_ids.append(potential_id)
                    
    return pdb_ids

def fetch_and_parse_pdb(pdb_id, download_dir="pdb_files"):
    """
    Fetches the PDB file from the RCSB database and parses its structure.
    """
    pdbl = PDBList()
    parser = PDBParser(QUIET=True) # Suppress standard structural warnings
    
    os.makedirs(download_dir, exist_ok=True)
    
    print(f"Fetching PDB: {pdb_id.upper()}...")
    file_path = pdbl.retrieve_pdb_file(pdb_id, pdir=download_dir, file_format='pdb')
    
    if not os.path.exists(file_path):
        print(f"Error: Could not retrieve {pdb_id}")
        return None
        
    structure = parser.get_structure(pdb_id, file_path)
    return structure

def calculate_phi_psi(structure):
    """
    Extracts the phi and psi angles for all residues in a given structure.
    Returns two lists: one for phi angles, one for psi angles (in degrees).
    """
    phi_angles = []
    psi_angles = []
    
    ppb = PPBuilder()
    
    for model in structure:
        for chain in model:
            polypeptides = ppb.build_peptides(chain)
            for poly in polypeptides:
                phi_psi_list = poly.get_phi_psi_list()
                
                for phi, psi in phi_psi_list:
                    if phi is not None and psi is not None:
                        phi_angles.append(math.degrees(phi))
                        psi_angles.append(math.degrees(psi))
                        
    return phi_angles, psi_angles

def plot_ramachandran(pdb_ids, output_filename):
    """
    Main function: downloading, calculating, plotting, and saving.
    """
    plt.rcParams.update({'font.size': 14})
    
    plt.figure(figsize=(10, 10))
    
    plt.title('Ramachandran Plot', fontsize=22, fontweight='bold', pad=15)
    plt.xlabel(r'$\phi$ (degrees)', fontsize=18, labelpad=10)
    plt.ylabel(r'$\psi$ (degrees)', fontsize=18, labelpad=10)
    
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    
    plt.axhline(0, color='black', linewidth=1.2)
    plt.axvline(0, color='black', linewidth=1.2)
    plt.grid(color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    
    for pdb_id in pdb_ids:
        structure = fetch_and_parse_pdb(pdb_id)
        
        if structure:
            phi, psi = calculate_phi_psi(structure)
            if phi and psi:
                plt.scatter(phi, psi, s=35, alpha=0.7, edgecolors='none', label=f"{pdb_id.upper()} (n={len(phi)})")
            else:
                print(f"Warning: No valid backbone angles found for {pdb_id.upper()}")
    
    plt.legend(loc='upper right', framealpha=1.0, fontsize=14)
    plt.tight_layout()
    
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"\nSuccess! Highly visible plot saved to: {output_filename}")
    
    plt.close()

if __name__ == "__main__":
    print("Starting Ramachandran Plot generation by Gesa Freimann...\n")
    
    parser = argparse.ArgumentParser(description="Plot and save a large-text Ramachandran plot from a CSV/TSV list.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input CSV or TSV file.")
    parser.add_argument("-o", "--output", default="ramachandran_plot.png", help="Path/filename for the output image.")
    
    args = parser.parse_args()
    
    my_pdb_list = read_pdb_ids(args.input)
    
    if not my_pdb_list:
        print("Error: No valid PDB IDs found in the provided file. Exiting.")
    else:
        print(f"Found {len(my_pdb_list)} valid PDB IDs to process: {', '.join(my_pdb_list)}\n")
        plot_ramachandran(my_pdb_list, args.output)