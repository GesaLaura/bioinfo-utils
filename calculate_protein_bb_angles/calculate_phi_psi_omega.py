from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.vectors import calc_dihedral
import math
import csv

def rad_to_deg(angle):
    return math.degrees(angle) if angle is not None else None

def calculate_omega(residue1, residue2):
    """
    Calculate omega torsion angle between two adjacent residues.
    Omega is the dihedral angle between:
    CA (residue1) - C (residue1) - N (residue2) - CA (residue2)
    """
    try:
        CA_prev = residue1['CA'].get_vector()
        C_prev = residue1['C'].get_vector()
        N_next = residue2['N'].get_vector()
        CA_next = residue2['CA'].get_vector()
        omega = calc_dihedral(CA_prev, C_prev, N_next, CA_next)
        return omega
    except KeyError:
        return None

def main(pdb_file, output_csv):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    ppb = PPBuilder()

    with open(output_csv, mode='w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['ResNum', 'ResName', 'Phi (°)', 'Psi (°)', 'Omega (°)'])
        
        for pp in ppb.build_peptides(structure):
            phi_psi = pp.get_phi_psi_list()
            residues = pp
            for i, (residue, (phi, psi)) in enumerate(zip(residues, phi_psi)):
                res_id = residue.get_id()[1]
                res_name = residue.get_resname()
                phi_deg = rad_to_deg(phi)
                psi_deg = rad_to_deg(psi)

                # Omega needs two residues, so check if next residue exists
                if i < len(residues) - 1:
                    omega = calculate_omega(residues[i], residues[i+1])
                    omega_deg = rad_to_deg(omega)
                else:
                    omega_deg = None

                writer.writerow([res_id, res_name, phi_deg, psi_deg, omega_deg])

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python3 calculate_phi_psi_omega.py your_structure.pdb output.csv")
        sys.exit(1)
    pdb_file = sys.argv[1]
    output_csv = sys.argv[2]
    main(pdb_file, output_csv)
