# Ramachandran Plot Generator

**Author:** Gesa Freimann

Tool to automatically fetch protein structures from the Protein Data Bank (PDB), calculate backbone $\phi$ (phi) and $\psi$ (psi) dihedral angles, and visualize them on a multi-protein Ramachandran plot.

## Features
- **Automated Fetching:** Automatically downloads and extracts structural data (`.ent` files) directly from the RCSB database.
- **Multi-Protein Support:** Processes a list of PDB IDs simultaneously, overlaying them on a single plot with distinct colors.
- **Parsing:** Accepts CSV or TSV inputs, processes the first column, and automatically filters out header rows or invalid IDs.

## Installation

Ensure you have Python 3 installed along with the required structural biology and visualization libraries:

```bash
pip install biopython matplotlib
```

## Usage
Run the script from your terminal using the -i (or --input) flag to provide your file containing PDB IDs.
```bash
python plot_ramachandran.py -i input_pdb_ids.csv
```

## Command-Line Arguments
* -i, --input (Required): Path to the input CSV or TSV file.
* -o, --output (Optional): Path/filename for the generated output image (Defaults to ramachandran_plot.png).

## Example: saving as pdf

```bash
python plot_ramachandran.py -i input_pdb_ids.csv -o high_quality_plot.pdf
```

## Input file format
The script looks at the first column of your file and extracts valid 4-character alphanumeric PDB IDs. Headers and other metadata columns are automatically ignored during parsing.

Example input_pdb_ids.csv

```bash
PDB_ID,Protein_Name,Notes
1crn,Crambin,Plant protein
1mbn,Myoglobin,Oxygen transport
1tup,p53,Tumor suppressor core domain
```

