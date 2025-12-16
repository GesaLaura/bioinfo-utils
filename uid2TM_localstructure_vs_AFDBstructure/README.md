#Automated Structural Alignment via TM-align

This Python script automates the one-to-many structural comparison of a specific local reference PDB file against a list of UniProt entries. It dynamically fetches predicted protein structures from the **AlphaFold Protein Structure Database**, calculates the **Radius of Gyration (Rg)** for both the reference and query structures, executes the **TM-align** algorithm, and parses the output for **RMSD**, **TM-Score**, and **Aligned Length**.

The script generates a consolidated CSV file for statistical analysis and saves the raw PDB files and TM-align output logs in a dedicated working directory.

##Requirements* Python 3
* [Biopython](https://biopython.org/) (for PDB parsing and Rg calculation)
* [Pandas](https://pandas.pydata.org/)
* [Requests](https://pypi.org/project/requests/)
* [NumPy](https://numpy.org/)
* **TM-align** (The `TMalign` executable must be installed and accessible).
* *Note: If `TMalign` is not in your system PATH, you must specify the full path in the script configuration.*
* Source: [Zhang Lab TM-align](https://zhanggroup.org/TM-align/)

#### For uid2TM_conda.py script
Install TM-align via Conda following these instructions :
> conda install -c schrodinger pymol
> conda install -c schrodinger pymol-psico
> conda install -c speleo3 tmalign

After installation, add the path to TM-align to your PATH. 
Replace /path/to/tmalign with the actual directory containing the TM-align executable:
> export PATH="$PATH:/path/to/tmalign"

Save and exit the editor, then reload the shell configuration:
> source ~/.bashrc

Confirm the installation:
> TMalign -h

##UsageThe script is configured by modifying the variables at the top of the file to match your specific experiment.

1. **Configure script variables:**
Open the script and edit the configuration section:
* `REFERENCE_PDB_PATH`: The absolute path to your local reference PDB file (the model you are comparing against).
* `UNIPROT_IDS_TO_COMPARE`: Define the list of UniProt Accession IDs you wish to fetch and align against the reference.
* `WORKING_DIR`: Specify the directory name where downloaded PDBs and alignment logs will be stored.
* `OUTPUT_CSV_NAME`: Specify the filename for the final summary statistics CSV.
* `TMALIGN_CMD`: The command to call TM-align (default is `"TMalign"`, but can be changed to a full path if necessary).


2. **Run the script:**
Once configured, run the script from your terminal:
```bash
python3 run_structure_comparison.py

```



##OutputThe script will generate a CSV file (defined in `OUTPUT_CSV_NAME`) containing the following columns for every processed ID:

* **Query_ID**: The UniProt ID of the fetched structure.
* **Rg_Query**: Radius of Gyration of the fetched structure.
* **Rg_Reference**: Radius of Gyration of your local reference.
* **RMSD**: Root-mean-square deviation of the alignment.
* **TM-Score**: The structural similarity score (0.0 - 1.0).
* **Aligned_Length**: Number of aligned residues.



