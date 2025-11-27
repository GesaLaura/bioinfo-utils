# Automated Sequence Alignment via EMBOSS Needle

This Python script automates the pairwise sequence alignment of a specific target protein against a list of UniProt entries. It dynamically fetches protein sequences from the **UniProt REST API**, executes the local **EMBOSS Needle** tool (Needleman-Wunsch algorithm), and parses the output for **Identity**, **Similarity**, and **Score** metrics.

The script creates a detailed CSV file for statistical analysis and saves the full, human-readable alignment text files for every processed ID.

## Requirements

  - Python 3
  - [Biopython](https://biopython.org/)
  - [Requests](https://pypi.org/project/requests/)
  - **EMBOSS Suite** (The `needle` command must be installed and accessible in your system PATH). It can be installed via: `conda install -c bioconda emboss`

## Usage

The script is configured by modifying the variables at the top of the file.

1.  **Configure script variables:**
    Open the script and edit the configuration section:

      * `TARGET_SEQUENCE`: Paste your reference protein sequence string here.
      * `uniprot_ids`: Define the list of UniProt Accession IDs you wish to align against the target.
      * `csv_filename`: Specify the name for the summary statistics CSV.
      * `txt_output_dir`: Specify the directory name where individual alignment text files will be saved.

2.  **Run the script:**
    Once configured, run the script from your terminal:

    ```bash
    python3 run_needle_alignment.py
    ```



