# Extract AlphaFold 3 Confidence Metrics

This Python script traverses a directory of AlphaFold 3 server outputs to extract key confidence metrics. It parses `*.json` files for **ipTM**, **pTM**, **ranking scores**, and **PAE**, and parses the corresponding `*.cif` files to calculate **global pLDDT** and **region-specific pLDDT** scores.

All extracted metrics are compiled into a single, comprehensive CSV file for analysis.

## Requirements
-   Python 3
-   [Pandas](https://pandas.pydata.org/)
-   [NumPy](https://numpy.org/)
-   [Biopython](https://biopython.org/)

## Usage

The script is configured by modifying the variables inside the `if __name__ == "__main__":` block at the bottom of the script.

1.  **Configure script variables:**
    Open the script file and edit the variables at the very bottom:
    * `BASE_DIR`: Set this to the root directory containing your AlphaFold 3 model folders.
    * `REGIONS_OF_INTEREST`: Define a list of dictionaries for each protein region you want to analyze.
    * `output_path`: Specify the name for the output CSV file.

2.  **Run the script:**
    Once configured, run the script from your terminal:
    ```bash
    python3 analysis_af3_models_from_afserver_plddtregion.py
    ```



