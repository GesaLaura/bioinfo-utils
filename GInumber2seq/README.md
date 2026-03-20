# NCBI Sequence Mapper

A Python utility to map protein identifiers (GI or Accession numbers) to their full FASTA sequences and updated Accession IDs.

## Features
- **Direct Lookup**: Uses NCBI Entrez `efetch` for fast retrieval via ID.
- **BLAST Backup**: If a fragment sequence is provided, the script can perform a BLASTp search to find the parent protein if the ID is deprecated or missing.
- **Rate-Limiting**: Built-in delays to comply with NCBI server policies.

## Installation of Requirements
```bash
pip install pandas biopython
```

## Usage
Run the script from the terminal:
```bash
python map_bioseq.py -i input.csv -o results.csv -e your.email@example.com --gi_col gi_id --seq_col fragment_sequence
```

### Arguments:
| Argument | Description | Required |
| :--- | :--- | :--- |
| `-i` / `--input` | Path to the input CSV. | Yes |
| `-o` / `--output` | Path for the output CSV. | Yes |
| `-e` / `--email` | Email address for NCBI identification. | Yes |
| `--gi_col` | The column name for GIs/Accessions. | No (Default: 'gi') |
| `--seq_col` | The column name for fragments (enables BLAST backup). | No |


### 3. Test Data (`test_data.csv`)
Save this as `test_data.csv` in your repository. It includes one valid GI, one old GI that might require BLAST (if you provide a sequence), and one placeholder.

```csv
gi_id,sequence_fragment
160879611,QKIERPLDILSNGVMQISNGNLEHRIDYEYQDEFAPVCADFNEMAARLKASVE
150016677,KKLIRPLELLSYGAEQIKNGNLDFEMNYESDDEFGQVCGDFDEMRLRLKHSVD
163725388,KRILAPIRMLAAGTKDIAAQRFQVRLPISRSDELGRLASDFNTMAGRLEKSTQ
```

### Notes
The BLAST runs can take a while. You can speeden that up by using local BLAST runs or running the script at times of low usage of NCBI. 

### Author
Gesa Laura Freimann,
March 2026
