import requests
import subprocess
import re
import os
import csv
from Bio import SeqIO
from io import StringIO

# 1. Your Target Sequence
TARGET_SEQUENCE = "MLNFYSLLPFSALLTNLVFGIFILYIDPKSRLNRLYSLFTLSFAFWALGDFMVFMAYDQGSSLFWAKAASVGSTLSAAFAINFFVLLTKNRLANSKLMLLFYVPAAAFSLLSLGTDLITKGTKPVPWGLLQLRGPLYIPMTLFIVACMVASIILCYMYYRRTEKQDERKQMRIMIIGTSIPLFGGIVTQIVPIIMGFEMIPLSSTLSIIIVVAAALGVMRYRLMTPVSFSIRKKIIVSIMVASIIPLLMLSHVSI" 

# 2. List of UniProt IDs to fetch
uniprot_ids = [
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

# 3. Output CSV for the statistics
csv_filename = "summary_statistics_alignments.csv"

# 4. Output Directory for the text alignment files
txt_output_dir = "alignment_files"

# Create the directory if it doesn't exist
os.makedirs(txt_output_dir, exist_ok=True)

def get_sequence_from_uniprot(uniprot_id):
    """Fetches sequence from UniProt."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    try:
        response = requests.get(url)
        response.raise_for_status()
        fasta_io = StringIO(response.text)
        record = next(SeqIO.parse(fasta_io, "fasta"))
        return str(record.seq)
    except Exception as e:
        print(f"Error fetching {uniprot_id}: {e}")
        return None

def run_needle_cli(seq_a, seq_b):
    """
    Runs needle, returns the stats AND the full text output.
    """
    file_a = "temp_target.fasta"
    file_b = "temp_query.fasta"
    
    with open(file_a, "w") as f:
        f.write(f">Target\n{seq_a}")
    with open(file_b, "w") as f:
        f.write(f">Query\n{seq_b}")
        
    cmd = [
        "needle",
        "-asequence", file_a,
        "-bsequence", file_b,
        "-gapopen", "10.0",
        "-gapextend", "0.5",
        "-endopen", "10.0",   
        "-endextend", "0.5",  
        "-outfile", "stdout", 
        "-auto"               
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        output_text = result.stdout # This is the full alignment text
        
        # Parse the statistics
        ident_match = re.search(r"Identity:\s+\d+/\d+\s+\(([\d.]+)%\)", output_text)
        sim_match = re.search(r"Similarity:\s+\d+/\d+\s+\(([\d.]+)%\)", output_text)
        score_match = re.search(r"Score:\s+([\d.]+)", output_text)
        
        identity = ident_match.group(1) if ident_match else "0"
        similarity = sim_match.group(1) if sim_match else "0"
        score = score_match.group(1) if score_match else "0"
        
        return identity, similarity, score, output_text
        
    except subprocess.CalledProcessError as e:
        print(f"Error running Needle: {e}")
        return None, None, None, None
    finally:
        if os.path.exists(file_a): os.remove(file_a)
        if os.path.exists(file_b): os.remove(file_b)

results = []
print(f"Processing {len(uniprot_ids)} sequences...")

for uid in uniprot_ids:
    print(f"Processing {uid}...", end=" ")
    
    seq_query = get_sequence_from_uniprot(uid)
    
    if seq_query:
        ident, sim, score, full_text = run_needle_cli(TARGET_SEQUENCE, seq_query)
        
        if ident:
            print(f"Ident: {ident}% | Sim: {sim}%")
            
            # 1. Add to List for CSV
            results.append({
                "UniProt_ID": uid,
                "Identity": ident,
                "Similarity": sim,
                "Score": score
            })
            
            # 2. SAVE THE FULL ALIGNMENT TEXT
            txt_filename = os.path.join(txt_output_dir, f"{uid}_alignment.txt")
            with open(txt_filename, "w") as f:
                f.write(full_text)
                
        else:
            print("Failed to align.")
    else:
        print("Failed to download.")

# Export CSV
try:
    with open(csv_filename, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["UniProt_ID", "Identity", "Similarity", "Score"])
        writer.writeheader()
        writer.writerows(results)
    print(f"\nDone! \n1. Statistics saved to '{csv_filename}'\n2. Alignments saved in '{txt_output_dir}/'")
except Exception as e:
    print(f"Could not save CSV: {e}")