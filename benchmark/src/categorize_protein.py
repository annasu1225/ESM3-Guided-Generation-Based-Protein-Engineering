import pandas as pd
from Bio.PDB import PDBList, PDBParser
from Bio.PDB.Polypeptide import is_aa
import os
import warnings
import shutil

# Suppress PDB construction warnings for cleaner output
warnings.filterwarnings("ignore")

def get_protein_length(pdb_id, chain_id, pdb_dir="pdb_files"):
    """
    Downloads the PDB file, ensures it has a [pdb_id].pdb filename format, 
    and calculates the length of the specified chain.
    """
    pdbl = PDBList()
    parser = PDBParser(QUIET=True)
    
    # Ensure directory exists
    os.makedirs(pdb_dir, exist_ok=True)
    
    try:
        # Download the PDB file. 
        # By default, Biopython saves it as pdbXXXX.ent (or similar depending on OS)
        original_filepath = pdbl.retrieve_pdb_file(pdb_id, pdir=pdb_dir, file_format="pdb")
        
        # --- NEW Renaming Logic: Force filename to be simply {pdb_id}.pdb ---
        # os.path.split gets the filename from the path
        # We ignore whatever Biopython named it and enforce our own standard name
        final_filename = f"{pdb_id}.pdb"
        final_filepath = os.path.join(pdb_dir, final_filename)
        
        # Rename the file if it doesn't match our desired name
        # (Biopython might return an absolute path, so we compare carefully)
        if os.path.abspath(original_filepath) != os.path.abspath(final_filepath):
            shutil.move(original_filepath, final_filepath)
            
        # Parse the structure using the final .pdb path
        structure = parser.get_structure(pdb_id, final_filepath)
        
        # Get the specific model (usually model 0) and chain
        model = structure[0]
        
        if chain_id not in model:
            # Clean up the file if the chain is invalid for a clean rerun
            # os.remove(final_filepath) # Optional: keep file even if chain is missing
            print(f"Warning: Chain {chain_id} not found in {pdb_id}. Skipping.")
            return None
            
        chain = model[chain_id]
        
        # Count only standard amino acid residues
        residue_count = 0
        for residue in chain:
            if is_aa(residue, standard=True):
                residue_count += 1
                
        return residue_count

    except Exception as e:
        print(f"Error processing {pdb_id} chain {chain_id}: {e}")
        return None

def categorize_length(length):
    """
    Categorizes proteins based on sequence length criteria:
    < 50: Very Small, 50-100: Small, 100-200: Medium, > 200: Large.
    """
    if length < 50:
        return "Very Small"
    elif 50 <= length <= 100:
        return "Small"
    elif 100 < length <= 200:
        return "Medium"
    else: # length > 200
        return "Large"

def main():
    # Note: Using the file name provided in the original request's context.
    csv_file = "consolidated_protddg_data_with_source.csv"
    
    print(f"Reading {csv_file}...")
    try:
        df = pd.read_csv(csv_file)
    except FileNotFoundError:
        print(f"Error: Could not find {csv_file}")
        return

    # 1. Extract unique PDB and Chain combinations
    # The snippet shows PDB IDs combined with pipes (e.g., 1PGA|1EM7|2GB1).
    # We take the first PDB ID and the corresponding CHAIN.
    print("Extracting unique protein chains...")
    
    # Clean up entries where PDB ID might be a list
    df['PDB_CLEAN'] = df['PDB'].astype(str).apply(lambda x: x.split('|')[0].strip())
    
    unique_proteins = df[['PDB_CLEAN', 'CHAIN']].drop_duplicates().reset_index(drop=True)
    
    # Rename columns for clarity in the loop
    unique_proteins.columns = ['PDB', 'CHAIN']
    
    print(f"Found {len(unique_proteins)} unique protein chains (using first PDB ID if multiple provided).")
    
    results = {
        "Very Small": [],
        "Small": [],
        "Medium": [],
        "Large": []
    }

    # 2. Process each unique protein
    print("\nFetching PDBs and calculating lengths...")
    print("-" * 50)
    
    for index, row in unique_proteins.iterrows():
        pdb_id = row['PDB'].strip()
        chain_id = row['CHAIN'].strip()
        
        # Skip invalid entries like 'nan'
        if pdb_id == 'nan' or pdb_id == 'None' or pd.isna(chain_id):
            continue
            
        length = get_protein_length(pdb_id, chain_id)
        
        if length is not None:
            category = categorize_length(length)
            results[category].append((pdb_id, chain_id, length))
            print(f"[{index+1}/{len(unique_proteins)}] {pdb_id} (Chain {chain_id}): {length} residues -> {category}")
        else:
            print(f"[{index+1}/{len(unique_proteins)}] {pdb_id} (Chain {chain_id}): Failed to calculate length/Chain not found")

    # 3. Output the summary and save results
    print("\n" + "="*30)
    print("FINAL SUMMARY")
    print("="*30)
    
    category_order = ["Very Small", "Small", "Medium", "Large"]
    
    for category in category_order:
        proteins = results[category]
        print(f"\n{category.upper()} PROTEINS ({len(proteins)}):")
        print("-" * 20)
        if not proteins:
            print("None found.")
        for p in proteins:
            print(f"{p[0]} (Chain {p[1]}): {p[2]} residues")

    # 4. Save to CSV
    output_data = []
    for cat, items in results.items():
        for item in items:
            output_data.append({'PDB': item[0], 'CHAIN': item[1], 'LENGTH': item[2], 'CATEGORY': cat})
    
    pd.DataFrame(output_data).to_csv("protein_categories.csv", index=False)
    print("\nDetailed results saved to 'protein_categories.csv'")

if __name__ == "__main__":
    main()