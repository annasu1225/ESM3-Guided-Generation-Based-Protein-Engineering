# -*- coding: utf-8 -*-
import pandas as pd
import os
from pathlib import Path

# Define paths
original_data_dir = Path("Original_Data")
output_file = "consolidated_protddg_data.csv"

# Files to process (excluding p53.csv)
csv_files = [
    "broom.csv",
    "korpm.csv", 
    "myoglobin.csv",
    "ptmul.csv",
    "s2648.csv",
    "ssym.csv",
    "vb1432.csv"
]

# Store all dataframes
all_dfs = []

print("="*80)
print("CONSOLIDATING PROTDDG-BENCH ORIGINAL DATA")
print("="*80)

for csv_file in csv_files:
    file_path = original_data_dir / csv_file
    
    if not file_path.exists():
        print(f"\nSkipping {csv_file} - file not found")
        continue
    
    print(f"\nReading: {csv_file}")
    
    # Read the CSV file
    # Handle the case where ptmul.csv has # in the header
    df = pd.read_csv(file_path, comment=None)
    
    # Clean column names (remove # if present)
    df.columns = df.columns.str.replace('#', '')
    
    print(f"   Original columns: {list(df.columns)}")
    print(f"   Original records: {len(df)}")
    
    # Standardize column names
    # Map CLUID to CLID
    if 'CLUID' in df.columns:
        df = df.rename(columns={'CLUID': 'CLID'})
    
    # Map MUTS to MUT
    if 'MUTS' in df.columns:
        df = df.rename(columns={'MUTS': 'MUT'})
    
    # Select only the required columns
    required_cols = ['CLID', 'PDB', 'CHAIN', 'MUT', 'DDG']
    
    # Check if all required columns exist
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"   WARNING: Missing columns: {missing_cols}")
        continue
    
    # Extract only the required columns
    df_selected = df[required_cols].copy()
    
    # Add source file information
    df_selected['SOURCE_FILE'] = csv_file.replace('.csv', '')
    
    print(f"   Extracted {len(df_selected)} records")
    
    all_dfs.append(df_selected)

# Combine all dataframes
if not all_dfs:
    print("\nNo data found!")
    exit(1)

print("\n" + "="*80)
print("COMBINING AND SORTING DATA")
print("="*80)

combined_df = pd.concat(all_dfs, ignore_index=True)
print(f"\nTotal records before grouping: {len(combined_df)}")

# Sort by PDB to group same PDB entries together
combined_df_sorted = combined_df.sort_values(by=['PDB', 'CHAIN', 'MUT']).reset_index(drop=True)

print(f"Total records after sorting: {len(combined_df_sorted)}")

# Display summary statistics
print("\n" + "="*80)
print("SUMMARY BY SOURCE FILE")
print("="*80)
source_summary = combined_df_sorted.groupby('SOURCE_FILE').size().sort_values(ascending=False)
for source, count in source_summary.items():
    print(f"{source:15s}: {count:5d} records")

print("\n" + "="*80)
print("SUMMARY BY PDB")
print("="*80)
pdb_summary = combined_df_sorted.groupby('PDB').size().sort_values(ascending=False).head(20)
print(f"\nTop 20 PDBs by number of mutations:")
for pdb, count in pdb_summary.items():
    print(f"  {pdb}: {count} mutations")

# Save to CSV
print("\n" + "="*80)
print("SAVING RESULTS")
print("="*80)

# Save with SOURCE_FILE column
output_with_source = "consolidated_protddg_data_with_source.csv"
combined_df_sorted.to_csv(output_with_source, index=False)
print(f"Saved to: {output_with_source}")
print(f"  Columns: {list(combined_df_sorted.columns)}")

# Save without SOURCE_FILE column (only requested columns)
combined_df_final = combined_df_sorted[['CLID', 'PDB', 'CHAIN', 'MUT', 'DDG']].copy()
combined_df_final.to_csv(output_file, index=False)
print(f"\nSaved to: {output_file}")
print(f"  Columns: {list(combined_df_final.columns)}")
print(f"  Total records: {len(combined_df_final)}")

print("\n" + "="*80)
print("CONSOLIDATION COMPLETE!")
print("="*80)

# Display first few rows
print("\nFirst 10 rows of consolidated data:")
print(combined_df_final.head(10).to_string(index=False))

print("\nLast 10 rows of consolidated data:")
print(combined_df_final.tail(10).to_string(index=False))
