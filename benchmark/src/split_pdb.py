"""
Script to organize PDB files into subdirectories based on protein categories.
Reads protein_categories.csv and organizes PDB files from pdb_files/ into
category-based subdirectories in a new directory called pdb_by_category/.

Author: Automated Script
Date: 2025-11-19
"""

import os
import shutil
import csv
from pathlib import Path
from collections import defaultdict

def organize_pdb_files():
    # Define base paths
    base_dir = Path(__file__).parent
    csv_file = base_dir / "protein_categories.csv"
    source_dir = base_dir / "pdb_files"
    target_dir = base_dir / "pdb_by_category"
    
    # Check if source files exist
    if not csv_file.exists():
        print(f"Error: CSV file not found at {csv_file}")
        return
    
    if not source_dir.exists():
        print(f"Error: Source directory not found at {source_dir}")
        return
    
    # Create target directory if it doesn't exist
    target_dir.mkdir(exist_ok=True)
    print(f"Created/verified target directory: {target_dir}")
    
    # Read CSV and build mapping of PDB -> Category
    pdb_to_category = {}
    categories = set()
    
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            pdb_id = row['PDB']
            category = row['CATEGORY']
            pdb_to_category[pdb_id] = category
            categories.add(category)
    
    print(f"\nFound {len(pdb_to_category)} protein entries across {len(categories)} categories:")
    print(f"Categories: {sorted(categories)}")
    
    # Create subdirectories for each category
    category_dirs = {}
    for category in categories:
        cat_dir = target_dir / category.replace(" ", "_")
        cat_dir.mkdir(exist_ok=True)
        category_dirs[category] = cat_dir
        print(f"  Created directory: {cat_dir}")
    
    # Copy/move PDB files to appropriate category directories
    stats = defaultdict(int)
    not_found = []
    
    print("\nOrganizing PDB files...")
    for pdb_id, category in pdb_to_category.items():
        source_file = source_dir / f"{pdb_id}.pdb"
        
        if source_file.exists():
            target_file = category_dirs[category] / f"{pdb_id}.pdb"
            shutil.copy2(source_file, target_file)
            stats[category] += 1
            if stats[category] % 50 == 0:
                print(f"  {category}: {stats[category]} files copied...")
        else:
            not_found.append(pdb_id)
    
    # Print summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    total_copied = 0
    for category in sorted(categories):
        count = stats[category]
        total_copied += count
        print(f"{category:15s}: {count:4d} files")
    
    print(f"\nTotal files copied: {total_copied}")
    
    if not_found:
        print(f"\nWarning: {len(not_found)} PDB files not found in source directory:")
        for pdb_id in sorted(not_found)[:10]:  # Show first 10
            print(f"  - {pdb_id}.pdb")
        if len(not_found) > 10:
            print(f"  ... and {len(not_found) - 10} more")
    
    print(f"\nAll PDB files organized in: {target_dir}")
    print("="*60)

if __name__ == "__main__":
    organize_pdb_files()