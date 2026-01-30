"""
Boltz Scoring Utilities for ESM3 Guided Generation

This module provides scoring functions for protein-ligand binding affinity
prediction using Boltz-1 model during ESM3 guided generation.

Author: Adapted from Anna Su's original implementation
Date: 2025-12-06
"""

import os
import json
import subprocess
import shutil
import hashlib
import uuid
from datetime import datetime
from pathlib import Path
from typing import Optional, Dict
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

from esm.sdk.api import ESMProtein
from guided_generation import GuidedDecodingScoringFunction


def get_cache_path(seq: str, smiles: str, cache_dir: str) -> str:
    """Generate cache path based on sequence and SMILES hash."""
    combined = f"{seq}_{smiles}"
    h = hashlib.sha256(combined.encode()).hexdigest()
    return os.path.join(cache_dir, f"boltz_affinity_{h[:16]}.json")


def compute_boltz_affinity(
    seq: str,
    smiles: str,
    cache_dir: str,
    timeout_sec: int = 3600,
    cleanup_tmp: bool = True,
    verbose_boltz: bool = False,
    step: int = 1,
    log_file_path: Optional[str] = None,
    score_type: str = "affinity"  # "affinity" or "binary"
) -> Optional[float]:
    """
    Compute protein-ligand binding affinity using Boltz-1.
    
    Args:
        seq: Protein sequence
        smiles: Ligand SMILES string
        cache_dir: Directory for caching results
        timeout_sec: Timeout for Boltz execution
        cleanup_tmp: Whether to clean up temporary files
        verbose_boltz: Whether to print verbose output
        step: Current generation step (for logging)
        log_file_path: Path to log file
    
    Returns:
        Predicted binding affinity score, or None if prediction fails
    """
    import joblib
    
    # Check cache first
    cache_path = get_cache_path(seq, smiles, cache_dir)
    if os.path.exists(cache_path):
        try:
            cached_data = joblib.load(cache_path)
            return cached_data
        except:
            pass
    
    # Create temporary directory for Boltz run
    tmp_dir = os.path.join(cache_dir, f"boltz_run_{uuid.uuid4().hex[:8]}")
    os.makedirs(tmp_dir, exist_ok=True)
    
    affinity = None
    
    try:
        # Create YAML file for Boltz input
        yaml_path = os.path.join(tmp_dir, "protein_ligand.yaml")
        yaml_content = (
            "version: 1\n"
            "sequences:\n"
            "  - protein:\n"
            "      id: A\n"
            f"      sequence: '{seq}'\n"
            "  - ligand:\n"
            "      id: B\n"
            f"      smiles: '{smiles}'\n"
            "properties:\n"
            "  - affinity:\n"
            "      binder: B\n"
        )
        
        with open(yaml_path, "w") as f:
            f.write(yaml_content)
        
        # Run Boltz prediction
        boltz_exe = shutil.which("boltz")
        if not boltz_exe:
            raise RuntimeError("Boltz executable not found in PATH")
        
        cmd = [
            boltz_exe, "predict", yaml_path,
            "--out_dir", tmp_dir,
            "--use_msa_server"
        ]
        
        proc = subprocess.run(
            cmd,
            cwd=tmp_dir,
            capture_output=True,
            text=True,
            timeout=timeout_sec
        )
        
        # Log Boltz output for first step if requested
        if step == 1 and log_file_path and verbose_boltz:
            with open(log_file_path, 'a') as f:
                f.write(f"\n--- BOLTZ Log for sequence {seq[:15]}... ---\n")
                f.write(f"STDOUT:\n{proc.stdout}\n")
                f.write(f"STDERR:\n{proc.stderr}\n")
                f.write("---------------------------------------------------\n")
        
        if proc.returncode == 0:
            # Parse affinity from JSON output
            json_files = list(Path(tmp_dir).rglob("affinity_*.json"))
            if json_files:
                json_path = max(json_files, key=lambda p: p.stat().st_mtime)
                with open(json_path, "r") as jf:
                    data = json.load(jf)
                
                # Choose which score to return based on score_type
                if score_type == "binary":
                    # affinity_probability_binary: higher = more likely binder (0 to 1)
                    affinity = float(data.get("affinity_probability_binary", None))
                else:
                    # affinity_pred_value: log10(IC50) where LOWER = tighter binding
                    # We NEGATE so that higher = better (consistent with guided generation)
                    raw_affinity = float(data.get("affinity_pred_value", None))
                    affinity = -raw_affinity if raw_affinity is not None else None
                
                if verbose_boltz and affinity is not None:
                    print(f"[BOLTZ] Predicted affinity score: {affinity:.4f} for seq: {seq[:20]}...")
        else:
            if log_file_path:
                with open(log_file_path, 'a') as f:
                    f.write(f"[ERROR] Boltz failed for seq {seq[:15]}...\n")
                    f.write(f"Return code: {proc.returncode}\n")
                    f.write(f"STDERR: {proc.stderr}\n")
    
    except subprocess.TimeoutExpired:
        if log_file_path:
            with open(log_file_path, 'a') as f:
                f.write(f"[ERROR] Boltz timeout for seq {seq[:15]}...\n")
    
    except Exception as e:
        if log_file_path:
            with open(log_file_path, 'a') as f:
                f.write(f"[ERROR] Exception during Boltz run for seq {seq[:15]}: {e}\n")
        if verbose_boltz:
            print(f"[ERROR] Boltz prediction failed: {e}")
    
    finally:
        # Clean up temporary directory
        if cleanup_tmp:
            shutil.rmtree(tmp_dir, ignore_errors=True)
    
    # Cache the result if successful
    if affinity is not None:
        joblib.dump(affinity, cache_path)
    
    return affinity


class BoltzScorer(GuidedDecodingScoringFunction):
    """
    Scoring function for ESM3 guided generation using Boltz-1 binding affinity.
    
    Higher affinity values indicate better binding.
    """
    
    def __init__(
        self,
        smiles: str,
        cache_dir: str,
        timeout_sec: int = 3600,
        cleanup_tmp: bool = True,
        verbose_boltz: bool = False,
        score_type: str = "affinity"  # "affinity" for optimization, "binary" for hit discovery
    ):
        """
        Initialize Boltz scorer.
        
        Args:
            smiles: SMILES string of the target ligand
            cache_dir: Directory for caching Boltz results
            timeout_sec: Timeout for Boltz execution
            cleanup_tmp: Whether to clean up temporary files
            verbose_boltz: Whether to print verbose output
        """
        self.smiles = smiles
        self.cache_dir = cache_dir
        self.timeout_sec = timeout_sec
        self.cleanup_tmp = cleanup_tmp
        self.verbose_boltz = verbose_boltz
        self.score_type = score_type
        
        # Define invalid amino acids (same as FoldX implementation)
        self.invalid_amino_acids = {'B', 'J', 'O', 'U', 'X', 'Z'}
    
    def __call__(self, protein: ESMProtein, step: int, log_file_path: str) -> float:
        """
        Score a protein based on predicted binding affinity to the target ligand.
        
        Args:
            protein: ESMProtein object to score
            step: Current generation step
            log_file_path: Path to log file
        
        Returns:
            Score (higher is better). Returns -inf for invalid sequences.
        """
        sequence = protein.sequence
        
        # Check for invalid amino acids
        if any(char in self.invalid_amino_acids for char in sequence):
            print(f"[SCORE] Invalid sequence with non-standard amino acid found. Score: -inf")
            return float("-inf")
        
        # Compute binding affinity using Boltz
        affinity = compute_boltz_affinity(
            seq=sequence,
            smiles=self.smiles,
            cache_dir=self.cache_dir,
            timeout_sec=self.timeout_sec,
            cleanup_tmp=self.cleanup_tmp,
            verbose_boltz=self.verbose_boltz,
            step=step,
            log_file_path=log_file_path,
            score_type=self.score_type
        )
        
        if affinity is None:
            return float("-inf")
        
        # Return affinity as score
        # For score_type="affinity": already negated in compute_boltz_affinity, so higher = better
        # For score_type="binary": higher probability = better
        return float(affinity)


def plot_affinity_history(all_scores_history: Dict, save_path: str):
    """
    Creates and saves a box plot of binding affinity scores over generation steps.
    
    Args:
        all_scores_history: Dictionary mapping step numbers to lists of scores
        save_path: Path to save the plot
    """
    plot_data = []
    for step, scores in all_scores_history.items():
        for score in scores:
            if score is not None and score != float('-inf'):
                plot_data.append({'step': step, 'affinity': score})
    
    if not plot_data:
        print("[WARN] No valid scores to plot.")
        return

    df = pd.DataFrame(plot_data)
    
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Create box plot
    sns.boxplot(x='step', y='affinity', data=df, ax=ax, palette="viridis")
    
    # Overlay individual points
    sns.stripplot(
        x='step', 
        y='affinity', 
        data=df, 
        ax=ax, 
        color='black', 
        alpha=0.6, 
        jitter=True
    )
    
    ax.set_title('Distribution of Predicted Binding Affinity Scores per Generation Step', fontsize=16)
    ax.set_xlabel('Generation Step', fontsize=12)
    ax.set_ylabel('Predicted Binding Affinity (pK)', fontsize=12)
    ax.grid(True, alpha=0.3)
    
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"\n[INFO] Plot saved to: {save_path}")
    plt.close()