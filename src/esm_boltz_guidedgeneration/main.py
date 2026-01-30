"""
ESM3-Boltz Guided Generation - Main Script

This script performs ESM3 guided generation to design protein sequences
with enhanced binding affinity to a specified ligand using Boltz-1 scoring.

Usage:
    python main.py --smiles "CCO" --seq_length 256 --masking_percentage 0.4 \\
                   --num_decoding_steps 32 --num_samples_per_step 20

Author: Adapted to match FoldX implementation structure
Date: 2025-12-06
"""

import os
import random
import time
import argparse
import warnings
from typing import Optional

import torch
from esm.models.esm3 import ESM3
from esm.sdk.api import ESMProtein
from esm.utils.structure.protein_chain import ProteinChain

# Import custom modules (must be in same directory or PYTHONPATH)
import sys
sys.path.insert(0, os.path.dirname(__file__))
from guided_generation import ESM3GuidedDecoding
from boltz_scoring_utils import BoltzScorer, plot_affinity_history


warnings.simplefilter(action='ignore', category=FutureWarning)


# --- Default Settings ---
TIMEOUT_SEC = 3600
CLEANUP_TMP = True
VERBOSE_BOLTZ = False

DEFAULT_LOG_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "logs", "boltz")
)

DEFAULT_RESULTS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "results", "boltz")
)


def get_masked_sequence(wildtype_pdb: str, wildtype_chain: str, masking_percentage: float) -> str:
    """
    Load a wildtype sequence from PDB and create a masked version.
    
    Args:
        wildtype_pdb: PDB ID (e.g., "1ABC")
        wildtype_chain: Chain ID (e.g., "A")
        masking_percentage: Fraction of residues to mask (0.0 to 1.0)
    
    Returns:
        Masked sequence string with '_' characters
    """
    print("[INFO] Loading wildtype from PDB: {}, Chain: {}".format(wildtype_pdb, wildtype_chain))
    
    wildtype_protein = ESMProtein.from_protein_chain(
        ProteinChain.from_rcsb(wildtype_pdb, chain_id=wildtype_chain)
    )
    
    wt_sequence = wildtype_protein.sequence
    maskable_indices = list(range(len(wt_sequence)))
    num_to_mask = int(len(maskable_indices) * masking_percentage)
    indices_to_mask = random.sample(maskable_indices, num_to_mask)
    
    refinement_template_list = list(wt_sequence)
    for i in indices_to_mask:
        refinement_template_list[i] = '_'
    
    masked_seq = "".join(refinement_template_list)
    
    print("[INFO] Wildtype sequence length: {}".format(len(wt_sequence)))
    print("[INFO] Masked {} positions ({:.1f}%)".format(num_to_mask, masking_percentage * 100))
    print("[INFO] Masked sequence: {}...".format(masked_seq[:80]))
    
    return masked_seq


def main():
    parser = argparse.ArgumentParser(
        description="Run ESM3 Guided Generation with Boltz-1 binding affinity scoring.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    parser.add_argument(
        "--smiles",
        type=str,
        required=True,
        help="SMILES string of the target ligand for binding affinity prediction."
    )
    
    # Sequence initialization options
    sequence_group = parser.add_mutually_exclusive_group(required=True)
    sequence_group.add_argument(
        "--seq_length",
        type=int,
        help="Length of de novo protein sequence to generate (fully masked start)."
    )
    sequence_group.add_argument(
        "--wildtype",
        nargs=2,
        metavar=("PDB_ID", "CHAIN_ID"),
        help="Start from wildtype PDB structure (e.g., --wildtype 1ABC A)."
    )
    
    # Generation parameters
    parser.add_argument(
        "--masking_percentage",
        type=float,
        default=0.40,
        help="Percentage of residues to mask (only used with --wildtype). Default: 0.40"
    )
    parser.add_argument(
        "--num_decoding_steps",
        type=int,
        default=32,
        help="Number of generation steps. Default: 32"
    )
    parser.add_argument(
        "--num_samples_per_step",
        type=int,
        default=20,
        help="Number of candidates to generate per step. Default: 20"
    )
    
    # Output options
    parser.add_argument(
        "--log_file",
        type=str,
        default=None,
        help="Path to generation log file. Defaults to timestamped file in logs/boltz/"
    )
    parser.add_argument(
        "--cache_dir",
        type=str,
        default=None,
        help="Directory for caching Boltz results. Default: ./boltz_cache"
    )
    parser.add_argument(
        "--score_type",
        type=str,
        choices=["affinity", "binary"],
        default="affinity",
        help="Score type: 'affinity' for lead optimization (log10 IC50), 'binary' for hit discovery (probability). Default: affinity"
    )
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.wildtype and not (0.0 < args.masking_percentage < 1.0):
        parser.error("--masking_percentage must be between 0.0 and 1.0")
    
    if args.seq_length and args.seq_length < 1:
        parser.error("--seq_length must be positive")
    
    # Setup
    os.environ["TOKENIZERS_PARALLELISM"] = "false"
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    
    # Setup cache directory
    if args.cache_dir:
        cache_dir = args.cache_dir
    else:
        cache_dir = os.path.join(os.getcwd(), "boltz_cache")
    os.makedirs(cache_dir, exist_ok=True)
    
    # --- 1. Create Starting Sequence ---
    if args.wildtype:
        wildtype_pdb, wildtype_chain = args.wildtype
        masked_seq = get_masked_sequence(wildtype_pdb, wildtype_chain, args.masking_percentage)
        seq_identifier = "{}_{}".format(wildtype_pdb, wildtype_chain)
    else:
        # De novo generation - fully masked sequence
        masked_seq = "_" * args.seq_length
        seq_identifier = "denovo_{}".format(args.seq_length)
        print("[INFO] Starting de novo generation with sequence length: {}".format(args.seq_length))
    
    starting_protein = ESMProtein(sequence=masked_seq)
    
    # --- 2. Setup Log File ---
    if args.log_file:
        log_filepath = args.log_file
    else:
        os.makedirs(DEFAULT_LOG_DIR, exist_ok=True)
        mask_str = "{}".format(int(args.masking_percentage * 100)) if args.wildtype else "100"
        log_filename = "{}_mask{}_steps{}_{}.txt".format(seq_identifier, mask_str, args.num_decoding_steps, timestamp)
        log_filepath = os.path.join(DEFAULT_LOG_DIR, log_filename)
    
    log_dirname = os.path.dirname(log_filepath)
    if log_dirname:
        os.makedirs(log_dirname, exist_ok=True)
    
    print("[INFO] Log file: {}".format(log_filepath))
    
    # Write header to log file
    header_text = (
        "======================================================================\n"
        "ESM3-Boltz Guided Generation Log\n"
        "Run Timestamp: {}\n".format(timestamp) +
        "======================================================================\n\n"
        "--- RUN PARAMETERS ---\n"
        "Target Ligand SMILES: {}\n".format(args.smiles)
    )
    
    if args.wildtype:
        header_text += (
            "Starting from: PDB {}, Chain {}\n".format(wildtype_pdb, wildtype_chain) +
            "Masking Percentage: {:.1f}%\n".format(args.masking_percentage * 100)
        )
    else:
        header_text += (
            "De novo generation\n"
            "Sequence Length: {}\n".format(args.seq_length)
        )
    
    header_text += (
        "Decoding Steps: {}\n".format(args.num_decoding_steps) +
        "Samples per Step: {}\n".format(args.num_samples_per_step) +
        "Boltz will use GPU for sequential scoring\n\n"
        "--- SEQUENCE SETUP ---\n"
        "Starting Masked Sequence:\n{}\n".format(masked_seq) +
        "======================================================================\n"
    )
    
    with open(log_filepath, 'w') as f:
        f.write(header_text)
    print(header_text)
    
    # --- 3. Initialize Model ---
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print("\nLoading ESM3 model to {}...".format(device))
    
    num_gpus = torch.cuda.device_count()
    print("Available GPUs: {}".format(num_gpus))
    
    model = ESM3.from_pretrained().float()
    
    # Use DataParallel for multi-GPU
    if num_gpus > 1:
        print("Using DataParallel across {} GPUs: {}".format(num_gpus, [torch.cuda.get_device_name(i) for i in range(num_gpus)]))
        model = torch.nn.DataParallel(model)
        model = model.to(device)
    else:
        model = model.to(device)
        if num_gpus == 1:
            print("Using single GPU: {}".format(torch.cuda.get_device_name(0)))
        else:
            print("Using CPU")
    
    # --- 4. Initialize Boltz Scorer ---
    print("\nInitializing Boltz scorer...")
    print("Cache directory: {}".format(cache_dir))
    print("Boltz will use GPU for scoring each candidate sequentially")
    
    scoring_function = BoltzScorer(
        smiles=args.smiles,
        cache_dir=cache_dir,
        timeout_sec=TIMEOUT_SEC,
        cleanup_tmp=CLEANUP_TMP,
        verbose_boltz=VERBOSE_BOLTZ,
        score_type=args.score_type
    )
    
    # --- 5. Run Guided Generation ---
    guided_decoding = ESM3GuidedDecoding(client=model, scoring_function=scoring_function)
    
    print("\n======================================================================")
    print("Starting ESM3-Boltz Guided Generation")
    print("======================================================================\n")
    
    generated_protein, all_scores, best_overall_score, best_overall_step = guided_decoding.guided_generate(
        protein=starting_protein,
        num_decoding_steps=args.num_decoding_steps,
        num_samples_per_step=args.num_samples_per_step,
        track="sequence",
        log_file_path=log_filepath
    )
    
    # --- 6. Display and Save Results ---
    if generated_protein:
        # Get best score from last step
        actual_last_step = max(all_scores.keys()) if all_scores else 0
        best_last_step_score = -float('inf')
        
        if actual_last_step > 0:
            valid_scores_last_step = [
                s for s in all_scores.get(actual_last_step, [])
                if s is not None and s != float('-inf')
            ]
            if valid_scores_last_step:
                best_last_step_score = max(valid_scores_last_step)
        
        # Write final summary
        final_log_entry = (
            "\n\n======================================================================\n"
            "   FINAL RESULTS SUMMARY\n"
            "======================================================================\n\n"
            "--- Best from Final Step (Step {}) ---\n".format(actual_last_step) +
            "Binding Affinity Score: {:.4f}\n\n".format(best_last_step_score) +
            "--- Best Overall Result ---\n"
            "Found in Step: {}\n".format(best_overall_step) +
            "Optimized Sequence:\n{}\n\n".format(generated_protein.sequence) +
            "Final Binding Affinity Score: {:.4f}\n".format(best_overall_score) +
            "======================================================================\n"
        )
        
        with open(log_filepath, 'a') as f:
            f.write(final_log_entry)
        print(final_log_entry)
        
        # Save plot
        os.makedirs(DEFAULT_RESULTS_DIR, exist_ok=True)
        plot_filename = "affinity_history_{}_{}.png".format(seq_identifier, timestamp)
        plot_filepath = os.path.join(DEFAULT_RESULTS_DIR, plot_filename)
        plot_affinity_history(all_scores, save_path=plot_filepath)
        
        print("\n✓ Generation complete!")
        print("✓ Log saved to: {}".format(log_filepath))
        print("✓ Plot saved to: {}".format(plot_filepath))
    else:
        print("[ERROR] Generation failed to produce a valid protein.")


if __name__ == "__main__":
    main()


# Example usage:
# De novo generation:
#   python main.py --smiles "CCO" --seq_length 256 --num_decoding_steps 32 --num_samples_per_step 20
#
# From wildtype:
#   python main.py --smiles "NC1=Nc2n(cnc2C(=O)N1)[C@@H]3O[C@H](CO)[C@@H](O)[C@H]3O[P](O)(O)=O" --wildtype 1RNT A --masking_percentage 0.4 --num_decoding_steps 32 --num_samples_per_step 20
