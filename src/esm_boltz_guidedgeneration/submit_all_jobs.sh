#!/bin/bash
# ====================
# Submit Multiple Guided Generation Jobs
# ====================
# This script reads protein-ligand pairs from a config file and submits
# a separate SLURM job for each pair.
#
# Usage: ./submit_all_jobs.sh [config_file]
#        Default config file: protein_ligand_pairs.conf

set -e

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${1:-${SCRIPT_DIR}/protein_ligand_pairs.conf}"
JOB_SCRIPT="${SCRIPT_DIR}/guided_generation.sh"

# Check if config file exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    echo ""
    echo "Please create a config file with the following format:"
    echo "# Comments start with #"
    echo "# Format: PDB_ID,CHAIN_ID,SMILES,MASKING_PCT,NUM_DECODING_STEPS,NUM_SAMPLES,TIME_LIMIT"
    echo "1RNT,A,NC1=Nc2n(cnc2C(=O)N1)[C@@H]3O[C@H](CO)[C@@H](O)[C@H]3O[P](O)(O)=O,0.4,32,10,15:00:00"
    echo "1DB1,A,CCO,0.3,64,20,24:00:00"
    exit 1
fi

# Check if job script exists
if [ ! -f "$JOB_SCRIPT" ]; then
    echo "ERROR: Job script not found: $JOB_SCRIPT"
    exit 1
fi

echo "=============================================="
echo "Submitting Multiple Guided Generation Jobs"
echo "=============================================="
echo "Config File: $CONFIG_FILE"
echo "Job Script: $JOB_SCRIPT"
echo ""

# Counter for jobs
JOB_COUNT=0
SUBMITTED_JOBS=()

# Read config file and submit jobs
while IFS= read -r line || [ -n "$line" ]; do
    # Skip empty lines and comments
    [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]] && continue
    
    # Parse the line (CSV format)
    IFS=',' read -r PDB_ID CHAIN_ID SMILES MASKING_PCT NUM_DECODING_STEPS NUM_SAMPLES TIME_LIMIT <<< "$line"
    
    # Trim whitespace
    PDB_ID=$(echo "$PDB_ID" | xargs)
    CHAIN_ID=$(echo "$CHAIN_ID" | xargs)
    SMILES=$(echo "$SMILES" | xargs)
    MASKING_PCT=$(echo "${MASKING_PCT:-0.4}" | xargs)
    NUM_DECODING_STEPS=$(echo "${NUM_DECODING_STEPS:-32}" | xargs)
    NUM_SAMPLES=$(echo "${NUM_SAMPLES:-10}" | xargs)
    TIME_LIMIT=$(echo "${TIME_LIMIT:-15:00:00}" | xargs)
    
    # Validate required fields
    if [ -z "$PDB_ID" ] || [ -z "$CHAIN_ID" ] || [ -z "$SMILES" ]; then
        echo "WARNING: Skipping invalid line: $line"
        continue
    fi
    
    echo "Submitting job for ${PDB_ID}_${CHAIN_ID} (time: ${TIME_LIMIT})..."
    
    # Submit job with exported variables
    JOB_ID=$(sbatch \
        --job-name="boltz_${PDB_ID}_${CHAIN_ID}" \
        --output="boltz_${PDB_ID}_${CHAIN_ID}_%j.out" \
        --error="boltz_${PDB_ID}_${CHAIN_ID}_%j.err" \
        --time="$TIME_LIMIT" \
        --export=ALL,SMILES="$SMILES",PDB_ID="$PDB_ID",CHAIN_ID="$CHAIN_ID",MASKING_PCT="$MASKING_PCT",NUM_DECODING_STEPS="$NUM_DECODING_STEPS",NUM_SAMPLES="$NUM_SAMPLES",TIME_LIMIT="$TIME_LIMIT" \
        "$JOB_SCRIPT" | awk '{print $4}')
    
    echo "  -> Submitted Job ID: $JOB_ID"
    SUBMITTED_JOBS+=("$JOB_ID:${PDB_ID}_${CHAIN_ID}")
    JOB_COUNT=$((JOB_COUNT + 1))
    
done < "$CONFIG_FILE"

echo ""
echo "=============================================="
echo "Submission Summary"
echo "=============================================="
echo "Total jobs submitted: $JOB_COUNT"
echo ""
echo "Submitted Jobs:"
for job_info in "${SUBMITTED_JOBS[@]}"; do
    IFS=':' read -r job_id job_name <<< "$job_info"
    echo "  Job $job_id: $job_name"
done
echo ""
echo "Monitor jobs with: squeue -u $USER"
echo "Cancel all jobs with: scancel ${SUBMITTED_JOBS[*]%%:*}"
