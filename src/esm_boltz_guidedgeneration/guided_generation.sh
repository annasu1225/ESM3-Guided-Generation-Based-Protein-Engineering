#!/bin/bash

#SBATCH --partition=gpu
#SBATCH --job-name=boltz_guided_gen
#SBATCH -t 15:00:00
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH --output=boltz_guided_generation_%j.out
#SBATCH --error=boltz_guided_generation_%j.err
#SBATCH --gres=gpu:1
#SBATCH --mail-user=4752279178@vtext.com
#SBATCH --mail-type=ALL

# ====================
# Environment Setup
# ====================
export ESM_API_TOKEN="6Zk4FIijMlhj5iMQqWUU2q"

module load CUDA/11.8.0
module load miniconda

# Source conda for non-interactive shell before activating
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate boltz_esm

# ====================
# Directory Setup
# ====================
# Use SLURM_SUBMIT_DIR (where sbatch was called from)
SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
PROJECT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"
LOG_DIR="${PROJECT_DIR}/logs/gpu_usage"
mkdir -p "$LOG_DIR"

# ====================
# GPU Logging Setup
# ====================
GPU_LOG="${LOG_DIR}/gpu_usage_${SLURM_JOB_ID}_$(date +%Y%m%d_%H%M%S).csv"

# Print job info
echo "=============================================="
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start Time: $(date)"
echo "Working Directory: $SCRIPT_DIR"
echo "GPU Log: $GPU_LOG"
echo "=============================================="

# Print GPU info
echo "GPU Information:"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
echo "=============================================="

# Start GPU monitoring in background
# Logs: timestamp, GPU utilization %, memory used MB, memory total MB, temperature
echo "timestamp,gpu_util_percent,memory_used_mb,memory_total_mb,temperature_c,power_draw_w" > "$GPU_LOG"

nvidia-smi --query-gpu=timestamp,utilization.gpu,memory.used,memory.total,temperature.gpu,power.draw \
    --format=csv,noheader,nounits -l 30 >> "$GPU_LOG" &
GPU_MONITOR_PID=$!

echo "Started GPU monitoring (PID: $GPU_MONITOR_PID)"

# ====================
# Main Execution
# ====================
echo ""
echo "Starting Guided Generation..."
echo ""

python "${SCRIPT_DIR}/main.py" \
    --smiles "NC1=Nc2n(cnc2C(=O)N1)[C@@H]3O[C@H](CO)[C@@H](O)[C@H]3O[P](O)(O)=O" \
    --wildtype 1RNT A \
    --masking_percentage 0.4 \
    --num_decoding_steps 32 \
    --num_samples_per_step 10

PYTHON_EXIT_CODE=$?

# ====================
# Cleanup & Summary
# ====================
# Stop GPU monitoring
kill $GPU_MONITOR_PID 2>/dev/null
wait $GPU_MONITOR_PID 2>/dev/null

echo ""
echo "=============================================="
echo "Job Complete"
echo "End Time: $(date)"
echo "Python Exit Code: $PYTHON_EXIT_CODE"
echo "=============================================="

# Generate GPU usage summary
if [ -f "$GPU_LOG" ]; then
    echo ""
    echo "GPU Usage Summary:"
    echo "------------------"
    # Skip header, calculate stats
    tail -n +2 "$GPU_LOG" | awk -F',' '
    BEGIN { 
        count=0; sum_util=0; sum_mem=0; max_util=0; max_mem=0 
    }
    {
        count++
        sum_util += $2
        sum_mem += $3
        if ($2 > max_util) max_util = $2
        if ($3 > max_mem) max_mem = $3
    }
    END {
        if (count > 0) {
            printf "  Samples: %d\n", count
            printf "  Avg GPU Utilization: %.1f%%\n", sum_util/count
            printf "  Max GPU Utilization: %.1f%%\n", max_util
            printf "  Avg Memory Used: %.0f MB\n", sum_mem/count
            printf "  Max Memory Used: %.0f MB\n", max_mem
        }
    }'
    echo "  Full log: $GPU_LOG"
fi

exit $PYTHON_EXIT_CODE
