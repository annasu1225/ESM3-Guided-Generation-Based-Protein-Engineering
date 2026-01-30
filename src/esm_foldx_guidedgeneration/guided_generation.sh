#!/bin/bash
#SBATCH -N 1 
##SBATCH -C gpu&hbm80g                     
#SBATCH -C gpu 
#SBATCH -G 1 
##SBATCH -q regular 
#SBATCH -q regular              
#SBATCH -J esm3-foldx-parallel_1    
#SBATCH -t 08:00:00              
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --account=m4929_g
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=ananda@ucsd.edu

#SBATCH --array=0

#SBATCH -o foldx_run_%j.log

# PDB_FILES=("1rnt.pdb" "1hii.pdb" "1fbm.pdb" "1db1.pdb" "6q30.pdb")
# CHAIN_IDS=("A" "A" "A" "A" "A")

# PDB_FILES=("1hii.pdb" "1fbm.pdb" "1db1.pdb")
PDB_FILES=(
    # PDB 1: 1db1.pdb, Chain A - Test 3 different mask rates
    "2LZM.pdb"    # Job 0
    # "1ARR.pdb"    # Job 1
    # "1ARR.pdb"    # Job 2
    # "1ARR.pdb"  
    # "1CSP.pdb"
    # "1CSP.pdb"
    # "1CSP.pdb"
    # "1CSP.pdb"
    # "1CSP.pdb"
    # # PDB 2: 6q30.pdb, Chain A - Test 3 different mask rates
    # "1K9Q.pdb"    # Job 3
    # "1K9Q.pdb"    # Job 4
    # "1K9Q.pdb" 
    # "1K9Q.pdb"    # Job 5
    
    # # PDB 3: 2lwy.pdb, Chain A - Test 3 different mask rates
    # "2lwy.pdb"    # Job 6
    # "2lwy.pdb"    # Job 7
    # "2lwy.pdb"    # Job 8
    
    # # PDB 4: 1hii.pdb, Chain B - Test 3 different mask rates (different chain example)
    # "1hii.pdb"    # Job 9
    # "1hii.pdb"    # Job 10
    # "1hii.pdb"    # Job 11
)

CHAIN_IDS=(
    # Corresponding chain IDs
    "A"   # Jobs 0-2: 1db1.pdb chain A
    # "A" "A" "A" "A" "A"   # Jobs 3-5: 6q30.pdb chain A
    # "A" "A" "A"    # Jobs 6-8: 2lwy.pdb chain A
    # "B" "B" "B"    # Jobs 9-11: 1hii.pdb chain B
)

MASK_RATES=(
    # Different mask rates for 1db1.pdb chain A
    0.25          # Job 0 - Conservative
    # 0.40          # Job 1 - Moderate
    0.30
    # 0.20          # Job 2 - Aggressive
    
    # # Different mask rates for 6q30.pdb chain A
    # 0.60
    # 0.50          # Job 0 - Conservative
    # 0.40          # Job 1 - Moderate
    # 0.30
    # 0.20         # Job 5
    
    # # Different mask rates for 2lwy.pdb chain A
    # 0.20          # Job 6
    # 0.25          # Job 7
    # 0.30          # Job 8
    
    # # Different mask rates for 1hii.pdb chain B
    # 0.25          # Job 9
    # 0.30          # Job 10
    # 0.35          # Job 11
)

DECODING_STEPS=(
    # 1db1.pdb (Large protein ~250aa)
    48            # Job 0 - More steps for lower mask
    48            # Job 1
    # 32
    # # 32            # Job 2
    
    # # 6q30.pdb (Large protein)
    # 32            # Job 0 - More steps for lower mask
    # 32            # Job 1
    # 32
    # 32  
    # 32         # Job 5
    
    # # 2lwy.pdb (Medium protein ~137aa)
    # 32            # Job 6
    # 32            # Job 7
    # 28            # Job 8
    
    # # 1hii.pdb (Small protein ~99aa)
    # 28            # Job 9
    # 24            # Job 10
    # 24            # Job 11
)

SAMPLES_PER_STEP=(
    20     # Jobs 0-2
    # 20 20 20 20 20     # Jobs 3-5
    # 20 20 20      # Jobs 6-8
    # 15 15 15      # Jobs 9-11 (fewer for smaller protein)
)

NUM_WORKERS=(
    20      # Jobs 0-2
    # 20 20 20 20 20     # Jobs 3-5
    # 20 20 20      # Jobs 6-8
    # 15 15 15      # Jobs 9-11
)

# MASK_RATES=(0.3 0.3 0.4 0.20 0.20)
# MASK_RATES=(0.2 0.2)

export SLURM_CPU_BIND=cores

CURRENT_PDB=${PDB_FILES[$SLURM_ARRAY_TASK_ID]}
CURRENT_CHAIN=${CHAIN_IDS[$SLURM_ARRAY_TASK_ID]}
CURRENT_MASK=${MASK_RATES[$SLURM_ARRAY_TASK_ID]}

CURRENT_STEPS=${DECODING_STEPS[$SLURM_ARRAY_TASK_ID]}
CURRENT_SAMPLES=${SAMPLES_PER_STEP[$SLURM_ARRAY_TASK_ID]}
CURRENT_WORKERS=${NUM_WORKERS[$SLURM_ARRAY_TASK_ID]}

echo "========================================================"
echo "STARTING JOB: ${SLURM_JOB_ID}, ARRAY_TASK: ${SLURM_ARRAY_TASK_ID}"
echo "Protein: ${CURRENT_PDB}, Chain: ${CURRENT_CHAIN}, Mask: ${CURRENT_MASK}"
echo "Steps: ${CURRENT_STEPS}, Samples: ${CURRENT_SAMPLES}, Workers: ${CURRENT_WORKERS}"
echo "========================================================"

module load conda

conda activate /pscratch/sd/a/ananda/conda_envs/proteinenv_new

cd /pscratch/sd/a/ananda/ESM3-Guided-Generation-Based-Protein-Engineering/src


srun python -m esm_foldx_guidedgeneration.main --pdb_filename "$CURRENT_PDB" --chain_id "$CURRENT_CHAIN" --masking_percentage "$CURRENT_MASK" --num_decoding_steps "$CURRENT_STEPS" --num_samples_per_step "$CURRENT_SAMPLES" --num_workers "$CURRENT_WORKERS"

echo "========================================================"
echo "FINISHED JOB: ${SLURM_JOB_ID}, ARRAY_TASK: ${SLURM_ARRAY_TASK_ID}"
echo "========================================================"


