# ESM3-Boltz Guided Generation - Usage Notes

This directory contains scripts for running ESM3 guided protein generation with Boltz-1 binding affinity scoring on HPC clusters using SLURM.

## Overview

The guided generation process:
1. Takes a protein (from PDB) and a target ligand (SMILES)
2. Masks a percentage of the protein sequence
3. Iteratively unmasks positions using ESM3
4. Scores each candidate using Boltz-1 for ligand binding affinity
5. Selects the best candidate at each step

---

## Scripts

### `guided_generation.sh`
The main SLURM job script that runs a single guided generation job.

**Parameters (passed via environment variables):**
| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `SMILES` | Yes | - | SMILES string of the target ligand |
| `PDB_ID` | Yes | - | PDB identifier (e.g., `1RNT`) |
| `CHAIN_ID` | Yes | - | Chain identifier (e.g., `A`) |
| `MASKING_PCT` | No | `0.4` | Fraction of residues to mask (0.0-1.0) |
| `NUM_DECODING_STEPS` | No | `32` | Number of generation steps |
| `NUM_SAMPLES` | No | `10` | Candidates generated per step |
| `TIME_LIMIT` | No | `15:00:00` | SLURM time limit (HH:MM:SS) |

---

### `submit_all_jobs.sh`
Batch submission script that reads protein-ligand pairs from a config file and submits multiple jobs.

**Usage:**
```bash
./submit_all_jobs.sh [config_file]
```

Default config file: `protein_ligand_pairs.conf`

---

### `protein_ligand_pairs.conf`
Configuration file defining protein-ligand pairs for batch submission.

**Format:**
```
PDB_ID,CHAIN_ID,SMILES,MASKING_PCT,NUM_DECODING_STEPS,NUM_SAMPLES,TIME_LIMIT
```

**Example:**
```conf
# Comments start with #
1RNT,A,NC1=Nc2n(cnc2C(=O)N1)[C@@H]3O[C@H](CO)[C@@H](O)[C@H]3O[P](O)(O)=O,0.4,32,10,12:00:00
1HII,A,CC(C)[C@H](NC(C)=O)C(=O)...,0.3,32,10,10:00:00
1DB1,A,C[C@H](CCCC(C)(C)O)...,0.3,64,15,20:00:00
```

---

## Usage Examples

### Submitting Multiple Jobs (Recommended)

1. Edit `protein_ligand_pairs.conf` with your protein-ligand pairs
2. Run:
```bash
cd /home/as4272/ESM3-Guided-Generation-Based-Protein-Engineering/src/esm_boltz_guidedgeneration
./submit_all_jobs.sh
```

Or with a custom config file:
```bash
./submit_all_jobs.sh my_custom_pairs.conf
```

---

### Submitting a Single Job

Use `sbatch` with `--export` to pass parameters directly:

```bash
sbatch \
  --job-name="boltz_1HII_A" \
  --output="boltz_1HII_A_%j.out" \
  --error="boltz_1HII_A_%j.err" \
  --time="10:00:00" \
  --export=ALL,SMILES="CC(C)[C@H](NC(C)=O)C(=O)N[C@@H](Cc1ccccc1)[C@@H](O)CN(CC2CCCCC2)NC(=O)[C@@H](NC(C)=O)C(C)C",PDB_ID="1HII",CHAIN_ID="A",MASKING_PCT="0.3",NUM_DECODING_STEPS="32",NUM_SAMPLES="10",TIME_LIMIT="10:00:00" \
  guided_generation.sh
```

**Another example (1RNT with GMP):**
```bash
sbatch \
  --job-name="boltz_1RNT_A" \
  --output="boltz_1RNT_A_%j.out" \
  --error="boltz_1RNT_A_%j.err" \
  --time="15:00:00" \
  --export=ALL,SMILES="NC1=Nc2n(cnc2C(=O)N1)[C@@H]3O[C@H](CO)[C@@H](O)[C@H]3O[P](O)(O)=O",PDB_ID="1RNT",CHAIN_ID="A",MASKING_PCT="0.4",NUM_DECODING_STEPS="32",NUM_SAMPLES="10" \
  guided_generation.sh
```

---

## Monitoring Jobs

**Check job status:**
```bash
squeue -u $USER
```

**View job output in real-time:**
```bash
tail -f boltz_1RNT_A_<job_id>.out
```

**Cancel a job:**
```bash
scancel <job_id>
```

**Cancel all your jobs:**
```bash
scancel -u $USER
```

---

## Output Files

| File | Description |
|------|-------------|
| `boltz_<PDB>_<CHAIN>_<jobid>.out` | Standard output (progress, results) |
| `boltz_<PDB>_<CHAIN>_<jobid>.err` | Standard error (warnings, errors) |
| `logs/boltz/<PDB>_<CHAIN>_mask<X>_steps<Y>_<timestamp>.txt` | Detailed generation log |
| `logs/gpu_usage/gpu_usage_<jobid>_<timestamp>.csv` | GPU utilization metrics |
| `results/boltz/affinity_history_<PDB>_<CHAIN>_<timestamp>.png` | Score progression plot |

---

## Troubleshooting

### ESM API Errors
If a job fails with `ESMProteinError`, the error message will be printed in the `.err` file. Common causes:
- API rate limiting
- Sequence validation issues
- Biosecurity screening (see below)

### Biosecurity Screening - Sequences of Concern

The ESM3 Forge API performs biosecurity screening and will **reject sequences flagged as potential biosecurity concerns**. This includes:
- **Viral proteins** (e.g., HIV-1 protease - PDB `1HII`)
- Toxins and other dangerous pathogens
- Select agents

**Error message:**
```
RuntimeError: ESM API error: Failure in forward_and_sample: {"status":"error","message":"We detected a potential sequence of concern. Please apply for elevated access and once access is approved you can self-disclose with potential_sequence_of_concern=True"}
```

**Solutions:**
1. **Use a different protein** that isn't flagged as a biosecurity concern
2. **Apply for elevated access** at [EvolutionaryScale Forge](https://forge.evolutionaryscale.ai) - once approved, you can add `potential_sequence_of_concern=True` to API calls

**Note:** This is an intentional safety measure by EvolutionaryScale, not a bug in the scripts.

### Job Timeout
Increase `TIME_LIMIT` in the config file or `--time` flag. Larger proteins or more decoding steps require more time.

### GPU Memory Issues
The default allocation requests 1 GPU with 80GB memory (A100). If running on smaller GPUs, reduce `NUM_SAMPLES` per step.

---

## File Structure

```
esm_boltz_guidedgeneration/
├── guided_generation.sh      # SLURM job script
├── submit_all_jobs.sh        # Batch submission script
├── protein_ligand_pairs.conf # Config file for batch jobs
├── main.py                   # Main Python entry point
├── guided_generation.py      # Core generation logic
├── boltz_scoring_utils.py    # Boltz-1 scoring utilities
├── boltz_cache/              # Cached Boltz results
└── notes.md                  # This file
```
