#!/bin/bash

#SBATCH --partition=commons
#SBATCH --time=0:30:00
#SBATCH --mem=8G
#SBATCH --job-name=7_mhcii_tcr_ensemble_scoring
#SBATCH --mail-user=ab215@rice.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-49

PDB_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" pdb_ids.txt)

for SETTING in "ensemble/annealing"; do
    for STAGE in "haddock3_rigidbody"; do
        python3 7_mhcii_tcr_ensemble_scoring.py $PDB_ID $SETTING $STAGE
    done
done