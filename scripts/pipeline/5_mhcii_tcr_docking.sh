#!/bin/bash

#SBATCH --partition=commons
#SBATCH --time=8:59:59
#SBATCH --mem=24G
#SBATCH --ntasks=16
#SBATCH --job-name=5_mhcii_tcr_docking
#SBATCH --mail-user=ab215@rice.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-49

PDB_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" pdb_ids.txt)
SETTING=static/crystal
# SETTING=static/model
# SETTING=ensemble/annealing

cd ../../haddock3_cdr3/${PDB_ID}/${SETTING}/
haddock3 dock.cfg