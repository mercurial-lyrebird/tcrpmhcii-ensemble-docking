#!/bin/bash

#SBATCH --partition=commons
#SBATCH --time=0:05:00
#SBATCH --job-name=4_mhcii_tcr_docking_setup
#SBATCH --mail-user=ab215@rice.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-49

PDB_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" pdb_ids.txt)

# python3 4_mhcii_tcr_docking_setup.py $PDB_ID
# python3 4_mhcii_tcr_docking_setup.py $PDB_ID -c ../../data/mhcii_models/pandora2/${PDB_ID}.pdb
python3 4_mhcii_tcr_docking_setup.py $PDB_ID -s annealing