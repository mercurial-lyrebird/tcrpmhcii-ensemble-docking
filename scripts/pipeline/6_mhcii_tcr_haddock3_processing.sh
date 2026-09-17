#!/bin/bash

#SBATCH --partition=commons
#SBATCH --time=0:30:00
#SBATCH --job-name=6_mhcii_tcr_haddock3_processing
#SBATCH --mail-user=ab215@rice.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-49

PDB_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" pdb_ids.txt)

for SETTING in "ensemble/annealing"; do
    gunzip ../../haddock3_cdr3/${PDB_ID}/$SETTING/results/1_rigidbody/rigidbody_*.pdb.gz
    python3 6_mhcii_tcr_haddock3_processing.py $PDB_ID ../../haddock3_cdr3/$PDB_ID/$SETTING/results/1_rigidbody/ ../../data/mhcii_tcr_ensembles/$PDB_ID/$SETTING/haddock3_rigidbody_unrestrained/
    gzip ../../haddock3_cdr3/${PDB_ID}/$SETTING/results/1_rigidbody/rigidbody_*.pdb
done