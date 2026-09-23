#!/bin/bash

PDB_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" pdb_ids.txt)
SETTING="ensemble/annealing"

gunzip ../../haddock3_cdr3/${PDB_ID}/$SETTING/results/1_rigidbody/rigidbody_*.pdb.gz
python3 6_mhcii_tcr_haddock3_processing.py $PDB_ID ../../haddock3_cdr3/$PDB_ID/$SETTING/results/1_rigidbody/ ../../data/mhcii_tcr_ensembles/$PDB_ID/$SETTING/haddock3_rigidbody/
gzip ../../haddock3_cdr3/${PDB_ID}/$SETTING/results/1_rigidbody/rigidbody_*.pdb
