#!/bin/bash

PDB_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" pdb_ids.txt)
SETTING="ensemble/annealing"
STAGE="haddock3_rigidbody"

python3 7_mhcii_tcr_ensemble_scoring.py $PDB_ID $SETTING $STAGE
