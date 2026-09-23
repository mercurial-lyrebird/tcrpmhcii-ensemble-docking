#!/bin/bash

PDB_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" pdb_ids.txt)

python3 4_mhcii_tcr_docking_setup.py $PDB_ID -s annealing
