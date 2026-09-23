#!/bin/bash

PDB_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" pdb_ids.txt)
SETTING=ensemble/annealing

cd ../../haddock3_cdr3/${PDB_ID}/${SETTING}/
rm -rf results/
haddock3 dock.cfg
