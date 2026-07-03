
import os
from pathlib import Path


DEFAULT_STORE = os.path.join(str(Path.home()), "/scratch/pmhclib/")
DEFAULT_PDB_STORE = os.path.join(DEFAULT_STORE, "pdb/")
DEFAULT_TRAJ_STORE = os.path.join(DEFAULT_STORE, "traj/")
DEFAULT_FASTA_STORE = os.path.join(DEFAULT_STORE, "fasta/")

IMGTHLA_STORE = os.path.join(str(Path.home()), "/scratch/imgthla/")
IMGTHLA_FASTA_STORE = os.path.join(IMGTHLA_STORE, "fasta/")
