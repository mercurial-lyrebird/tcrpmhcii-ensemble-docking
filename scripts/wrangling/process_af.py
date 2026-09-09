
from glob import glob
import os
from pymol import cmd
import sys


if __name__ == "__main__":

    ensemble_dir = sys.argv[1]

    pdb_id = ensemble_dir[-4:]
    for fp in glob(os.path.join(ensemble_dir, "*.cif")):
        id = os.path.basename(fp)[:-4]
        cmd.load(fp, id)
        if cmd.get_chains(id) == ["A", "C", "D", "E"]:
            cmd.delete(id)
            continue
        cmd.alter(f"({id} & chain D)", "chain='E'")
        cmd.alter(f"({id} & chain C)", "chain='D'")
        cmd.alter(f"({id} & chain B)", "chain='C'")
        cmd.save(fp.replace(".cif", ".pdb"), id)
        cmd.delete(id)