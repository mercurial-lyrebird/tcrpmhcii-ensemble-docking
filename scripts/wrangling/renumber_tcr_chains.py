
from glob import glob
import os

from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Structure import Structure


def renumber(structure: Structure, chain_id: str) -> Structure:

    model = structure[0]
    for chain in model:
        if chain.id == chain_id:
            for i, res in enumerate(chain.get_residues(), start=1):
                old_res_id = res.id
                res.id = (' ', i, ' ')
                print(f"{old_res_id} -> {res.id}")
            return structure

    return


if __name__ == "__main__":

    parser = PDBParser(QUIET=True)
    io = PDBIO()

    for pdb_fp in glob("./data/mhcii_tcr_templates/*.pdb"):
        pdb_fn = os.path.basename(pdb_fp)
        pdb_id = pdb_fn[:4]
        structure = parser.get_structure(pdb_id, pdb_fp)
        structure = renumber(structure, "A")
        structure = renumber(structure, "B")
        io.set_structure(structure)
        io.save(pdb_fp)