
from glob import glob
import os

from pymol import cmd


CHAIN_DICT = {
    "1D9K": ["A", "B"],
    "1YMM": ["D", "E"],
    "1ZGL": ["D", "E"],
    "2IAN": ["D", "E"],
    "2Z31": ["A", "B"],
    "3MBE": ["C", "D"],
    "3QIW": ["C", "D"],
    "4GG6": ["E", "F"],
    "4H1L": ["I", "J"],
    "4OZG": ["E", "F"],
    "4OZH": ["E", "F"],
    "4OZI": ["E", "F"],
    "4P2Q": ["D", "E"],
    "4P2R": ["D", "E"],
    "6DFS": ["A", "B"],
    "6PY2": ["D", "E"]
}


if __name__ == "__main__":

    for pdb_fp in glob("./full_templates/*.pdb"):
        pdb_fn = os.path.basename(pdb_fp)
        pdb_id = pdb_fn[:4]
        cmd.load(pdb_fp, pdb_id)
        cmd.alter(f"({pdb_id} & chain {CHAIN_DICT[pdb_id][0]})", "chain='A'")
        cmd.alter(f"({pdb_id} & chain {CHAIN_DICT[pdb_id][1]})", "chain='B'")
        print(cmd.get_chains(pdb_id))
        cmd.save(os.path.join("./data/mhcii_tcr_templates/", pdb_fn))
        cmd.delete(pdb_id)
