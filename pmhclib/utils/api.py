
import os
from urllib import request


PDB_URL = "https://files.rcsb.org/download/"


def query_pdb(pdb_id: str, pdb_dir: str = None, pdb_fp: str = None):

    pdb_fn = pdb_id + ".pdb"
    if pdb_fp is None:
        pdb_fp = os.path.join(pdb_dir, pdb_fn)

    if not os.path.exists(pdb_fp):
        try:
            request.urlretrieve(PDB_URL + pdb_fn, pdb_fp)
        except Exception as e:
            # TODO: logger warning
            print(f"Couldn't query PDB entry {pdb_id} ({e})")

    return
