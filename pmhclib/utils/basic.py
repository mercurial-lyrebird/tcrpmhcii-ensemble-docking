
import os


def get_pdb_id_from_fp(pdb_fp: str) -> str:

    return ".".join(os.path.basename(pdb_fp).split(".")[:-1])
