
import os
import subprocess
from typing import List

from ..components.pmhc import MHCIPeptide, MHCIIPeptide
from ..components.pmhc.structure.pmhc import pMHCI, pMHCII
from ..utils import get_pdb_id_from_fp
from .eval import rmsd


def run_pepglad(
    mhc_pdb_fp: str, mhc_class: str, allotypes: List[str],
    receptor_chain_ids: List[str], peptide_chain_id: str, pepglad_src_dir: str
):

    init_wd = os.getcwd()
    os.chdir(pepglad_src_dir)

    pmhc_class = pMHCII if mhc_class == "II" else pMHCI
    pmhc = pmhc_class(mhc_pdb_fp, allotype=None)
    peptide_class = MHCIIPeptide if mhc_class == "II" else MHCIPeptide
    peptide = peptide_class.from_pmhc(pmhc, chain_id=peptide_chain_id)
    peptide.infer_anchors(allotypes)

    subprocess.call(
        [
            "python", "-m", "api.detect_pocket",
            "--pdb", mhc_pdb_fp,
            "--target_chains", *receptor_chain_ids,
            "--ligand_chains", peptide_chain_id,
            "--out", "pocket.json"
        ]
    )

    output_dir = "./output/struct_pred"
    subprocess.call(
        [
            # "CUDA_VISIBLE_DEVICES=0",
            "python", "-m", "api.run",
            "--mode", "struct_pred",
            "--pdb", mhc_pdb_fp,
            "--pocket", "pocket.json",
            "--out_dir", output_dir,
            "--peptide_seq", str(peptide.sequence),
            "--n_samples", str(10)
        ]
    )

    pmhc_model = pMHC(
        f"./output/struct_pred/{get_pdb_id_from_fp(mhc_pdb_fp)}_0.pdb",
        allotype=None
    )
    peptide_model = peptide_class.from_pmhc(pmhc_model, chain_id="O")

    rmsd(peptide, peptide_model, reindex_model=True)

    os.chdir(init_wd)

    return
