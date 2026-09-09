
import os
import shutil
import sys

if ".." not in sys.path:
    sys.path.append("..")

import pandas as pd

from pmhclib.components.pmhc import TCRpMHCII
from pmhclib.modeling import run_pandora


if __name__ == "__main__":

    pdb_id = sys.argv[1].upper()

    template_df = pd.read_csv("../../data/mhcii_tcr_templates.csv")
    template = template_df[template_df.pdb_id == pdb_id].iloc[0]
    allotypes = template.allotype.split(";")
    anchors = [int(x) for x in template.anchors.split(";")]

    pdb_fp = os.path.join(f"../../data/mhcii_tcr_templates/{pdb_id}.pdb")
    pmhc = TCRpMHCII(pdb_fp, mhc_allotype=allotypes, anchors=anchors).pmhc

    run_id = f"{pdb_id}"
    if os.path.exists(run_id):
        shutil.rmtree(run_id)

    kwargs = {
        "n_loop_models": 100,
        "benchmark": True
    }
    run_pandora(run_id, pmhc, allotypes, anchors, **kwargs)

    dest_dir = f"../../data/mhcii_ensembles/{pdb_id}/pandora2/"
    if os.path.exists(dest_dir):
        shutil.rmtree(dest_dir)
    else:
        os.makedirs(dest_dir)

    for item in os.listdir(run_id):
        shutil.move(os.path.join(run_id, item), os.path.join(dest_dir, item))

    shutil.rmtree(run_id)
