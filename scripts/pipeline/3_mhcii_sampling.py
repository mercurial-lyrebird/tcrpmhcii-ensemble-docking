
import os
import sys

import pandas as pd

from pmhclib.components.pmhc import pMHCII
from pmhclib.modeling.sampling import sample, SimulatedAnnealingSampler


if __name__ == "__main__":

    pdb_id = sys.argv[1]
    init_model_dir = sys.argv[2]

    init_model = pMHCII(fp=os.path.join(init_model_dir, f"{pdb_id}.pdb"))

    df = pd.read_csv("../../data/mhcii_tcr_templates.csv", index_col="pdb_id")
    record = df.loc[pdb_id]
    anchors = [int(x) for x in record.anchors.split(";")]

    sample(
        100,
        SimulatedAnnealingSampler,
        [init_model],
        {
            "chain_id": "P",
            "region_lens": (anchors[0] - 1, len(record.peptide_seq) - anchors[-1]),
            "n_iter": 1000
        },
        os.path.join(
            "../../data/mhcii_ensembles/", pdb_id, "annealing"
        ),
        f"{pdb_id}__annealing", n_proc=4
    )
