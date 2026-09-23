
from itertools import product
from multiprocessing import Pool
import os
import sys
from typing import Tuple

import pandas as pd

from pmhclib.components.pmhc import pMHCII
from pmhclib.modeling import Ensemble
from pmhclib.modeling.scoring import RMSDScorer, DockQScorer, RosettaScorer


SAMPLING_METHODS = ["amber", "annealing", "pandora2"]
FN_PATTERNS = {
    "amber": "*_frame_*.pdb",
    "annealing": "conf_*.pdb",
    "pandora2": "*BL*.pdb"
}


def process_ensemble(x: Tuple[str, str]):

    (pdb_id, sampling_method) = x
    ensemble_dir = os.path.join("../../data/mhcii_ensembles/", os.path.join(*x))

    template_df = pd.read_csv(
        "../../data/mhcii_tcr_templates.csv", index_col="pdb_id"
    )
    record = template_df.loc[pdb_id]
    anchors = (int(record.anch_1), int(record.anch_4))
    template = pMHCII(
        os.path.join(f"../../data/mhcii_templates/{pdb_id}.pdb"),
        anchors=anchors
    )
    ensemble_id = "__".join([pdb_id, sampling_method])
    ensemble = Ensemble.from_pdb_dir(
        ensemble_id, ensemble_dir, fn_pattern=FN_PATTERNS[sampling_method]
    )

    if not os.path.exists(ensemble_dir):
        print(f"{ensemble_dir} not found")
        return

    print(f"Working on {ensemble_id}")

    stats_csv_fp = os.path.join(ensemble_dir, "stats.csv")
    if os.path.exists(stats_csv_fp):
        stats_df = pd.read_csv(stats_csv_fp, index_col="model_id")
    else:
        stats_df = pd.DataFrame(index=ensemble.ids)
        stats_df.index.name = "model_id"

    stats_df["pdb_fp"] = ensemble.pdb_fps

    added = 0

    if "rmsd" not in stats_df.columns:
        rmsd_scorer = RMSDScorer(template)
        stats_df["rmsd"] = ensemble.score(rmsd_scorer)
        added += 1

    if "dockq" not in stats_df.columns:
        dockq_scorer = DockQScorer(template, ["P", ("M", "N")])
        stats_df["dockq"] = ensemble.score(dockq_scorer)
        added += 1

    if "dockq_fnat" not in stats_df.columns:
        dockq_scorer = DockQScorer(template, ["P", ("M", "N")], f_nat=True)
        stats_df["dockq_fnat"] = ensemble.score(dockq_scorer)
        added += 1

    if "rosetta" not in stats_df.columns:
        rosetta_scorer = RosettaScorer()
        stats_df["rosetta"] = ensemble.score(rosetta_scorer)
        added += 1

    if added > 0:
        print(f"Added {added} to {ensemble_id}")
        stats_df.to_csv(stats_csv_fp)
    else:
        print(f"Skipped {ensemble_id}")

    return


if __name__ == "__main__":

    pdb_id = sys.argv[1]
    cases = list(product([pdb_id], SAMPLING_METHODS))
    print(cases)

    with Pool(processes=8) as pool:
        pool.map(process_ensemble, cases)
