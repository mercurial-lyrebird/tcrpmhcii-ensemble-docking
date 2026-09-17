
from glob import glob
from itertools import product
from multiprocessing import Pool
import os
import sys
from typing import Tuple

import pandas as pd

from pmhclib.components.pmhc import TCRpMHCII
from pmhclib.modeling import Ensemble
from pmhclib.modeling.scoring import RMSDScorer, DockQScorer, RosettaScorer


DOCKING_METHODS = ["haddock3_rigidbody", "haddock3_flexref"]


def process_ensemble(x: Tuple[str, str, str]):

    (pdb_id, sampling_method, docking_method) = x
    ensemble_dir = os.path.join("../../data/mhcii_tcr_ensembles/", os.path.join(*x))

    template_df = pd.read_csv(
        "../../data/mhcii_tcr_templates.csv", index_col="pdb_id"
    )
    record = template_df.loc[pdb_id]
    anchors = (int(record.anch_1), int(record.anch_4))
    template = TCRpMHCII(
        os.path.join(f"../../data/mhcii_tcr_templates/{pdb_id}.pdb"),
        anchors=anchors
    )
    ensemble_id = "__".join([pdb_id, sampling_method.replace("/", "_"), docking_method])
    ensemble = Ensemble.from_pdb_dir(ensemble_id, ensemble_dir)

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
        dockq_scorer = DockQScorer(template, ["P", ("A", "B")])
        stats_df["dockq"] = ensemble.score(dockq_scorer)
        added += 1

    if "dockq_fnat" not in stats_df.columns:
        dockq_scorer = DockQScorer(template, ["P", ("A", "B")], f_nat=True)
        stats_df["dockq_fnat"] = ensemble.score(dockq_scorer)
        added += 1

    def compute_docking_angles(record):

        complex = TCRpMHCII(record.pdb_fp, anchors=template.pmhc.anchors)
        angles = complex.compute_docking_angles()
        complex.delete(rec=True)

        return angles

    if "crossing_angle" not in stats_df.columns or "incident_angle" not in stats_df.columns:
        angles = stats_df.apply(compute_docking_angles, axis=1)
        stats_df["crossing_angle"] = angles.apply(lambda x: x[0])
        stats_df["incident_angle"] = angles.apply(lambda x: x[1])
        added += 2

    if "rosetta" not in stats_df.columns:
        rosetta_scorer = RosettaScorer()
        stats_df["rosetta"] = ensemble.score(rosetta_scorer)
        added += 1

    def compute_binding_core_sasa(record):

        complex = TCRpMHCII(
            record.pdb_fp, id=record.name, anchors=template.pmhc.anchors
        )
        sasa = complex.peptide_sasa(core_only=True).mean()
        complex.delete(rec=True)

        return sasa

    if "binding_core_sasa" not in stats_df.columns:
        stats_df["binding_core_sasa"] = stats_df.apply(
            compute_binding_core_sasa, axis=1
        )
        added += 1

    if added > 0:
        print(f"Added {added} to {ensemble_id}")
        stats_df.to_csv(stats_csv_fp)
    else:
        print(f"Skipped {ensemble_id}")

    return


if __name__ == "__main__":

    pdb_id = sys.argv[1]
    sampling_method = sys.argv[2]
    docking_method = sys.argv[3]
    process_ensemble((pdb_id, sampling_method, docking_method))
