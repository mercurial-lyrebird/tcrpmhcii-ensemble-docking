
import math
import os
import sys
import time
from tqdm import tqdm

import pandas as pd

import pyrosetta
pyrosetta.init("-out:levels basic:error core:error")

from pyrosetta import pose_from_pdb, FoldTree, get_fa_scorefxn
from pyrosetta.rosetta.core.pose import deep_copy
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.protocols.simple_moves import SmallMover
from pyrosetta.rosetta.protocols.moves import MonteCarlo
from pyrosetta.rosetta.protocols.minimization_packing import MinMover

from pmhclib.utils import get_pdb_id_from_fp
from pmhclib.components.pmhc import pMHC


if __name__ == "__main__":

    input_arg = sys.argv[1]
    if input_arg.endswith(".pdb"):
        pdb_fp = input_arg
        pdb_id = get_pdb_id_from_fp(pdb_fp)[:4]
        job_id = ".".join(os.path.basename(pdb_fp).split(".")[:-1])
    else:
        pdb_id = input_arg
        pdb_fp = f"./data/mhcii_templates/{pdb_id}.pdb"
        job_id = pdb_id

    pmhc = pMHC(pdb_fp, allotype=None)

    print(pdb_id)
    template_df = pd.read_csv("./data/mhcii_templates.csv")
    record = template_df[template_df.pdb_id == pdb_id].iloc[0]
    anchors = [int(i) for i in record.anchors.split(";")]
    (left_anchor, right_anchor) = anchors[0], anchors[-1]
    print(left_anchor, right_anchor)

    crystal_pose = pose_from_pdb(pmhc.fp)
    res_chain_ids = [crystal_pose.chain(i + 1) for i in range(len(crystal_pose.residues))]
    chain_ids = set(res_chain_ids)

    peptide_chain_id = None
    min_chain_len = math.inf
    for chain_id in chain_ids:
        chain_len = res_chain_ids.count(chain_id)
        if chain_len < min_chain_len:
            min_chain_len = chain_len
            peptide_chain_id = chain_id

    print(f"peptide: chain {peptide_chain_id} with {min_chain_len} residues")
    peptide_start_i = res_chain_ids.index(peptide_chain_id)
    peptide_len = res_chain_ids.count(peptide_chain_id)

    peptide_start_res = peptide_start_i + 1
    peptide_left_flank_end_res = peptide_start_i + left_anchor - 1
    peptide_left_anchor_res = peptide_left_flank_end_res + 1
    peptide_right_flank_start_res = peptide_start_i + right_anchor + 1
    peptide_right_anchor_res = peptide_right_flank_start_res - 1
    peptide_end_res = peptide_start_i + peptide_len

    print(
        peptide_start_res, peptide_left_flank_end_res, peptide_left_anchor_res,
        peptide_right_anchor_res, peptide_right_flank_start_res, peptide_end_res
    )

    init_pose = deep_copy(crystal_pose)

    left_ft = FoldTree()
    left_ft.add_edge(1, 179, -1)
    left_ft.add_edge(179, 180, 1)
    left_ft.add_edge(180, peptide_start_res -1, -1)
    left_ft.add_edge(peptide_start_res - 1, peptide_end_res, 2)
    left_ft.add_edge(peptide_end_res, peptide_start_res, -1)

    init_pose.fold_tree(left_ft)

    for i in range(peptide_left_flank_end_res, peptide_start_res - 1, -1):
        init_pose.set_phi(i, 180.0)
        init_pose.set_psi(i, 180.0)

    right_ft = FoldTree()
    right_ft.add_edge(1, 179, -1)
    right_ft.add_edge(179, 180, 1)
    right_ft.add_edge(180, peptide_start_res - 1, -1)
    right_ft.add_edge(peptide_start_res - 1, peptide_start_res, 2)
    right_ft.add_edge(peptide_start_res, peptide_end_res, -1)

    init_pose.fold_tree(right_ft)

    for i in range(peptide_right_flank_start_res, peptide_end_res + 1):
        init_pose.set_phi(i, 180.0)
        init_pose.set_psi(i, 180.0)

    angle_max = 360
    n_iter = int(sys.argv[2] if len(sys.argv) > 2 else 1000)

    left_move_map = MoveMap()
    left_move_map.set_bb_true_range(peptide_start_res, peptide_left_flank_end_res)
    left_move_map.set_chi_true_range(peptide_start_res, peptide_left_flank_end_res)
    left_mover = SmallMover(left_move_map, 1.0, 10)
    left_mover.angle_max(angle_max)
    # left_minimizer = MinMover(left_move_map, get_fa_scorefxn(), "dfpmin", 0.01, True)

    right_move_map = MoveMap()
    right_move_map.set_bb_true_range(peptide_right_flank_start_res, peptide_end_res)
    right_move_map.set_chi_true_range(peptide_right_flank_start_res, peptide_end_res)
    right_mover = SmallMover(right_move_map, 1.0, 10)
    right_mover.angle_max(angle_max)
    # right_minimizer = MinMover(right_move_map, get_fa_scorefxn(), "dfpmin", 0.01, True)

    init_temp = 1.0

    dest_dir = os.path.join("./data/mhcii_ensembles/", pdb_id, f"annealing_max={angle_max}_iter={n_iter}")
    os.makedirs(dest_dir)
    ts = str(time.time())

    for i in tqdm(range(100)):

        q = deep_copy(init_pose)
        mc = MonteCarlo(q, get_fa_scorefxn(), init_temp)

        for j in tqdm(range(n_iter)):

            q.fold_tree(left_ft)
            left_mover.apply(q)

            q.fold_tree(right_ft)
            right_mover.apply(q)

            mc.boltzmann(q)

            mc.set_temperature(init_temp * 0.9)

        mc.recover_low(q)
        print(get_fa_scorefxn()(q))

        q.dump_pdb(os.path.join(dest_dir, f"{job_id}_{ts}_{i}.pdb"))
