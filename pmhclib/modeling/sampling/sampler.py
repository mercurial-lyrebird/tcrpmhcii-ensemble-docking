
import abc
from multiprocessing import set_start_method, Pool
import os
from tqdm import tqdm
from typing import Dict, List, Tuple, Type

import pyrosetta
pyrosetta.init("-out:levels basic:error core:error")

from pyrosetta import Pose, pose_from_pdb, FoldTree, get_fa_scorefxn
from pyrosetta.rosetta.core.pose import deep_copy
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.protocols.simple_moves import SmallMover
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.protocols.moves import MonteCarlo
from pyrosetta.rosetta.protocols.loops import Loop, Loops
from pyrosetta.rosetta.protocols.loops.loop_mover.refine import (
    LoopMover_Refine_CCD
)

from ...components.structure import Chain, Structure
from ..ensemble import Ensemble


class Sampler(abc.ABC):

    def __init__(self, init_structure: Structure, chain_id: str = None):

        self.init_structure = init_structure
        self.chain_id = chain_id
        self.movable_ranges = None

        if isinstance(init_structure, Chain):
            self.chain_len = len(init_structure)
        else:
            if chain_id is None:
                raise ValueError(
                    "Must specify chain for sampling multi-chain structure"
                )
            self.chain_len = len(
                Chain.from_complex(init_structure, chain_id=self.chain_id)
            )

        return

    @abc.abstractmethod
    def sample(self) -> Pose:

        raise NotImplementedError


class LoopSampler(Sampler):

    def __init__(
        self, init_structure: Structure, chain_id: str = None,
        loops: List[Tuple[int, int]] = None, cut_ratios: List[float] = None,
        **kwargs
    ):

        super().__init__(init_structure, chain_id=chain_id, **kwargs)

        if loops is None:
            raise ValueError("Must specify loop residues")
        elif any([loop[1] - loop[0] < 0 for loop in loops]) < 0:
            raise ValueError("loop_start > loop_end")
        elif any(
            [loop[0] < 0 or loop[1] > self.chain_len - 1 for loop in loops]
        ):
            raise ValueError("invalid loop bounds")

        self.movable_ranges = loops

        if cut_ratios is None:
            cut_ratios = [0.5] * len(loops)

        self.cutpoints = [
            loops[i][0] + int(ratio * (loops[i][1] - loops[i][0]))
            for i, ratio in enumerate(cut_ratios)
        ]

        return


class TerminalRegionSampler(Sampler):

    def __init__(
        self, init_structure: Structure, chain_id: str = None,
        region_lens: Tuple[int, int] = None, **kwargs
    ):

        super().__init__(init_structure, chain_id=chain_id, **kwargs)

        if region_lens is None:
            raise ValueError("Must specify loop residues")

        self.movable_ranges = [
            [0, region_lens[0] - 1],
            [self.chain_len - region_lens[1], self.chain_len - 1]
        ]

        return


class RosettaSampler(Sampler):

    def __init__(
        self, init_structure: Structure, chain_id: str = None, **kwargs
    ):

        super().__init__(init_structure, chain_id=chain_id, **kwargs)

        self.init_pose = pose_from_pdb(self.init_structure.fp)
        self.res_chain_map = [
            self.init_pose.pdb_info().chain(i + 1)
            for i in range(self.init_pose.size())
        ]

        self.chain_id = chain_id
        self.chain_resis = [
            i + 1 for i in range(len(self.res_chain_map))
            if self.res_chain_map[i] == self.chain_id
        ]

        self.scoring_fxn = get_fa_scorefxn()

        return

    def map_movable_ranges(self) -> List[List[int]]:

        mapped_ranges = list()
        for movable_range in self.movable_ranges:
            mapped_ranges.append([self.chain_resis[i] for i in movable_range])

        return mapped_ranges

    def build_fold_tree(self, reverse_chain: bool = False) -> FoldTree:

        chain_start_resis = {
            chain_id: self.res_chain_map.index(chain_id)
            for chain_id in set(self.res_chain_map)
        }
        chain_ids_ordered = list(
            sorted(chain_start_resis.keys(), key=lambda x: chain_start_resis[x])
        )
        inactive_chain_ids_ordered = [
            chain_id for chain_id in chain_ids_ordered
            if chain_id != self.chain_id
        ]
        chain_ranges = {
            chain_ids_ordered[i]: [
                chain_start_resis[chain_ids_ordered[i]],
                (
                    len(self.res_chain_map) if i == len(chain_ids_ordered) - 1
                    else chain_start_resis[chain_ids_ordered[i + 1]]
                ) - 1
            ] for i in range(len(chain_ids_ordered))
        }

        ft = FoldTree()

        jump_i = 1
        for i in range(len(inactive_chain_ids_ordered)):
            this_chain = inactive_chain_ids_ordered[i]
            if i > 0:
                prev_chain = inactive_chain_ids_ordered[i - 1]
                ft.add_edge(
                    chain_ranges[prev_chain][1] + 1,
                    chain_ranges[this_chain][0] + 1,
                    jump_i
                )
                jump_i += 1
            ft.add_edge(
                chain_ranges[this_chain][0] + 1,
                chain_ranges[this_chain][1] + 1,
                -1
            )

        if reverse_chain:
            x = chain_ranges[self.chain_id][0]
            chain_ranges[self.chain_id][0] = chain_ranges[self.chain_id][1]
            chain_ranges[self.chain_id][1] = x

        ft.add_edge(
            chain_ranges[this_chain][1] + 1,
            chain_ranges[self.chain_id][0] + 1,
            jump_i
        )
        ft.add_edge(
            chain_ranges[self.chain_id][0] + 1,
            chain_ranges[self.chain_id][1] + 1,
            -1
        )

        return ft


class SequentialSampler(Sampler):

    @abc.abstractmethod
    def sample(self, n: int, ensemble_dir: str, ensemble_id: str) -> Ensemble:

        raise NotImplementedError


class RosettaLoopSampler(LoopSampler, RosettaSampler):

    def __init__(
        self, init_structure: Structure, chain_id: str = None,
        loops: List[Tuple[int, int]] = None, **kwargs
    ):

        super().__init__(
            init_structure, chain_id=chain_id, loops=loops, **kwargs
        )

        loop_ranges = self.map_movable_ranges()
        cutpoints = self.map_cutpoints()

        loops = Loops()
        for loop_range, cutpoint in zip(loop_ranges, cutpoints):
            loops.add_loop(Loop(loop_range[0], loop_range[1], cutpoint))

        self.loop_modeler = LoopMover_Refine_CCD()
        self.loop_modeler.loops(loops)
        self.loop_modeler.set_scorefxn(self.scoring_fxn)
        self.loop_modeler.set_fold_tree_from_loops(True)
        self.loop_modeler.max_inner_cycles(10)

        return

    def map_cutpoints(self) -> List[int]:

        return [self.chain_resis[cutpoint] for cutpoint in self.cutpoints]

    def sample(self) -> Pose:

        q = deep_copy(self.init_pose)
        self.loop_modeler.apply(q)

        return q


class BasinHoppingSampler(
    TerminalRegionSampler, RosettaSampler, SequentialSampler
):

    def sample(self, n: int, ensemble_dir: str, ensemble_id: str) -> Ensemble:

        ensemble_dir = os.path.join(ensemble_dir, ensemble_id)
        os.makedirs(ensemble_dir)

        [left_range, right_range] = self.map_movable_ranges()

        right_ft = self.build_fold_tree()
        left_ft = self.build_fold_tree(reverse_chain=True)

        left_move_map = MoveMap()
        left_move_map.set_bb_true_range(*left_range)
        left_move_map.set_chi_true_range(*left_range)
        left_mover = SmallMover(left_move_map, 1.0, 10)
        left_mover.angle_max(360.0)
        left_minimizer = MinMover(
            left_move_map, get_fa_scorefxn(), "dfpmin", 0.01, True
        )

        right_move_map = MoveMap()
        right_move_map.set_bb_true_range(*right_range)
        right_move_map.set_chi_true_range(*right_range)
        right_mover = SmallMover(right_move_map, 1.0, 10)
        right_mover.angle_max(360.0)
        right_minimizer = MinMover(
            right_move_map, get_fa_scorefxn(), "dfpmin", 0.01, True
        )

        q = deep_copy(self.init_pose)

        for i in tqdm(range(n)):

            q.fold_tree(left_ft)
            left_mover.apply(q)
            left_minimizer.apply(q)

            q.fold_tree(right_ft)
            right_mover.apply(q)
            right_minimizer.apply(q)

            q.dump_pdb(os.path.join(ensemble_dir, f"conf_{i}.pdb"))

        return Ensemble.from_pdb_dir(
            ensemble_id, ensemble_dir, fn_pattern="conf_*.pdb"
        )


class SimulatedAnnealingSampler(TerminalRegionSampler, RosettaSampler):

    def __init__(
        self, init_structure: Structure, chain_id: str = None,
        region_lens: Tuple[int, int] = None, n_iter: int = 1000,
        n_proc: int = None
    ):

        super().__init__(
            init_structure, region_lens=region_lens, chain_id=chain_id
        )

        self.n_iter = n_iter
        self.n_proc = n_proc

        [left_range, right_range] = self.map_movable_ranges()

        self.right_ft = self.build_fold_tree()
        self.left_ft = self.build_fold_tree(reverse_chain=True)

        left_move_map = MoveMap()
        left_move_map.set_bb_true_range(*left_range)
        left_move_map.set_chi_true_range(*left_range)
        self.left_mover = SmallMover(left_move_map, 1.0, 10)
        self.left_mover.angle_max(360.0)

        right_move_map = MoveMap()
        right_move_map.set_bb_true_range(*right_range)
        right_move_map.set_chi_true_range(*right_range)
        self.right_mover = SmallMover(right_move_map, 1.0, 10)
        self.right_mover.angle_max(360.0)

        return

    def sample(self) -> Pose:

        init_temp = 1.0

        q = deep_copy(self.init_pose)
        mc = MonteCarlo(q, self.scoring_fxn, init_temp)

        for _ in tqdm(range(self.n_iter)):

            q.fold_tree(self.left_ft)
            self.left_mover.apply(q)

            q.fold_tree(self.right_ft)
            self.right_mover.apply(q)

            mc.boltzmann(q)
            mc.set_temperature(init_temp * 0.9)

        mc.recover_low(q)

        return q


def sample(
    n: int, sampler_class: Type[Sampler], sampler_args: Dict,
    sampler_kwargs: Dict, ensemble_dir: str, ensemble_id: str,
    n_proc: int = None
) -> Ensemble:

    os.makedirs(ensemble_dir)

    args = [
        (
            sampler_class, sampler_args, sampler_kwargs,
            os.path.join(ensemble_dir, f"conf_{i}.pdb")
        )
        for i in range(n)
    ]

    set_start_method("spawn", force=True)
    with Pool(processes=n_proc) as pool:
        pool.map(_sample, args)

    return Ensemble.from_pdb_dir(
        ensemble_id, ensemble_dir, fn_pattern="conf_*.pdb"
    )


def _sample(args: Tuple[Type[Sampler], Dict, Dict, str]) -> Ensemble:

    import pyrosetta
    pyrosetta.init("-mute all")

    (sampler_class, sampler_args, sampler_kwargs, pdb_fp) = args

    sampler = sampler_class(*sampler_args, **sampler_kwargs)
    q = sampler.sample()
    q.dump_pdb(pdb_fp)

    return
