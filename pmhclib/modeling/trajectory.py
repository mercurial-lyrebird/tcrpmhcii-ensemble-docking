
from __future__ import annotations

import os
import shutil

import numpy as np

import MDAnalysis as mda
from MDAnalysis import Universe
from MDAnalysis.analysis.rms import RMSD

from ..constants import DEFAULT_TRAJ_STORE
from ..components.structure import Structure
from ..components import Chain


class Trajectory:

    def __init__(
        self, topology: Structure, trajectory_fp: str, id: str = None
    ):

        self.topology = topology
        self.trajectory_fp = trajectory_fp
        self.universe = mda.Universe(self.topology.fp, self.trajectory_fp)

        [_id, self.ext] = os.path.basename(trajectory_fp).split(".")
        if id is not None:
            self.id = id
        else:
            self.id = _id
        self.store = os.path.join(DEFAULT_TRAJ_STORE, self.id)

        if os.path.exists(self.store):
            shutil.rmtree(self.store)
        os.mkdir(self.store)

        return

    def __len__(self) -> int:

        return len(self.universe.trajectory)

    def slice(self, t_0: int = None, t_f: int = None) -> Trajectory:

        slice_dir = os.path.join(self.store, "slices/")
        if not os.path.exists(slice_dir):
            os.mkdir(slice_dir)

        if t_0 is None:
            t_0 = 0
        if t_f is None:
            t_f = len(self.universe.trajectory)

        slice_id = self.id + f"_{t_0}-{t_f}"
        slice_fp = os.path.join(slice_dir, f"{slice_id}.{self.ext}")

        if not os.path.exists(slice_fp):
            self.universe.atoms.write(
                slice_fp, frames=self.universe.trajectory[t_0:t_f]
            )

        return Trajectory(self.topology, slice_fp)

    def select_chain(self, chain_id: str) -> Trajectory:

        selection_dir = os.path.join(self.store, "selections/")
        if not os.path.exists(selection_dir):
            os.mkdir(selection_dir)

        # get selection topology
        topology = Chain.from_complex(self.topology, chain_id=chain_id)

        # get selection trajectory
        traj_id = self.id + f"_chain{chain_id}"
        traj_fp = os.path.join(selection_dir, f"{traj_id}.xtc")
        if not os.path.exists(traj_fp):
            self.universe.select_atoms(f"chainID {chain_id}").write(
                traj_fp, frames="all"
            )

        return Trajectory(topology, traj_fp)

    def rmsd(self, ref: Structure, **kwargs) -> float:

        ref_universe = Universe(ref.fp)
        rmsd = RMSD(self.universe, ref_universe, **kwargs)
        rmsd.run()

        return rmsd.results.rmsd


class TrajectorySampler:

    def __init__(self, trajectory: Trajectory, pdb_dir: str):

        self.trajectory = trajectory
        self.pdb_dir = pdb_dir

        return

    def sample(self, mode: str, n: int = None, interval: int = None):

        if mode == "uniform":
            ts = np.sort(
                np.random.choice(np.arange(len(self.trajectory)), size=n)
            )
        elif mode == "interval":
            ts = np.arange(len(self.trajectory), step=interval)

        for i, _ in enumerate(self.trajectory.universe.trajectory[ts]):
            self.trajectory.universe.atoms.write(
                os.path.join(self.pdb_dir, f"t{ts[i]}.pdb")
            )

        return
