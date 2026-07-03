
from math import pi
from typing import List

import numpy as np
import torch
from torch import Tensor, cos

import pyrosetta
pyrosetta.init("-out:levels basic:error core:error")

from pymol import cmd
from pyrosetta import pose_from_pdb

from ._representation import Representation
from ..components import Chain


# TODO: implement other proxies (e.g. geometric center, center of mass)

class CoordinateVector(Representation):

    def __init__(self, proxy: str = "c_alpha", flatten: bool = True):

        super(CoordinateVector, self).__init__()

        self.proxy = proxy
        self.flatten = flatten

        return

    def forward(self, chain: Chain) -> Tensor:

        cmd.load(chain.fp, chain.id)

        if self.proxy == "c_alpha":
            X = cmd.get_coords(f"{chain.id} & name CA")
        else:
            raise ValueError

        cmd.delete(chain.id)

        return Tensor(X)


class TorsionVector(Representation):

    def __init__(
        self, angles: List[str] = ["psi", "phi"], cosine: bool = False,
        flatten: bool = True
    ):

        super(TorsionVector, self).__init__()

        self.angles = angles
        self.cosine = cosine
        self.flatten = flatten

        return

    # TODO: probably a sleeker way to support chi
    # TODO: likely buggy for chi_n > 1

    def forward(self, chain: Chain) -> Tensor:

        pose = pose_from_pdb(chain.fp)
        torsion_vector = list()
        for i in range(len(pose)):
            res_torsion_vector = list()
            for angle in self.angles:
                if angle.startswith("chi"):
                    if pose.aa(i + 1) in (1, 6):  # is alanine or glycine
                        val = 0.0
                    else:
                        suffix = angle[3:]
                        if len(angle) < 4 or not suffix.isnumeric():
                            raise ValueError(angle)
                        chi_n = int(suffix)
                        val = pose.chi(chi_n, i + 1)
                else:
                    val = getattr(pose, angle)(i + 1)
                res_torsion_vector.append(val)
            torsion_vector.append(res_torsion_vector)

        X = Tensor(torsion_vector)
        X += (360 * (X < 0))
        if self.cosine:
            X *= (pi / 180)
            X = cos(X)

        if self.flatten:
            X = X.flatten()

        return X


class BFactorVector(Representation):

    def __init__(self):

        super(BFactorVector, self).__init__()

        return

    def forward(self, chain: Chain) -> Tensor:

        pose = pose_from_pdb(chain.fp)
        pdb_info = pose.pdb_info()

        return Tensor(
            [
                np.mean(
                    [
                        pdb_info.bfactor(i + 1, j + 1)
                        for j in range(len(pose[i + 1].atoms()))
                    ]
                )
                for i in range(len(pose))
            ]
        )
