
from __future__ import annotations

from glob import glob
import os
from tqdm import tqdm
from typing import List

import numpy as np
import torch
from torch import Tensor

from pymol import cmd

from ...representations._representation import Representation
from ...components.structure import Structure
from ..scoring.scorer import Scorer, RosettaScorer


class Ensemble:

    def __init__(self, id: str, pdb_fps: List[str]):

        self.id = id

        self.pdb_fps = pdb_fps
        self.ids = [
            self.id + "_" + ".".join(os.path.basename(pdb_fp).split(".")[:-1])
            for pdb_fp in self.pdb_fps
        ]
        self.id_to_fp = {
            conf_id: pdb_fp for conf_id, pdb_fp in zip(self.ids, self.pdb_fps)
        }
        self.n_conf = len(self.id_to_fp)

        return

    @classmethod
    def from_pdb_dir(
        cls, id: str, pdb_dir: str, fn_pattern: str = None
    ) -> Ensemble:

        if fn_pattern is None:
            fn_pattern = "*.pdb"

        return Ensemble(id, sorted(glob(os.path.join(pdb_dir, fn_pattern))))

    def apply_representation(
        self, model: Representation, selection: str = "all"
    ) -> Tensor:

        X = list()
        for conf_id in tqdm(self.ids):
            structure = Structure(self.id_to_fp[conf_id], id=conf_id)
            structure = structure.select(selection)
            X.append(model(structure))
            structure.delete()

        return torch.stack(X)

    def rmsd(self, ref_pdb_fp: str, selection: str = "all") -> Tensor:

        rmsds = list()
        for conf_id in tqdm(self.ids):
            structure = Structure(fp=self.id_to_fp[conf_id], id=conf_id)
            rmsds.append(structure.rmsd(ref_pdb_fp, selection=selection))
            structure.delete()

        return Tensor(rmsds)

    def score(self, scorer: Scorer = None) -> List:

        if scorer is None:
            scorer = RosettaScorer()

        structures = [
            Structure(fp=self.id_to_fp[conf_id], id=conf_id)
            for conf_id in self.ids
        ]
        scores = [scorer(structure) for structure in tqdm(structures)]

        for structure in structures:
            structure.delete()

        return scores

    def filter(
        self, id: str, scores: List[float] = None, scorer: Scorer = None,
        cutoff: float = None, percentile: float = None
    ) -> Ensemble:

        if scores is not None:
            pass
        elif scorer is not None:
            scores = self.score(scorer)
        else:
            raise ValueError

        scores = np.array(scores)

        if cutoff is not None:
            idx = np.where(scores <= cutoff)[0]
        elif percentile is not None:
            raise NotImplementedError
        else:
            raise ValueError

        return Ensemble(id, [self.pdb_fps[i] for i in idx])

    def to_pymol(self, pymol_fp: str):

        for conf_id, pdb_fp in self.id_to_fp.items():
            cmd.load(pdb_fp, conf_id)

        cmd.save(pymol_fp)

        for conf_id, pdb_fp in self.id_to_fp.items():
            cmd.delete(conf_id)

        return
