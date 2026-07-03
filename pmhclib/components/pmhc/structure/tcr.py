
from __future__ import annotations

import os
from typing import List, Tuple

import numpy as np

from ...structure import Chain, Multimer
from ..sequence import TCRSequence


class TCR(Multimer):

    def __init__(
        self, fp: str, id: str = None, chain_ids: Tuple[str] = ("D", "E")
    ):

        super().__init__(fp, id=id)

        self._numbering = None
        self.chain_ids = chain_ids

        self.alpha_chain = TCRChain.from_multimer(self, chain_id=self.chain_ids[0])
        self.beta_chain = TCRChain.from_multimer(self, chain_id=self.chain_ids[1])

        return

    def get_ab_vector(self) -> np.ndarray:

        vector = self.beta_chain.center_of_mass() - self.alpha_chain.center_of_mass()

        return vector / np.linalg.norm(vector)

    def delete(self, rec: bool = False):

        os.remove(self.fp)

        if rec:
            self.alpha_chain.delete()
            self.beta_chain.delete()

        return


class TCRChain(Chain):

    def __init__(self, fp: str, id: str = None):

        super().__init__(fp, id=id)

        self.sequence = TCRSequence(seq=self._load_sequence_str())

        return

    @property
    def numbering(self) -> List[int]:

        return self.sequence.numbering

    @property
    def cdr_loops(self) -> List[List[int]]:

        return self.sequence.cdr_loops
