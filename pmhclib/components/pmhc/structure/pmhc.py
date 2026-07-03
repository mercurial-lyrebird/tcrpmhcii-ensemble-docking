
import abc
import os
from typing import List, Tuple, Union

import numpy as np

from ...structure import Complex
from .peptide import MHCPeptide, MHCIPeptide, MHCIIPeptide
from .mhc import MHC, MHCI, MHCII
from ..sequence import MHCIAllotype, MHCIIAllotype


class pMHC(abc.ABC, Complex):

    def __init__(self, fp: str, id: str = None):

        super().__init__(fp, id=id)

        self.mhc: MHC = None
        self.peptide: MHCPeptide = None

        return

    @property
    def allotype(self) -> Union[MHCIAllotype, MHCIIAllotype]:

        return self.mhc.allotype

    def infer_allotype(self):

        self.mhc.infer_allotype()

        return

    def set_allotype(self, allotype: Union[MHCIAllotype, MHCIIAllotype]):

        self.mhc.set_allotype(allotype)

        return

    @property
    def anchors(self) -> List[int]:

        return self.peptide.anchors

    def infer_anchors(self):

        self.peptide.infer_anchors(self.mhc.allotype)

        return

    def set_anchors(self, anchors: List[int]):

        self.peptide.set_anchors(anchors)

        return

    def peptide_sasa(self, core_only: bool = False) -> np.ndarray:

        if core_only:
            res_i = range(self.anchors[0], self.anchors[-1] + 1)
        else:
            res_i = None

        return super().sasa(chain_id=self.peptide_chain_id, res_i=res_i)

    def estimate_anchors(self):

        peptide_sasa = self.peptide_sasa()
        spacing_sasas = [sum([peptide_sasa[i + j] for j in [0, 3, 5, 8]]) for i in range(len(self.peptide) - 9)]

        self.set_anchors([np.array(spacing_sasas).argmin() + 1 + j for j in [0, 3, 5, 8]])

        return

    def delete(self, rec: bool = False):

        os.remove(self.fp)

        if rec:
            self.mhc.delete(rec=True)
            self.peptide.delete()

        return


class pMHCI(pMHC):

    def __init__(
        self, fp: str, mhc_chain_id: str = "A", peptide_chain_id: str = "C",
        id: str = None, allotype: MHCIAllotype = None, anchors: List[int] = None
    ):

        super().__init__(fp, id=id)

        self.mhc = MHCI.from_complex(
            self, chain_ids=[mhc_chain_id], chain_id=mhc_chain_id,
            allotype=allotype
        )

        self.peptide_chain_id = peptide_chain_id
        self.peptide = MHCIPeptide.from_complex(
            self, chain_id=peptide_chain_id, anchors=anchors
        )

        return


class pMHCII(pMHC):

    def __init__(
        self, fp: str, mhc_chain_ids: Tuple[str, str] = ("M", "N"),
        peptide_chain_id: str = "P", id: str = None,
        allotype: MHCIIAllotype = None, anchors: List[int] = None
    ):

        super().__init__(fp, id=id)

        self.mhc = MHCII.from_complex(
            self, chain_ids=mhc_chain_ids, allotype=allotype
        )

        self.peptide_chain_id = peptide_chain_id
        self.peptide = MHCIIPeptide.from_complex(
            self, chain_id=peptide_chain_id, anchors=anchors
        )

        return
