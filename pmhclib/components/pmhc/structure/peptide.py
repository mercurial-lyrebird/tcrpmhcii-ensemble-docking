
from __future__ import annotations

import abc
from typing import List, Union

import numpy as np

from ...structure import Chain
from ..sequence import (
    MHCIAllotype, MHCIIAllotype, MHCPeptideSequence, MHCIPeptideSequence,
    MHCIIPeptideSequence
)


class MHCPeptide(abc.ABC, Chain):

    def __init__(self, fp: str, id: str = None):

        super().__init__(fp, id=id)

        self.sequence: MHCPeptideSequence = None

        return

    @property
    def anchors(self) -> List[int]:

        return self.sequence.anchors

    def infer_anchors(self, allotype: Union[MHCIAllotype, MHCIIAllotype]):

        self.sequence.infer_anchors(allotype)

        return

    def set_anchors(self, anchors: List[int]):

        self.sequence.set_anchors(anchors)

        return

    @property
    def binding_core(self) -> List[int]:

        return self.sequence.binding_core


class MHCIPeptide(MHCPeptide):

    def __init__(
        self, fp: str, id: str = None, anchors: List[int] = None
    ):

        super().__init__(fp, id=id)

        self.sequence = MHCIPeptideSequence(
            seq=self._load_sequence_str(), anchors=anchors
        )

        return


class MHCIIPeptide(MHCPeptide):

    def __init__(
        self, fp: str, id: str = None, anchors: List[int] = None
    ):

        super().__init__(fp, id=id)

        self.sequence = MHCIIPeptideSequence(
            seq=self._load_sequence_str(), anchors=anchors
        )

        return

    @property
    def flanking_regions(self) -> List[List[int]]:

        return self.sequence.flanking_regions
