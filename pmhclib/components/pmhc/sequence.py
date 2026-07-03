
from abc import ABC, abstractmethod
import os
import random as rd
from typing import List, Tuple, Union

import pandas as pd

from anarci import anarci

from PANDORA import PANDORA_data
from PANDORA.Pandora.Modelling_functions import (
    blast_mhc_seq,
    predict_anchors_netMHCpan,
    predict_anchors_netMHCIIpan
)

from .. import Sequence


MHCAllotype = str
MHCIAllotype = str
MHCIIAllotype = Tuple[str, str]


REFSEQ_BLAST_PATH = "BLAST_databases/refseq_blast_db/refseq_blast_db"
IMGT_LOOPS = [[27, 38], [56, 66], [105, 117]]


class MHCSequence(ABC, Sequence):

    def __init__(self, seq: str, allotype: MHCAllotype = None):

        super().__init__(seq)

        self.allotype = allotype

        return

    def infer_allotype(self, chain_id: str):

        hits = blast_mhc_seq(
            str(self), chain=chain_id,
            blastdb=os.path.join(PANDORA_data, REFSEQ_BLAST_PATH)
        )

        self.allotype = hits[0][0]

        return

    def set_allotype(self, allotype: MHCAllotype):

        self.allotype = allotype

        return

    # @abstractmethod
    # def get_domain(self, domain: str) -> List[int]:

    #     # TODO

    #     raise NotImplementedError


class MHCISequence(MHCSequence):

    pass


class MHCIISequence(MHCSequence):

    def __init__(self, seq: str):

        super().__init__(seq)

        return

    def infer_allotype(self):

        raise NotImplementedError

    def set_allotype(self, allotype: MHCIIAllotype):

        raise NotImplementedError


class MHCPeptideSequence(ABC, Sequence):

    def __init__(
        self, seq: str, anchors: List[int] = None
    ):

        super().__init__(seq)

        self.anchors = anchors

        return

    @abstractmethod
    def infer_anchors(self, allotype: Union[MHCIAllotype, MHCIIAllotype]):

        raise NotImplementedError

    def set_anchors(self, anchors: List[int]):

        self.anchors = anchors

        return


class MHCIPeptideSequence(MHCPeptideSequence):

    def __init__(self, seq: str, anchors: List[int] = None):

        super().__init__(seq, anchors=anchors)

        return

    def infer_anchors(self, allotype: MHCIAllotype):

        self.anchors = predict_anchors_netMHCpan(
            self.sequence, [] if allotype is None else [allotype], "."
        )

        return

    @property
    def binding_core(self) -> List[int]:

        return [self.anchors[0], self.anchors[1]]


class MHCIIPeptideSequence(MHCPeptideSequence):

    def __init__(self, seq: str, anchors: List[int] = None):

        super().__init__(seq, anchors=anchors)

        return

    def infer_anchors(self, allotype: MHCIAllotype):

        self.anchors = predict_anchors_netMHCIIpan(
            self.sequence, [*allotype], "."
        )

        return

    @property
    def binding_core(self) -> List[int]:

        return [self.anchors[0], self.anchors[-1]]

    @property
    def flanking_regions(self) -> List[List[int]]:

        return [[1, self.anchors[0] - 1], [self.anchors[-1] + 1, len(self)]]


class TCRSequence(Sequence):

    def __init__(self, seq: str):

        super().__init__(seq)

        self._numbering = None

        return

    @property
    def numbering(self) -> List[int]:

        if self._numbering is None:
            self._numbering = get_anarci_numbering(self.sequence)

        return self._numbering

    @property
    def cdr_loops(self) -> List[List[int]]:

        return [
            [
                self.numbering.index(loop[0]) + 1,
                self.numbering.index(loop[1]) + 1
            ]
            for loop in IMGT_LOOPS
        ]


def get_anarci_numbering(tcr_seq: str):

    id = "".join([str(rd.randint(0, 9)) for _ in range(10)])
    anarci_fp = f"anarci_{id}.txt"
    anarci_stripped_fp = f"anarci_stripped_{id}.txt"
    anarci_stripped_clean_fp = f"anarci_stripped_clean_{id}.txt"

    TCR_seq = [("seq", tcr_seq)]
    anarci(TCR_seq, scheme="imgt", output=True, outfile=anarci_fp)
    os.system(f"tr -s '[:blank:]' ',' < {anarci_fp} > {anarci_stripped_fp}")

    # remove cases with additional letters on rows...
    with (
        open(anarci_stripped_fp, "r") as infile,
        open(anarci_stripped_clean_fp, "w") as outfile
    ):
        for line in infile:
            if line.count(",") > 2:
                line = line[:-4]+line[-2:]
            outfile.write(line)

    renumbered_TCR = pd.read_csv(
        anarci_stripped_clean_fp, names=["chain", "idx", "acid"], comment="#"
    )
    renumbered_TCR = renumbered_TCR[renumbered_TCR["acid"] != "-"].dropna()

    numbering = list()

    prefix = "".join(renumbered_TCR.acid.values[:5])
    offset = tcr_seq.index(prefix)

    for _ in range(0, offset):
        numbering.append(-1)

    numbering.extend([int(row.idx) for _, row in renumbered_TCR.iterrows()])

    for _ in range(len(numbering), len(tcr_seq)):
        numbering.append(-1)

    os.remove(anarci_fp)
    os.remove(anarci_stripped_fp)
    os.remove(anarci_stripped_clean_fp)

    return numbering
