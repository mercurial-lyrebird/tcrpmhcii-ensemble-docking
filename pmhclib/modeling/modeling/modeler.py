
from abc import ABC, abstractmethod

from pmhclib.components import Sequence
from pmhclib.components.structure import Structure
from pmhclib.components.pmhc import (
    MHCIPeptideSequence, MHCIIPeptideSequence, MHCISequence, MHCIISequence,
    pMHCI, pMHCII, TCR
)
from pmhclib.components.pmhc.sequence import MHCSequence, MHCPeptideSequence
from pmhclib.components.pmhc.structure.pmhc import pMHC


class Modeler(ABC):

    def __init__(self):

        return

    @abstractmethod
    def model(self, *seqs: str) -> Structure:

        raise NotImplementedError


class TCRModeler(Modeler):

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

        return

    @abstractmethod
    def model(self, alpha_seq: Sequence, beta_seq: Sequence) -> TCR:

        raise NotImplementedError


class pMHCModeler(Modeler):

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

        return

    @abstractmethod
    def model(
        self, *seqs: MHCSequence, peptide_seq: MHCPeptideSequence
    ) -> pMHC:

        raise NotImplementedError


class pMHCIModeler(pMHCModeler):

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

        return

    @abstractmethod
    def model(
        self, alpha_seq: MHCISequence, peptide_seq: MHCIPeptideSequence
    ) -> pMHCI:

        raise NotImplementedError


class pMHCIIModeler(pMHCModeler):

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

        return

    @abstractmethod
    def model(
        self, alpha_seq: MHCIISequence, beta_seq: MHCIISequence,
        peptide_seq: MHCIIPeptideSequence
    ) -> pMHCII:

        raise NotImplementedError


class tFoldModeler(TCRModeler):

    def model(self, alpha_seq: Sequence, beta_seq: Sequence) -> TCR:

        raise NotImplementedError
