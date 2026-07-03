
from __future__ import annotations

from abc import ABC, abstractmethod
import os
from typing import Tuple, Union

from ...structure import Chain, Multimer
from ..sequence import (
    MHCAllotype, MHCIAllotype, MHCIIAllotype, MHCSequence, MHCISequence,
    MHCIISequence
)


class MHC(ABC, Multimer):

    def __init__(self, fp: str, id: str = None):

        super().__init__(fp, id=id)

        return

    @property
    @abstractmethod
    def allotype(self) -> Union[MHCIAllotype, MHCIIAllotype]:

        raise NotImplementedError

    @abstractmethod
    def infer_allotype(self):

        raise NotImplementedError

    @abstractmethod
    def set_allotype(self, allotype: Union[MHCIAllotype, MHCIIAllotype]):

        raise NotImplementedError


class MHCChain(ABC, Chain):

    def __init__(self, fp: str, id: str = None, chain_id: str = "M"):

        super().__init__(fp, id=id)

        self.sequence: MHCSequence = None
        self.chain_id = chain_id

        return

    @property
    def allotype(self) -> MHCAllotype:

        return self.sequence.allotype

    def infer_allotype(self):

        self.sequence.infer_allotype(self.chain_id)

        return

    def set_allotype(self, allotype: Union[MHCIAllotype, MHCIIAllotype]):

        self.sequence.set_allotype(allotype)

        return


class MHCIChain(MHCChain):

    def __init__(
        self, fp: str, id: str = None, chain_id: str = "A",
        allotype: MHCIAllotype = None
    ):

        super().__init__(fp, id=id, chain_id=chain_id)

        self.sequence = MHCISequence(
            seq=self._load_sequence_str(), allotype=allotype
        )
        self.chain_id = chain_id

        return


class MHCIIChain(MHCChain):

    def __init__(self, fp: str, id: str = None, chain_id: str = "M"):

        super().__init__(fp, id=id, chain_id=chain_id)

        self.sequence = MHCIISequence(seq=self._load_sequence_str())
        self.chain_id = chain_id

        return

    def infer_allotype(self):

        raise NotImplementedError

    def set_allotype(self, allotype: MHCIIAllotype):

        raise NotImplementedError


class MHCI(MHC):

    def __init__(self, fp: str, id: str = None, chain_id: str = "A", allotype: MHCIAllotype = None):

        super().__init__(fp, id=id)

        self.chain_id = chain_id
        self.alpha_chain: MHCIChain = MHCIChain.from_complex(
            self, chain_id=self.chain_id, allotype=allotype
        )

        return

    @property
    def allotype(self) -> MHCIAllotype:

        return self.alpha_chain.allotype

    def infer_allotype(self):

        self.alpha_chain.infer_allotype()

        return

    def set_allotype(self, allotype: Union[MHCIAllotype, MHCIIAllotype]):

        self.alpha_chain.set_allotype(allotype)

        return

    def delete(self, rec: bool = False):

        os.remove(self.fp)

        if rec:
            self.alpha_chain.delete()

        return


class MHCII(MHC):

    def __init__(
        self, fp: str, id: str = None, chain_ids: Tuple[str] = ("M", "N"),
        allotype: MHCIIAllotype = None
    ):

        super().__init__(fp, id=id)

        self.chain_ids = chain_ids
        self.alpha_chain: MHCIIChain = MHCIIChain.from_multimer(
            self, chain_id=self.chain_ids[0]
        )
        self.beta_chain: MHCIIChain = MHCIIChain.from_multimer(
            self, chain_id=self.chain_ids[1]
        )

        self._allotype = allotype

        return

    @property
    def allotype(self) -> MHCIIAllotype:

        return self._allotype

    def infer_allotype(self):

        raise NotImplementedError

    def set_allotype(self, allotype: MHCIIAllotype):

        self._allotype = allotype

        return

    def delete(self, rec: bool = False):

        os.remove(self.fp)

        if rec:
            self.alpha_chain.delete()
            self.beta_chain.delete()

        return
