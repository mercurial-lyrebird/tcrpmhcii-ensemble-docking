
from typing import List

from PANDORA import Target
from PANDORA import Pandora
from PANDORA import Database

from ..components.pmhc.structure.pmhc import pMHC, pMHCI, pMHCII


def run_pandora(
    case_id: str, pmhc: pMHC, alleles: List[str], anchors: List[int], **kwargs
):

    if isinstance(pmhc, pMHCI):
        mhc_class = "I"
    elif isinstance(pmhc, pMHCII):
        mhc_class = "II"

    db = Database.load()
    target = Target(
        id=case_id, MHC_class=mhc_class, peptide=str(pmhc.peptide.sequence),
        allele_type=alleles, M_chain_seq=str(pmhc.mhc.alpha_chain.sequence),
        N_chain_seq=str(pmhc.mhc.beta_chain.sequence), anchors=anchors
    )

    run = Pandora.Pandora(target, db)
    run.model(**kwargs)

    return
