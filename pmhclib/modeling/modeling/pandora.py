
import os
import shutil
from typing import List

from PANDORA import Target, Database
from PANDORA.Pandora.Pandora import Pandora

from pmhclib.modeling.modeling.modeler import (
    pMHCModeler, pMHCIModeler, pMHCIIModeler
)
from pmhclib.components.pmhc import (
    MHCISequence, MHCIISequence, MHCIPeptideSequence, MHCIIPeptideSequence,
    pMHCI, pMHCII
)
from pmhclib.components.pmhc.sequence import MHCIAllotype, MHCIIAllotype


PANDORA_WORKDIR = "/home/ab215/scratch/pandora/targets/"
os.makedirs(PANDORA_WORKDIR, exist_ok=True)


class PANDORApMHCModeler(pMHCModeler):

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

        self.database = Database.load()

        return

    def _model(self, target: Target) -> str:

        model = Pandora(target, database=self.database)
        model.find_template(best_n_templates=1, benchmark=True, verbose=False)
        model.copy_template()
        model.clip_C_domain = False
        model.align(verbose=False)
        model.write_ini_script()
        model.create_initial_model(verbose=False)

        return os.path.join(PANDORA_WORKDIR, target.id, f"{target.id}.ini")


class PANDORApMHCIModeler(pMHCIModeler, PANDORApMHCModeler):

    def model(
        self, target_id: str, alpha_seq: MHCISequence,
        peptide_seq: MHCIPeptideSequence, output_fp: str,
        allotype: MHCIAllotype = None, anchors: List[int] = None
    ) -> pMHCI:

        if allotype is None:
            allotype = list()
        if anchors is None:
            anchors = list()

        target = Target(
            target_id, peptide=str(peptide_seq), MHC_class="I",
            M_chain_seq=str(alpha_seq), allele_type=allotype, anchors=anchors,
            output_dir=PANDORA_WORKDIR
        )
        model_fp = self._model(target)
        shutil.copy(model_fp, output_fp)

        return pMHCI(
            output_fp, mhc_chain_id="M", peptide_chain_id="P",
            allotype=allotype, anchors=anchors
        )


class PANDORApMHCIIModeler(pMHCIIModeler, PANDORApMHCModeler):

    def model(
        self, target_id: str, alpha_seq: MHCIISequence, beta_seq: MHCIISequence,
        peptide_seq: MHCIIPeptideSequence, output_fp: str,
        allotype: MHCIIAllotype = None, anchors: List[int] = None
    ) -> pMHCII:

        if allotype is None:
            allotype = list()
        if anchors is None:
            anchors = list()

        target = Target(
            target_id, peptide=str(peptide_seq), MHC_class="II",
            M_chain_seq=str(alpha_seq), N_chain_seq=str(beta_seq),
            allele_type=allotype, anchors=anchors,
            output_dir=PANDORA_WORKDIR
        )
        model_fp = self._model(target)
        shutil.copy(model_fp, output_fp)

        return pMHCII(output_fp, allotype=allotype, anchors=anchors)