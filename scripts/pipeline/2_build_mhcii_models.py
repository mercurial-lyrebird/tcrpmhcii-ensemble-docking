
import sys

import pandas as pd

from pmhclib.components.pmhc import TCRpMHCII
from pmhclib.modeling.modeling import PANDORApMHCIIModeler


if __name__ == "__main__":

    modeler = PANDORApMHCIIModeler()
    library_df = pd.read_csv("../../data/mhcii_tcr_templates.csv", index_col="pdb_id")

    pdb_id = sys.argv[1]
    record = library_df.loc[pdb_id]
    complex = TCRpMHCII(
        f"../../data/mhcii_tcr_templates/{pdb_id}.pdb", id=pdb_id,
        anchors=[int(record[f"anch_{i + 1}"]) for i in range(4)]
    )
    modeler.model(
        pdb_id, complex.pmhc.mhc.alpha_chain.sequence,
        complex.pmhc.mhc.beta_chain.sequence, complex.pmhc.peptide.sequence,
        f"./{pdb_id}.pdb",
        anchors=complex.pmhc.anchors
    )
