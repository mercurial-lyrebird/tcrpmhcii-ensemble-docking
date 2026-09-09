
import sys

import pandas as pd

from tfold.model import esm_ppi_650m_tcr, tfold_tcr_trunk
from tfold.deploy import TCRPredictor


if __name__ == "__main__":

    pdb_id = sys.argv[1]

    ppi_model_path = esm_ppi_650m_tcr()
    tfold_model_path = tfold_tcr_trunk()

    model = TCRPredictor.restore_from_module(ppi_model_path, tfold_model_path)

    seq_df = pd.read_csv("./data/mhcii_tcr_sequences.csv", index_col="pdb_id")
    model.infer_pdb(
        [
            {"id": "A", "sequence": seq_df.loc[pdb_id].A},
            {"id": "B", "sequence": seq_df.loc[pdb_id].B}
        ],
        f"./data/tcr_models/tfold-tcr/{pdb_id}.pdb"
    )
