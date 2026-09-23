import os
import requests
from tqdm import tqdm
from typing import List

import pandas as pd

from pymol import cmd

from pmhclib.components.pmhc import TCRpMHCII


canonical_chain_ids = ["M", "N", "P", "A", "B"]
temporary_chain_ids = ["V", "W", "X", "Y", "Z"]


ORIG_CHAIN_IDS = {
    "1D9K": ["C", "D", "P", "A", "B"],
    "1FYT": ["A", "B", "C", "D", "E"],
    "1J8H": ["A", "B", "C", "D", "E"],
    "1U3H": ["C", "D", "P", "A", "B"],
    "1YMM": ["A", "B", "C", "D", "E"],
    "1ZGL": ["A", "B", "C", "M", "P"],
    "2IAM": ["A", "B", "P", "C", "D"],
    "2IAN": ["A", "B", "C", "D", "E"],
    "2PXY": ["C", "D", "P", "A", "B"],
    "2Z31": ["C", "D", "P", "A", "B"],
    "3MBE": ["A", "B", "P", "C", "D"],
    "3QIB": ["A", "B", "P", "C", "D"],
    "3QIU": ["A", "B", "E", "C", "D"],
    "3QIW": ["A", "B", "E", "C", "D"],
    "4E41": ["A", "B", "C", "D", "E"],
    "4GG6": ["A", "B", "J", "G", "H"],
    "4H1L": ["A", "B", "C", "I", "J"],
    "4OZF": ["A", "B", "J", "G", "H"],
    "4OZG": ["A", "B", "J", "G", "H"],
    "4OZH": ["A", "B", "J", "G", "H"],
    "4OZI": ["A", "B", "J", "G", "H"],
    "4P2O": ["A", "B", "P", "C", "D"],
    "4P2Q": ["A", "B", "C", "D", "E"],
    "4P2R": ["A", "B", "C", "D", "E"],
    "4Y19": ["A", "B", "C", "D", "E"],
    "4Y1A": ["A", "B", "C", "D", "E"],
    "4Z7U": ["A", "B", "J", "G", "H"],
    "4Z7V": ["A", "B", "J", "G", "H"],
    "4Z7W": ["A", "B", "J", "G", "H"],
    "5KS9": ["A", "B", "J", "G", "H"],
    "5KSA": ["A", "B", "J", "C", "D"],
    "5KSB": ["A", "B", "J", "G", "H"],
    "6BGA": ["A", "B", "P", "C", "D"],
    "6CQL": ["A", "B", "C", "D", "E"],
    "6CQN": ["A", "B", "C", "D", "E"],
    "6CQQ": ["A", "B", "C", "D", "E"],
    "6CQR": ["A", "B", "C", "D", "E"],
    "6PX6": ["A", "B", "C", "D", "E"],
    "6PY2": ["A", "B", "C", "D", "E"],
    "6R0E": ["A", "B", "C", "D", "E"],
    "6V13": [], # ["A", "B", "C", "D", "E"],            # PTMs in peptide
    "6V15": [], # ["A", "B", "C", "D", "E"],            # PTMs in peptide
    "6V18": [], # ["A", "B", "C", "D", "E"],            # PTMs in peptide
    "6V19": [], # ["A", "B", "C", "D", "E"],            # PTMs in peptide
    "6V1A": [], # ["A", "B", "C", "D", "E"],            # PTMs in peptide
    "7SG1": ["A", "B", "C", "D", "E"],
    "7SG2": ["A", "B", "C", "D", "E"],
    "7T2B": ["A", "B", "C", "D", "E"],
    "7T2C": ["A", "B", "C", "D", "E"],
    "7T2D": ["A", "B", "C", "D", "E"],
    "7Z50": ["A", "B", "T", "H", "E"],
    "8PJG": ["A", "B", "C", "D", "E"],
    "8TRL": [], # ["A", "B", "C", "I", "J"],            # PTMs in peptide
    "8TRR": [], # ["A", "B", "C", "D", "E"],            # PTMs in peptide
    "8VCX": ["A", "B", "C", "D", "E"],
    "8VCY": ["A", "B", "C", "D", "E"],
    "8VD2": ["A", "B", "C", "D", "E"]
}


def handle_6BGA(raw_pdb_fp: str):

    cmd.load(raw_pdb_fp, "6BGA")
    cmd.alter("chain B & resi \-24-\-5", "chain='P'")
    cmd.save(raw_pdb_fp, "6BGA")
    cmd.delete("6BGA")

    return


if __name__ == "__main__":

    stcrdab_df = pd.read_csv("../../data/stcrdab_mhcii_tcr.tsv", sep="\t")
    pandora_df = pd.read_csv("../../data/mhcii_templates.csv", index_col="pdb_id")
    pdb_ids = stcrdab_df.pdb.unique()

    # os.makedirs("../../data/mhcii_tcr_templates/")
    # os.makedirs("../../data/mhcii_tcr_templates_raw/")

    for _pdb_id in tqdm(pdb_ids):
        pdb_id = _pdb_id.upper()
        # if pdb_id not in pandora_df.index:
        #     continue
        if len(ORIG_CHAIN_IDS[pdb_id]) == 0:
            continue
        url = f"https://files.rcsb.org/download/{_pdb_id}.pdb"
        response = requests.get(url)
        if response.status_code == 404:
            continue
        raw_pdb_fp = f"../../data/mhcii_tcr_templates_raw/{pdb_id}.pdb"
        with open(raw_pdb_fp, "wb") as f:
            f.write(response.content)
        if pdb_id == "6BGA":
            handle_6BGA(raw_pdb_fp)
        cmd.load(raw_pdb_fp, "complex")
        cmd.create(
            pdb_id, "complex & chain " + "+".join(ORIG_CHAIN_IDS[pdb_id])
        )
        cmd.delete("complex")
        for i in range(5):
            cmd.alter(
                f"{pdb_id} & chain {ORIG_CHAIN_IDS[pdb_id][i]}",
                f"chain='{temporary_chain_ids[i]}'"
            )
        for i in range(5):
            cmd.alter(
                f"{pdb_id} & chain {temporary_chain_ids[i]}",
                f"chain='{canonical_chain_ids[i]}'"
            )
        temp_pdb_fp = f"../../data/mhcii_tcr_templates/_{pdb_id}.pdb"
        clean_pdb_fp = f"../../data/mhcii_tcr_templates/{pdb_id}.pdb"
        cmd.save(temp_pdb_fp, pdb_id)
        cmd.delete(pdb_id)
        os.system(
            f"pdb_tidy {temp_pdb_fp} | pdb_delhetatm | grep -v '^ANISOU' > {clean_pdb_fp}"
        )
        os.remove(temp_pdb_fp)
        os.system(f"pdb_splitchain {clean_pdb_fp}")
        for chain_id in canonical_chain_ids:
            os.system(f"pdb_reres -1 {pdb_id}_{chain_id}.pdb > {pdb_id}_{chain_id}_reres.pdb")
        os.system("pdb_merge " + " ".join([f"{pdb_id}_{chain_id}_reres.pdb" for chain_id in canonical_chain_ids]) + f" > {temp_pdb_fp}")
        for chain_id in canonical_chain_ids:
            os.remove(f"{pdb_id}_{chain_id}.pdb")
            os.remove(f"{pdb_id}_{chain_id}_reres.pdb")
        os.remove(clean_pdb_fp)
        os.system(f"pdb_reatom {temp_pdb_fp} | pdb_tidy > {clean_pdb_fp}")
        os.remove(temp_pdb_fp)

    pdb_ids = sorted(
        [
            pdb_id.upper() for pdb_id in pdb_ids
            if os.path.exists(f"../../data/mhcii_tcr_templates/{pdb_id.upper()}.pdb")
        ]
    )

    seq_configs = list()

    for pdb_id in pdb_ids:
        seq_config = dict(pdb_id=pdb_id)
        cmd.load(f"../../data/mhcii_tcr_templates/{pdb_id}.pdb", pdb_id)
        for chain_id in canonical_chain_ids:
            seq_config[chain_id] = cmd.get_fastastr(
                f"{pdb_id} & chain {chain_id}"
            )[(1 + len(pdb_id) + 2):].replace("\n", "")
        cmd.delete(pdb_id)
        seq_configs.append(seq_config)

    seq_df = pd.DataFrame.from_dict(seq_configs).set_index("pdb_id")
    seq_df.to_csv("../../data/mhcii_tcr_template_seqs.csv")

    meta_df = pd.read_csv(
        "../../data/mhcii_tcr_template_metadata.csv", index_col="pdb_id"
    )
    records = list()

    for pdb_id in tqdm(pdb_ids):
        record = dict(pdb_id=pdb_id)
        record["peptide_seq"] = seq_df.loc[pdb_id].P
        record["peptide_len"] = len(record["peptide_seq"])
        record["species"] = meta_df.loc[pdb_id].species
        record["allotype"] = meta_df.loc[pdb_id].allotype
        template = TCRpMHCII(f"../../data/mhcii_tcr_templates/{pdb_id}.pdb")
        template.pmhc.estimate_anchors()
        anchors = template.pmhc.anchors
        record["anchors"] = ";".join([str(anch) for anch in anchors])
        for i in range(len(anchors)):
            record[f"anch_{i + 1}"] = anchors[i]
        record["left_pfr_len"] = anchors[0] - 1
        record["right_pfr_len"] = record["peptide_len"] - anchors[-1]
        record["resolution"] = meta_df.loc[pdb_id].resolution
        records.append(record)

    df = pd.DataFrame.from_dict(records).set_index("pdb_id")
    df.to_csv("../../data/mhcii_tcr_templates.csv")
