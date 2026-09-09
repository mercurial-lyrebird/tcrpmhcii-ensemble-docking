
from argparse import ArgumentParser
from glob import glob
import os
import shutil

import pandas as pd

from pmhclib.components import Chain
from pmhclib.components.pmhc import TCRpMHCII, pMHCII


PDB_PATTERNS = {
    "annealing": "conf_*.pdb",
    "pandora2": "*.BL0*.pdb"
}


if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("pdb_id", default=None)
    parser.add_argument("-s", "--sampling_method", default=None)
    parser.add_argument("-c", "--model_fp", default=None)
    parser.add_argument("-o", "--output_dir", default=None)
    parser.add_argument("-n", "--n_rigidbody", default=None)
    parser.add_argument("-m", "--n_refine", default=None)
    parser.add_argument("-x", "--model_id", default="model")
    parser.add_argument("-a", "--anchors", default=None)
    arguments = parser.parse_args()

    pdb_id = arguments.pdb_id
    id = pdb_id
    output_dir_prefix = os.path.join("../../haddock3_cdr3/", pdb_id)

    template_df = pd.read_csv("../../data/mhcii_tcr_templates.csv", index_col="pdb_id")
    template_record = template_df.loc[pdb_id]
    if arguments.anchors is None:
        anchors = [int(i) for i in template_record.anchors.split(";")]
    else:
        anchors = [int(i) for i in arguments.anchors.split(",")]

    template_pdb_fp = os.path.join(
        "../../data/mhcii_tcr_templates/", f"{pdb_id}.pdb"
    )
    template = TCRpMHCII(template_pdb_fp, id=pdb_id, anchors=anchors)
    tcr = template.tcr

    if arguments.sampling_method is None:
        if arguments.model_fp is None:
            # crystal-on-crystal docking
            id = f"{id}__crystal"
            pmhc = template.pmhc
            output_dir = os.path.join(output_dir_prefix, "static", "crystal")
        else:
            # model-on-crystal docking
            pmhc = pMHCII(
                arguments.model_fp, id=f"{pdb_id}_{arguments.model_id}",
                anchors=anchors
            )
            id = pmhc.id
            output_dir = os.path.join(
                output_dir_prefix, "static", arguments.model_id
            )
    else:
        # ensemble-on-crystal docking
        sampling_method = arguments.sampling_method
        id = f"{id}__{sampling_method}"
        pmhc = None
        output_dir = os.path.join(
            output_dir_prefix, "ensemble", sampling_method
        )

    if arguments.output_dir is not None:
        output_dir = arguments.output_dir

    os.makedirs(output_dir)
    shutil.copy("dock.cfg", output_dir)
    shutil.copy(tcr.fp, os.path.join(output_dir, "_tcr.pdb"))

    if arguments.sampling_method is not None:
        for i, pdb_fp in enumerate(glob(os.path.join("../../data/mhcii_ensembles/", pdb_id, sampling_method,PDB_PATTERNS[sampling_method]))):
            shutil.copy(pdb_fp, os.path.join(output_dir, f"_pmhc_{i}.pdb"))
            pwd = os.getcwd()
            os.chdir(output_dir)
            os.system(f"pdb_chain -R _pmhc_{i}.pdb | pdb_reres -1 > pmhc_{i}.pdb")
            if i == 0:
                pmhc_model = pMHCII("_pmhc_0.pdb", id=f"{id}__pmhc")
                os.system(f"pdb_chain -L _tcr.pdb | pdb_reres -1 > tcr.pdb")
                pmhc_merged = Chain("pmhc_0.pdb", id=f"{id}__pmhc_merged")
                tcr_merged = Chain("tcr.pdb", id=f"{id}__tcr_merged")
            os.chdir(pwd)
        os.chdir(output_dir)
        os.system(" ".join(["pdb_mkensemble", "pmhc_*.pdb", ">", "pmhc.pdb"]))
        pmhc = template.pmhc
    else:
        shutil.copy(pmhc.fp, os.path.join(output_dir, "_pmhc.pdb"))
        pwd = os.getcwd()
        os.chdir(output_dir)
        pmhc_model = pMHCII("_pmhc.pdb", id=f"{id}__pmhc")
        os.system("pdb_chain -R _pmhc.pdb | pdb_reres -1 > pmhc.pdb")
        os.system("pdb_chain -L _tcr.pdb | pdb_reres -1 > tcr.pdb")
        pmhc_merged = Chain("pmhc.pdb", id=f"{id}__pmhc_merged")
        tcr_merged = Chain("tcr.pdb", id=f"{id}__tcr_merged")

    print(pmhc_model.peptide.sequence)
    print(len(pmhc_merged))
    print(pmhc_merged.sequence)
    print(len(tcr_merged))
    print(tcr_merged.sequence)

    print(pmhc_merged.sequence.find(pmhc_model.mhc.alpha_chain.sequence))
    print(pmhc_merged.sequence.find(pmhc_model.mhc.beta_chain.sequence))
    print(pmhc_merged.sequence.find(pmhc_model.peptide.sequence))
    print(tcr.alpha_chain.cdr_loops)
    print(tcr.beta_chain.cdr_loops)

    print(pmhc.peptide.anchors)

    binding_core_seq = pmhc.peptide.sequence[
        (pmhc.peptide.anchors[0] - 1):(pmhc.peptide.anchors[-1] - 1)
    ]
    print(binding_core_seq)
    distal_anchors = pmhc_merged.sequence.find(binding_core_seq)
    print("weewoo", distal_anchors)

    with open("int_restraints.tbl", "w") as f:

        for j in range(2, 3):

            # distal_anchors[i] + 3 to account for 1-indexing and TER residues

            alpha_cdr_resi = tcr.alpha_chain.cdr_loops[j]
            alpha_cdr_seq = tcr.alpha_chain.sequence[alpha_cdr_resi[0]:alpha_cdr_resi[1]]
            alpha_cdr = tcr_merged.sequence.find(alpha_cdr_seq)
            f.write("assign ")
            f.write("(\n")
            for i in range(alpha_cdr[0], alpha_cdr[1]):
                f.write(f"\t(segid L and resi {i}) or\n")
            f.write(f"\t(segid L and resi {alpha_cdr[1]})\n")
            f.write(")\n")
            f.write("(\n")
            for i in range(distal_anchors[0] + 3, distal_anchors[1] + 3):
                f.write(f"\t(segid R and resi {i}) or\n")
            f.write(f"\t(segid R and resi {distal_anchors[1] + 3})\n")
            f.write(")\n")
            f.write("3.0 3.0 0.0\n")

            beta_cdr_resi = tcr.beta_chain.cdr_loops[j]
            beta_cdr_seq = tcr.beta_chain.sequence[beta_cdr_resi[0]:beta_cdr_resi[1]]
            beta_cdr = tcr_merged.sequence.find(beta_cdr_seq)
            f.write("assign ")
            f.write("(\n")
            for i in range(beta_cdr[0], beta_cdr[1]):
                f.write(f"\t(segid L and resi {i}) or\n")
            f.write(f"\t(segid L and resi {beta_cdr[1]})\n")
            f.write(")\n")
            f.write("(\n")
            for i in range(distal_anchors[0] + 3, distal_anchors[1] + 3):
                f.write(f"\t(segid R and resi {i}) or\n")
            f.write(f"\t(segid R and resi {distal_anchors[1] + 3})\n")
            f.write(")\n")
            f.write("3.0 3.0 0.0\n")

    tcr_restr_fp = "tcr_restraints.tbl"
    os.system(f"haddock3-restraints restrain_bodies tcr.pdb > {tcr_restr_fp}")

    pmhc_restr_fp = "pmhc_restraints.tbl"
    os.system(f"haddock3-restraints restrain_bodies pmhc.pdb > {pmhc_restr_fp}")

    chain_restr_fp = "chain_restraints.tbl"
    with open(chain_restr_fp, "w") as f_out:
        for restr_fp in [tcr_restr_fp, pmhc_restr_fp]:
            with open(restr_fp, "r") as f_in:
                shutil.copyfileobj(f_in, f_out)
