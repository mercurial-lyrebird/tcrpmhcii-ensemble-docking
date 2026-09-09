
import pandas as pd

from pymol import cmd

from pmhclib.components.pmhc import pMHC, MHCIIPeptide


if __name__ == "__main__":

    templates = pd.read_csv("../data/mhcii_templates.csv", index_col="pdb_id")
    templates.allotype.fillna("", inplace=True)

    with_tcr = list()
    without_tcr = list()
    for id, template in templates.iterrows():
        pmhc = pMHC(fp=f"../data/mhcii_templates/{id}.pdb", allotype=template.allotype.split(";"))
        peptide = MHCIIPeptide.from_pmhc(pmhc, chain_id="P")
        peptide.set_anchors([int(i) for i in template.anchors.split(";")])
        cmd.load(peptide.fp, peptide.id)
        cmd.label(f"{peptide.id} & resi {peptide.anchors[0]}-{peptide.anchors[-1]}", "\"core\"")
        if template.with_tcr:
            with_tcr.append(peptide.id)
        else:
            without_tcr.append(peptide.id)

    cmd.group("peptides", " ".join(with_tcr + without_tcr))
    cmd.create("ref", with_tcr[0])
    cmd.extra_fit("peptides & label \"core\"", "ref & label \"core\"", "super")
    cmd.label("label \"core\"", "")
    cmd.ungroup("peptides")
    cmd.group("with_tcr", " ".join(with_tcr))
    cmd.group("without_tcr", " ".join(without_tcr))
    cmd.show("cartoon", "with_tcr")
    cmd.show("cartoon", "without_tcr")
    cmd.spectrum("b", "green orange")
    cmd.save("mhcii_peptides_aligned.pse")
