
from pymol import cmd

from ..components.pmhc.structure.peptide import MHCPeptide


def rmsd(peptide: MHCPeptide, peptide_model: MHCPeptide, reindex_model: bool = False):

    cmd.load(peptide.fp, "peptide_actual")
    cmd.load(peptide_model.fp, "peptide_model")

    actual_chain_id = cmd.get_chains("peptide_actual")[0]
    cmd.alter("peptide_model", f"chain=\"{actual_chain_id}\"")
    cmd.alter("all", "segi=\"\"")
    if reindex_model:
        cmd.alter("peptide_model", "resi=str(int(resi) + 1)")

    cmd.create("peptide_actual_backbone", "peptide_actual & name C+CA+N+O")
    cmd.create("peptide_model_backbone", "peptide_model & name C+CA+N+O")

    core_selection = f"resi {min(peptide.anchors)}-{max(peptide.anchors)}"
    lrmsd = cmd.rms("peptide_model_backbone", "peptide_actual_backbone")
    lrmsd_core = cmd.rms(
        f"peptide_model_backbone & {core_selection}",
        f"peptide_actual_backbone & {core_selection}",
    )

    rmsd = cmd.rms("peptide_model", "peptide_actual")
    rmsd_core = cmd.rms(
        f"peptide_model & {core_selection}",
        f"peptide_actual & {core_selection}"
    )

    print(f"L-RMSD: {lrmsd}")
    print(f"L-RMSD (core): {lrmsd_core}")
    print(f"RMSD: {rmsd}")
    print(f"RMSD (core): {rmsd_core}")

    return
