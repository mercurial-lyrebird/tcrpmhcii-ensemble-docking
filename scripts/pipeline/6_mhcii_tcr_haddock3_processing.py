
import os
import shutil
import sys
from tqdm import tqdm

from pymol import cmd, stored

from pmhclib.components import Complex, Chain
from pmhclib.components.pmhc import TCRpMHCII
from pmhclib.modeling import Ensemble


def clip(chain_seq: str, full_seq: str, tol: int = 10):

    for i_p in range(0, tol):
        for j_p in range(len(chain_seq), len(chain_seq) - tol, -1):
            if full_seq.find(chain_seq[i_p:j_p]) != -1:
                return chain_seq[i_p:j_p]


def rechain_haddock3_models(model: Complex, reference: TCRpMHCII) -> TCRpMHCII:

    mhc_a_seq = str(reference.pmhc.mhc.alpha_chain.sequence)
    mhc_b_seq = str(reference.pmhc.mhc.beta_chain.sequence)
    pep_seq = str(reference.pmhc.peptide.sequence)
    tcr_a_seq = str(reference.tcr.alpha_chain.sequence)
    tcr_b_seq = str(reference.tcr.beta_chain.sequence)

    l_chain = Chain.from_complex(model, chain_id="L")
    r_chain = Chain.from_complex(model, chain_id="R")
    l_seq = str(l_chain.sequence)
    r_seq = str(r_chain.sequence)

    mhc_a_seq_clipped = clip(mhc_a_seq, r_seq)
    mhc_b_seq_clipped = clip(mhc_b_seq, r_seq)
    pep_seq_clipped = clip(pep_seq, r_seq)
    tcr_a_seq_clipped = clip(tcr_a_seq, l_seq)
    tcr_b_seq_clipped = clip(tcr_b_seq, l_seq)

    mhc_a_seq_i = r_seq.find(mhc_a_seq_clipped)
    mhc_b_seq_i = r_seq.find(mhc_b_seq_clipped)
    pep_seq_i = r_seq.find(pep_seq_clipped)
    tcr_a_seq_i = l_seq.find(tcr_a_seq_clipped)
    tcr_b_seq_i = l_seq.find(tcr_b_seq_clipped)

    cmd.load(model.fp, model.id)
    for chain_id in ["L", "R"]:
        stored.resi = 1
        stored.dict = dict()
        cmd.iterate(f"(chain {chain_id} & name CA)", "stored.dict[(chain, resi, resn)] = stored.resi; stored.resi += 1")
        cmd.alter(f"chain {chain_id}", "resi=str(stored.dict[(chain, resi, resn)])")

    cmd.alter(f"(chain R & resi {mhc_a_seq_i + 1}-{mhc_a_seq_i + len(mhc_a_seq_clipped)})", "chain='M'")
    cmd.alter(f"(chain R & resi {mhc_b_seq_i + 1}-{mhc_b_seq_i + len(mhc_b_seq_clipped)})", "chain='N'")
    cmd.alter(f"(chain R & resi {pep_seq_i + 1}-{pep_seq_i + len(pep_seq_clipped)})", "chain='P'")
    cmd.alter(f"(chain L & resi {tcr_a_seq_i + 1}-{tcr_a_seq_i + len(tcr_a_seq_clipped)})", "chain='A'")
    cmd.alter(f"(chain L & resi {tcr_b_seq_i + 1}-{tcr_b_seq_i + len(tcr_b_seq_clipped)})", "chain='B'")

    for chain_id in ["M", "N", "P", "A", "B"]:
        stored.resi = 1
        stored.dict = dict()
        cmd.iterate(f"(chain {chain_id} & name CA)", "stored.dict[(chain, resi, resn)] = stored.resi; stored.resi += 1")
        cmd.alter(f"chain {chain_id}", "resi=str(stored.dict[(chain, resi, resn)])")

    cmd.save(model.fp, model.id)
    cmd.delete(model.id)

    l_chain.delete()
    r_chain.delete()

    return TCRpMHCII(model.fp)


if __name__ == "__main__":

    pdb_id = sys.argv[1]
    template_pdb_fp = os.path.join(
        "../../data/mhcii_tcr_templates/", f"{pdb_id}.pdb"
    )
    haddock3_dir = sys.argv[2]
    dest_dir = sys.argv[3]

    os.makedirs(dest_dir, exist_ok=True)

    template = TCRpMHCII(template_pdb_fp)
    models = Ensemble.from_pdb_dir(pdb_id, haddock3_dir)

    for id in tqdm(models.ids):
        model = Complex(models.id_to_fp[id], id=id)
        model = rechain_haddock3_models(model, template)
        shutil.copy(model.fp, dest_dir)
        model.delete(rec=True)

    template.delete(rec=True)
