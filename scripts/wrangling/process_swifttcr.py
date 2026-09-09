
from glob import glob
import os
import sys
from tqdm import tqdm

from pymol import cmd, stored

from pmhclib.components.pmhc import TCRpMHCII
from pmhclib.components import Complex, Chain


def clip(chain_seq: str, full_seq: str, tol: int = 10):

    for i_p in range(0, tol):
        for j_p in range(len(chain_seq), len(chain_seq) - tol, -1):
            if full_seq.find(chain_seq[i_p:j_p]) != -1:
                return chain_seq[i_p:j_p]


if __name__ == "__main__":

    template_pdb_fp = sys.argv[1]
    swifttcr_dir = sys.argv[2]
    dest_dir = sys.argv[3]

    os.makedirs(dest_dir, exist_ok=True)

    template = TCRpMHCII(template_pdb_fp)
    mhc_a_seq = str(template.pmhc.mhc.alpha_chain.sequence)
    mhc_b_seq = str(template.pmhc.mhc.beta_chain.sequence)
    pep_seq = str(template.pmhc.peptide.sequence)
    tcr_a_seq = str(template.tcr.alpha_chain.sequence)
    tcr_b_seq = str(template.tcr.beta_chain.sequence)

    centroid_fps = list()
    with open(os.path.join(swifttcr_dir, "clustering.txt"), "r") as f:
        for line in f.readlines():
            centroid_fn = line.strip().split()[2]
            if not centroid_fn.endswith(".pdb"):
                continue
            centroid_fps.append(os.path.join(swifttcr_dir, "merged", centroid_fn))

    for pdb_fp in tqdm(centroid_fps):
        model = Complex(pdb_fp)
        l_chain = Chain.from_complex(model, chain_id="A")
        r_chain = Chain.from_complex(model, chain_id="D")
        l_seq = str(l_chain.sequence)
        r_seq = str(r_chain.sequence)
        mhc_a_seq_clipped = clip(mhc_a_seq, l_seq)
        mhc_b_seq_clipped = clip(mhc_b_seq, l_seq)
        pep_seq_clipped = clip(pep_seq, l_seq)
        tcr_a_seq_clipped = clip(tcr_a_seq, r_seq)
        tcr_b_seq_clipped = clip(tcr_b_seq, r_seq)
        mhc_a_seq_i = l_seq.find(mhc_a_seq_clipped)
        mhc_b_seq_i = l_seq.find(mhc_b_seq_clipped)
        pep_seq_i = l_seq.find(pep_seq_clipped)
        tcr_a_seq_i = r_seq.find(tcr_a_seq_clipped)
        tcr_b_seq_i = r_seq.find(tcr_b_seq_clipped)
        cmd.load(model.fp, model.id)
        cmd.alter("(chain A)", "chain='L'")
        cmd.alter("(chain D)", "chain='R'")
        for chain_id in ["L", "R"]:
            stored.resi = 1
            stored.dict = dict()
            cmd.iterate(f"(chain {chain_id} & name CA)", "stored.dict[(chain, resi, resn)] = stored.resi; stored.resi += 1")
            cmd.alter(f"chain {chain_id}", "resi=str(stored.dict[(chain, resi, resn)])")
        print(f"{mhc_a_seq_i + 1}-{mhc_a_seq_i + len(mhc_a_seq_clipped)})")
        print(f"{mhc_b_seq_i + 1}-{mhc_b_seq_i + len(mhc_b_seq_clipped)}")
        print(f"{pep_seq_i + 1}-{pep_seq_i + len(pep_seq_clipped)}")
        cmd.alter(f"(chain L & resi {mhc_a_seq_i + 1}-{mhc_a_seq_i + len(mhc_a_seq_clipped)})", "chain='M'")
        cmd.alter(f"(chain L & resi {mhc_b_seq_i + 1}-{mhc_b_seq_i + len(mhc_b_seq_clipped)})", "chain='N'")
        cmd.alter(f"(chain L & resi {pep_seq_i + 1}-{pep_seq_i + len(pep_seq_clipped)})", "chain='P'")
        cmd.alter(f"(chain R & resi {tcr_a_seq_i + 1}-{tcr_a_seq_i + len(tcr_a_seq_clipped)})", "chain='A'")
        cmd.alter(f"(chain R & resi {tcr_b_seq_i + 1}-{tcr_b_seq_i + len(tcr_b_seq_clipped)})", "chain='B'")
        for chain_id in ["M", "N", "P", "A", "B"]:
            stored.resi = 1
            stored.dict = dict()
            cmd.iterate(f"(chain {chain_id} & name CA)", "stored.dict[(chain, resi, resn)] = stored.resi; stored.resi += 1")
            cmd.alter(f"chain {chain_id}", "resi=str(stored.dict[(chain, resi, resn)])")
        cmd.save(os.path.join(dest_dir, os.path.basename(model.fp)), model.id)
        cmd.delete(model.id)
