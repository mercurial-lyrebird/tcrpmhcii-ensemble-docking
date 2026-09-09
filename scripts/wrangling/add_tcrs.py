
import json
import os

from pymol import cmd


if __name__ == "__main__":

    with open("../data/mhcii_tcr_chains.json", "rb") as f:
        chain_map = json.load(f)

    for pdb_id, [crystal_a_chain, crystal_b_chain] in chain_map.items():
        cmd.load(os.path.join("../data/mhcii_templates/", f"{pdb_id}.pdb"), "mhcii_template")
        cmd.fetch(pdb_id, "mhcii_tcr_crystal")
        cmd.create("mhcii_tcr_template", f"mhcii_template or (mhcii_tcr_crystal & chain {crystal_a_chain}+{crystal_b_chain})")
        cmd.alter(f"(mhcii_tcr_template & chain {crystal_a_chain})", "chain='A'")
        cmd.alter(f"(mhcii_tcr_template & chain {crystal_b_chain})", "chain='B'")
        print(cmd.get_chains("mhcii_tcr_template"))
        cmd.delete("mhcii_template")
        cmd.delete("mhcii_tcr_crystal")
        cmd.save(f"../data/mhcii_tcr_templates_raw/{pdb_id}.pdb", "mhcii_tcr_template")
        cmd.delete("mhcii_tcr_template")
        os.remove(f"{pdb_id.lower()}.cif")
        print("-----")
