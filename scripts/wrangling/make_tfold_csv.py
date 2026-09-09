
import sys

from pymol import cmd


if __name__ == "__main__":

    pdb_ids = sys.argv[1:]

    f = open("models.csv", "w")
    f.write("pep,MHC sequence,class,pmhc_id\n")

    for pdb_id in pdb_ids:
        cmd.load(f"../data/mhcii_templates/{pdb_id}.pdb", pdb_id)
        fasta_str = cmd.get_fastastr(pdb_id)
        chain_seqs = dict()
        for segment in fasta_str.split(">")[1:]:
            chain_id, subseqs = segment.split("\n")[0], segment.split("\n")[1:]
            chain_seqs[chain_id.split("_")[-1]] = "".join(subseqs)
        f.write(",".join([chain_seqs["P"], chain_seqs["M"] + "/" + chain_seqs["N"], "II", pdb_id]) + "\n")

    f.close()
