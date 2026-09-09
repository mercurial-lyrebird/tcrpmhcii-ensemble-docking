
import os
import sys

import pandas as pd

import MDAnalysis as mda
from pymol import cmd


if __name__ == "__main__":

    [pdb_id, topology_file, trajectory_file, output_dir] = sys.argv[1:]

    print(f"Loading topology: {topology_file}")
    print(f"Loading trajectory: {trajectory_file}")
    u = mda.Universe(topology_file, trajectory_file)
    protein = u.select_atoms("protein")

    df = pd.read_csv("../../data/mhcii_templates.csv", index_col="pdb_id")
    record = df.loc[pdb_id]

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    print(f"Writing frames to {output_dir}")
    for ts in u.trajectory:
        frame_number = ts.frame
        filename = os.path.join(output_dir, f"frame_{frame_number:05d}.pdb")
        protein.write(filename)
        print(f"Written: {filename}")
        cmd.load(filename, f"frame_{frame_number}")
        cmd.alter(f"frame_{frame_number} & resi {protein.n_residues - len(record.peptide_seq) + 1}-{protein.n_residues}", "chain='P'")
        cmd.save(filename, f"frame_{frame_number}")

    print("Done.")