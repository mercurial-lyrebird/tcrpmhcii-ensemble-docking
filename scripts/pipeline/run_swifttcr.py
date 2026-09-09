
import os
import subprocess
import sys


if __name__ == "__main__":

    pmhc_fp = os.path.abspath(sys.argv[1])
    tcr_fp = os.path.abspath(sys.argv[2])
    swifttcr_path = os.path.abspath(sys.argv[3])

    [pdb_id, method] = pmhc_fp.split("/")[-3:-1]
    id = "_".join([pdb_id, method, os.path.basename(pmhc_fp)[:-4]])

    os.chdir(swifttcr_path)
    subprocess.run(
        ["python3", "scripts/swift_tcr.py", "-r", pmhc_fp, "-l", tcr_fp, "-m",
         "20", "-mhc2", "-c", "8", "-o", "results", "-op", id]
    )
