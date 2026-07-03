
from __future__ import annotations

import shutil
import os
from typing import List, Union

import numpy as np

from Bio.PDB import PDBParser
from biopandas.pdb import PandasPdb
from freesasa import Structure as FSStructure
from freesasa import calc
from pymol import cmd, stored

from . import Sequence
from ..utils import get_pdb_id_from_fp
from ..constants import DEFAULT_PDB_STORE


ATOM_MASSES = {"H": 1.008, "C": 12.01, "N": 14.01, "O": 16.00, "S": 32.06, "P": 30.97}


class Structure:

    def __init__(self, fp: str, id: str = None, copy: bool = True):

        if id is None:
            id = os.path.basename(fp)[:-4]

        if copy and DEFAULT_PDB_STORE not in fp:
            dest_fp = os.path.join(DEFAULT_PDB_STORE, f"{id}.pdb")
            shutil.copy(fp, dest_fp)
            fp = dest_fp

        self.fp = fp
        self.id = id

        return

    @classmethod
    def from_cif(cls, cif_fp: str, **kwargs):

        cmd.load(cif_fp, "tmp")
        pdb_fp = os.path.join(
            DEFAULT_PDB_STORE, os.path.basename(cif_fp).replace(".cif", ".pdb")
        )
        cmd.save(pdb_fp, "tmp")
        cmd.delete("tmp")

        return cls(pdb_fp, **kwargs)

    @classmethod
    def from_sdf(cls, sdf_fp: str, **kwargs):

        cmd.load(sdf_fp, "tmp")
        pdb_fp = os.path.join(
            DEFAULT_PDB_STORE, os.path.basename(sdf_fp).replace(".sdf", ".pdb")
        )
        cmd.save(pdb_fp, "tmp")
        cmd.delete("tmp")

        return cls(pdb_fp, **kwargs)

    def select(self, selection: Union[str, List[str]]):

        if isinstance(selection, str):
            selection = [selection]

        cmd.load(self.fp, self.id)
        selection_id = "__".join(
            [self.id] + [
                "_".join(
                    [x.replace(" ", "") for x in y.split(" & ")]
                )
                for y in selection
            ]
        )
        cmd.create(selection_id, " | ".join(selection))
        pdb_fp = os.path.join(DEFAULT_PDB_STORE, f"{selection_id}.pdb")
        cmd.save(pdb_fp, selection_id)
        cmd.delete(selection_id)
        cmd.delete(self.id)

        return Structure(pdb_fp)

    def get_coordinates(self, ca_only: bool = False) -> np.ndarray:

        ppdb = PandasPdb().read_pdb(self.fp)
        df = ppdb.df["ATOM"]
        if ca_only:
            df = df[df["atom_name"] == "CA"]

        return df[["x_coord", "y_coord", "z_coord"]].to_numpy()

    def center_of_mass(self) -> np.ndarray:

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(self.id, self.fp)

        total_mass = 0
        x_sum, y_sum, z_sum = 0, 0, 0

        for atom in structure.get_atoms():
            element = atom.element.strip()
            mass = ATOM_MASSES.get(element, 0)
            if mass == 0:
                continue
            x, y, z = atom.coord
            x_sum += mass * x
            y_sum += mass * y
            z_sum += mass * z
            total_mass += mass

        x_com = x_sum / total_mass
        y_com = y_sum / total_mass
        z_com = z_sum / total_mass

        return np.array([x_com, y_com, z_com])

    def sasa(self, chain_id: str, res_i: List = None) -> np.ndarray:

        structure = FSStructure(self.fp)
        result = calc(structure).residueAreas()[chain_id]

        if res_i is None:
            res_i = range(1, len(result) + 1)

        return np.array([result[str(i)].total for i in res_i])

    def rmsd(self, ref_pdb_fp: str, selection: str = "all") -> float:

        cmd.load(self.fp, self.id)
        cmd.load(ref_pdb_fp, "ref")

        rmsd = cmd.align(f"{self.id} & {selection}", f"ref & {selection}")[0]

        cmd.delete(self.id)
        cmd.delete("ref")

        return rmsd

    def delete(self):

        os.remove(self.fp)


class Complex(Structure):

    @classmethod
    def from_complex(
        cls, complex: Complex, chain_ids: List[str], **kwargs
    ) -> Complex:

        cmd.load(complex.fp, complex.id)
        chain_ids_str = "".join(chain_ids)
        id = f"{complex.id}_chains{chain_ids_str}"
        selection = " or ".join([f"chain {chain_id}" for chain_id in chain_ids])
        cmd.create(id, f"{complex.id} & {selection}")
        pdb_fp = os.path.join(DEFAULT_PDB_STORE, f"{id}.pdb")
        cmd.save(pdb_fp, id)
        cmd.delete(id)
        cmd.delete(complex.id)

        return cls(pdb_fp, id=id, **kwargs)


class Multimer(Complex):

    def merge(
        self, new_chain_id: str, add_bonds: bool = False, fmt: str = "pdb",
        **kwargs
    ) -> Chain:

        if add_bonds:
            raise NotImplementedError

        merged_id = f"{self.id}__merged_chain{new_chain_id}"
        merged_pdb_fp = os.path.join(DEFAULT_PDB_STORE, f"{merged_id}.pdb")
        os.system(f"pdb_merge {self.fp} | pdb_chain -{new_chain_id} | pdb_tidy -strict | pdb_reres > {merged_pdb_fp}")

        return Chain(merged_pdb_fp, **kwargs)


class Chain(Structure):

    def __init__(self, fp: str, id: str = None):

        super().__init__(fp, id=id)

        self.sequence = Sequence(seq=self._load_sequence_str())

        return

    def _load_sequence_str(self):

        cmd.load(self.fp, self.id)
        seq_str = "".join(cmd.get_fastastr(self.id).split("\n")[1:])
        cmd.delete(self.id)

        return seq_str

    def __len__(self) -> int:

        return len(self.sequence)

    def __getitem__(self, idx: int) -> Union[str, List[str]]:

        return self.sequence[idx]

    @classmethod
    def from_multimer(
        cls, multimer: Multimer, chain_id: str, **kwargs
    ) -> Chain:

        cmd.load(multimer.fp, multimer.id)
        id = f"{multimer.id}_chain{chain_id}"
        cmd.create(id, f"{multimer.id} & chain {chain_id}")
        pdb_fp = os.path.join(DEFAULT_PDB_STORE, f"{id}.pdb")
        cmd.save(pdb_fp, id)
        cmd.delete(id)
        cmd.delete(multimer.id)

        return cls(pdb_fp, id=id, **kwargs)

    @classmethod
    def from_complex(cls, complex: Complex, chain_id: str, **kwargs) -> Chain:

        cmd.load(complex.fp, complex.id)
        id = f"{complex.id}_chain{chain_id}"
        cmd.create(id, f"{complex.id} & chain {chain_id}")
        pdb_fp = os.path.join(DEFAULT_PDB_STORE, f"{id}.pdb")
        cmd.save(pdb_fp, id)
        cmd.delete(id)
        cmd.delete(complex.id)

        return cls(pdb_fp, id=id, **kwargs)

    def subchain(self, i: int, j: int) -> Chain:

        cmd.load(self.fp, self.id)
        id = f"{self.id}_res{i}-{j}"
        cmd.create(id, f"{self.id} & resi {i}-{j}")
        pdb_fp = os.path.join(DEFAULT_PDB_STORE, f"{id}.pdb")
        cmd.save(pdb_fp, id)
        cmd.delete(id)
        cmd.delete(self.id)

        return Chain(pdb_fp, id=id)
