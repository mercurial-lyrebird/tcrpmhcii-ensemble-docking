
import abc
import os
from typing import List, Tuple, Union

import numpy as np

from MDAnalysis import Universe
from MDAnalysis.analysis.distances import distance_array

from pymol import cmd
from pymol.cgo import CYLINDER, CONE

from ...structure import Complex, Multimer
from . import pMHCI, pMHCII, TCR


class TCRpMHC(abc.ABC, Complex):

    def __init__(
        self, fp: str, id: str = None,
        tcr_chain_ids: Tuple[str, str] = ("D", "E")
    ):

        super().__init__(fp, id=id)

        tcr_multimer = Multimer.from_complex(self, tcr_chain_ids)
        self.tcr = TCR(
            tcr_multimer.fp, id=tcr_multimer.id, chain_ids=tcr_chain_ids
        )
        self.pmhc: Union[pMHCI, pMHCII] = None

        return

    def get_peptide_tcr_contacts(
        self, distance: float = 5.0, cdr3_only: bool = False
    ) -> List[int]:

        u = Universe(self.fp)

        (alpha_cdrs, beta_cdrs) = self.tcr.alpha_chain.cdr_loops, self.tcr.beta_chain.cdr_loops

        if cdr3_only:
            alpha_cdrs = [alpha_cdrs[2]]
            beta_cdrs = [beta_cdrs[2]]

        cdr_sel_str = "protein and " + " or ".join([f"((chainID {self.tcr.chain_ids[0]} and resid {alpha_cdr[0]}:{alpha_cdr[1]}) or (chainID {self.tcr.chain_ids[1]} and resid {beta_cdr[0]}:{beta_cdr[1]}))" for (alpha_cdr, beta_cdr) in zip(alpha_cdrs, beta_cdrs)])

        cdr_atoms = u.select_atoms(cdr_sel_str)
        peptide_atoms = u.select_atoms(
            f"protein and chainID {self.pmhc.peptide_chain_id}"
        )

        contacts = [0 for _ in peptide_atoms.residues]
        for i, res in enumerate(peptide_atoms.residues):
            distances = distance_array(res.atoms.positions, cdr_atoms.positions)
            if (distances <= distance).any():
                contacts[i] = 1

        return contacts

    def compute_vectors(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:

        X_peptide = self.pmhc.peptide.get_coordinates(ca_only=True)
        X_binding_core = X_peptide[(self.pmhc.anchors[0] - 1):self.pmhc.anchors[-1]]

        G = X_binding_core[-1] - X_binding_core[0]
        G_hat = G / np.linalg.norm(G)

        N = self.pmhc.peptide.center_of_mass() - self.pmhc.mhc.center_of_mass()
        N_hat = N / np.linalg.norm(N)

        T = self.tcr.beta_chain.center_of_mass() - self.tcr.alpha_chain.center_of_mass()
        T_hat = T / np.linalg.norm(T)

        return (G_hat, N_hat, T_hat)

    def compute_docking_angles(self) -> Tuple[float, float]:

        (G_hat, N_hat, T_hat) = self.compute_vectors()

        T_proj = T_hat - np.dot(T_hat, N_hat) * N_hat
        T_proj_hat = T_proj / np.linalg.norm(T_proj)

        crossing = np.degrees(np.arccos(np.clip(np.dot(T_proj_hat, G_hat), -1.0, 1.0)))
        incident = np.degrees(np.arcsin(np.clip(abs(np.dot(T_hat, N_hat)), -1.0, 1.0)))

        return (crossing, incident)

    def peptide_sasa(self, core_only: bool = False) -> np.ndarray:

        if core_only:
            res_i = range(self.pmhc.anchors[0], self.pmhc.anchors[-1] + 1)
        else:
            res_i = None

        return super().sasa(chain_id=self.pmhc.peptide_chain_id, res_i=res_i)

    def create_pymol_session(self):

        (G_hat, N_hat, T_hat) = self.compute_vectors()

        def draw_vector(name, origin, direction, length=15.0,
                radius=0.4, color=(1,0,0)):

            ox, oy, oz = origin
            dx, dy, dz = direction

            ex, ey, ez = ox + dx*length, oy + dy*length, oz + dz*length
            cx, cy, cz = ox + dx*length*0.85, oy + dy*length*0.85, oz + dz*length*0.85

            obj = [
                CYLINDER,
                ox, oy, oz,
                cx, cy, cz,
                radius,
                *color, *color,
                CONE,
                cx, cy, cz,
                ex, ey, ez,
                radius*1.8, 0.0,
                *color, *color,
                1.0, 0.0
            ]
            cmd.load_cgo(obj, name)

        draw_vector("G", self.pmhc.peptide.center_of_mass(), G_hat, length=25.0, color=(0.2, 0.2, 1.0))
        draw_vector("N", self.pmhc.peptide.center_of_mass(), N_hat, length=15.0, color=(0.7, 0.2, 0.7))
        draw_vector("T", self.pmhc.peptide.center_of_mass(), T_hat, length=25.0, color=(1.0, 0.2, 0.2))

        cmd.save(f"{self.id}.pse")

        return

    def delete(self, rec: bool = False):

        os.remove(self.fp)

        if rec:
            self.tcr.delete(rec=True)
            self.pmhc.delete(rec=True)

        return


class TCRpMHCI(TCRpMHC):

    def __init__(
        self, fp: str, mhc_chain_id: str = "A", peptide_chain_id: str = "C",
        tcr_chain_ids: Tuple[str, str] = ("D", "E"), id: str = None,
        mhc_allotype: str = None, anchors: List[int] = None
    ):

        super().__init__(fp, id=id, tcr_chain_ids=tcr_chain_ids)

        pmhc_complex = Complex.from_complex(
            self, [mhc_chain_id, peptide_chain_id]
        )
        self.pmhc = pMHCI(
            pmhc_complex.fp, mhc_chain_id=mhc_chain_id,
            peptide_chain_id=peptide_chain_id, id=pmhc_complex.id,
            allotype=mhc_allotype, anchors=anchors
        )

        return


class TCRpMHCII(TCRpMHC):

    def __init__(
        self, fp: str, mhc_chain_ids: Tuple[str, str] = ("M", "N"),
        peptide_chain_id: str = "P",
        tcr_chain_ids: Tuple[str, str] = ("A", "B"), id: str = None,
        mhc_allotype: str = None, anchors: List[int] = None
    ):

        super().__init__(fp, id=id, tcr_chain_ids=tcr_chain_ids)

        pmhc_complex = Complex.from_complex(
            self, [mhc_chain_ids[0], mhc_chain_ids[1], peptide_chain_id]
        )
        self.pmhc = pMHCII(
            pmhc_complex.fp, mhc_chain_ids=mhc_chain_ids,
            peptide_chain_id=peptide_chain_id, id=pmhc_complex.id,
            allotype=mhc_allotype, anchors=anchors
        )

        return
