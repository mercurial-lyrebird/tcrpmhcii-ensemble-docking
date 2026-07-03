
import abc
import math
import os
import re
import shutil
import subprocess
from typing import Dict, List, Tuple, Union

from scipy.spatial.distance import euclidean

import pyrosetta
pyrosetta.init("-out:levels basic:error core:error")

from Bio.PDB import Structure as BioPDBStructure
from DockQ.DockQ import load_PDB, run_on_all_native_interfaces
from pymol import cmd
from pyrosetta import pose_from_pdb, get_fa_scorefxn, Pose, ScoreFunction

from ...components.structure import Structure, Complex, Chain, Multimer


class Scorer(abc.ABC):

    def __init__(self, return_terms: bool = False):

        self.return_terms = return_terms

        return

    def __call__(self, structure: Structure, **kwargs) -> float:

        return self.score(structure, **kwargs)

    @abc.abstractmethod
    def score(self, structure: Structure, **kwargs) -> float:

        raise NotImplementedError


class MinChainDistScorer(Scorer):

    name = "min_chain_dist"

    def __init__(self, chains: Tuple[str], **kwargs):

        super(MinChainDistScorer, self).__init__(**kwargs)

        self.chains = chains
        self.name += "".join(chains)

        return

    def score(self, structure: Structure) -> float:

        (chain_1_id, chain_2_id) = self.chains

        cmd.load(structure.fp, structure.id)
        cmd.select("chain_1", f"{structure.id} & chain {chain_1_id} & name CA")
        chain_1 = cmd.get_model("chain_1")
        cmd.select("chain_2", f"{structure.id} & chain {chain_2_id} & name CA")
        chain_2 = cmd.get_model("chain_2")

        min_dist = math.inf
        for atom_1 in chain_1.atom:
            for atom_2 in chain_2.atom:
                dist = euclidean(atom_1.coord, atom_2.coord)
                if dist < min_dist:
                    min_dist = dist

        return min_dist


class DockQScorer(Scorer):

    name = "dockq"

    def __init__(
        self, reference: Structure, chain_config: List[Union[str, Tuple[str]]],
        f_nat: bool = False, **kwargs
    ):

        super(DockQScorer, self).__init__(**kwargs)

        self.reference = reference
        self.chain_config = chain_config
        self.f_nat = f_nat

        return

    def score(self, structure: Structure) -> float:

        reference = load_PDB(self.reference.fp)
        model = load_PDB(structure.fp)

        chain_map = dict()
        chain_ids = list()
        for obj in self.chain_config:
            if isinstance(obj, Tuple):
                new_chain_id = "".join(obj)
                reference = self.__merge_chains(reference, obj, new_chain_id)
                model = self.__merge_chains(model, obj, new_chain_id)
                chain_id = new_chain_id
            elif isinstance(obj, str):
                chain_id = obj
            else:
                raise TypeError
            chain_map[chain_id] = chain_id
            chain_ids.append(chain_id)

        result = run_on_all_native_interfaces(
            model, reference, chain_map=chain_map
        )

        if self.f_nat:
            score = result[0]["".join(chain_ids)]["fnat"]
        else:
            score = result[0]["".join(chain_ids)]["DockQ"]

        return score

    @staticmethod
    def __merge_chains(
        structure: BioPDBStructure, chains: Tuple[str], new_chain_id: str
    ) -> BioPDBStructure:

        chain_1 = chains[0]
        for chain_2 in chains[1:]:
            for i, res in enumerate(structure[chain_2]):
                res.id = (chain_2, res.id[1], res.id[2])
                structure[chain_1].add(res)
            structure.detach_child(chain_2)
        structure[chain_1].id = new_chain_id

        return structure


class RMSDScorer(Scorer):

    name = "rmsd"

    def __init__(self, reference: Structure, selection: str = "all", **kwargs):

        super(RMSDScorer, self).__init__(**kwargs)

        self.reference = reference
        self.selection = selection

        return

    def score(self, structure: Structure) -> float:

        return structure.rmsd(self.reference.fp, selection=self.selection)


class RosettaScorer(Scorer):

    name = "rosetta"

    def __init__(self, **kwargs):

        super(RosettaScorer, self).__init__(**kwargs)

        self.weights = None

        return

    def set_weights(self, weights: List[float]):

        self.weights = weights

        return

    def reset_weights(self):

        self.weights = None

        return

    def score(self, structure: Structure) -> float:

        pose = pose_from_pdb(structure.fp)
        scoring_fxn = get_fa_scorefxn()
        terms = scoring_fxn.get_nonzero_weighted_scoretypes()

        if self.weights is not None:
            for term, weight in zip(terms, self.weights):
                scoring_fxn.set_weight(term, weight)

        result = scoring_fxn(pose)

        if self.return_terms:
            result = self.__get_energy_terms(pose, scoring_fxn)

        return result

    def __get_energy_terms(
        self, pose: Pose, scoring_fxn: ScoreFunction
    ) -> List:

        terms = scoring_fxn.get_nonzero_weighted_scoretypes()
        scoring_fxn(pose)

        return [(term, pose.energies().total_energies()[term]) for term in terms]


class ReceptorLigandScorer(Scorer):

    def __init__(
        self, receptor_chain: Union[str, List[str]],
        ligand_chain: Union[str, List[str]], **kwargs
    ):

        super(ReceptorLigandScorer, self).__init__(**kwargs)

        self.receptor_chain = receptor_chain
        self.ligand_chain = ligand_chain

        return

    def split_components(self, complex: Complex) -> Tuple[Structure, Structure]:

        if isinstance(self.receptor_chain, str):
            receptor = Chain.from_complex(complex, chain_id=self.receptor_chain)
        else:
            receptor = Multimer.from_complex(
                complex, chain_ids=self.receptor_chain
            )

        if isinstance(self.ligand_chain, str):
            ligand = Chain.from_complex(complex, chain_id=self.ligand_chain)
        else:
            ligand = Multimer.from_complex(complex, chain_ids=self.ligand_chain)

        return receptor, ligand

    @abc.abstractmethod
    def score_rl(self, receptor: Structure, ligand: Structure) -> float:

        raise NotImplementedError

    def score(self, structure: Structure) -> float:

        receptor, ligand = self.split_components(structure)
        result = self.score_rl(receptor, ligand)
        receptor.delete()
        ligand.delete()

        return result


class SMINAScorer(ReceptorLigandScorer):

    name = "smina"

    def _preprocess(
        self, receptor: Structure, ligand: Structure, ligand_fmt: str = "pdb"
    ) -> Tuple[str, str]:

        receptor_pdbqt_fp = receptor.fp.replace(".pdb", ".pdbqt")
        ligand_pdbqt_fp = ligand.fp.replace(f".{ligand_fmt}", ".pdbqt")

        subprocess.call(
            ["prepare_receptor4.py", "-r", receptor.fp, "-o", receptor_pdbqt_fp]
        )
        shutil.copy(ligand.fp, os.getcwd())
        subprocess.call(
            [
                "prepare_ligand4.py", "-l", os.path.basename(ligand.fp), "-o",
                ligand_pdbqt_fp
            ]
        )
        os.remove(os.path.basename(ligand.fp))

        return (receptor_pdbqt_fp, ligand_pdbqt_fp)

    def score_rl(
        self, receptor: Structure, ligand: Structure, ligand_fmt: str = "pdbqt"
    ) -> float:

        (receptor_pdbqt_fp, ligand_pdbqt_fp) = self._preprocess(
            receptor, ligand, ligand_fmt=ligand_fmt
        )

        result = subprocess.run(
            [
                "smina", "-r", receptor_pdbqt_fp, "-l", ligand_pdbqt_fp,
                "--score_only"
            ],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        output = result.stdout + result.stderr

        if self.return_terms:
            result = self.__extract_terms(output)
        else:
            result = self.__extract_energy(output)

        return result

    @staticmethod
    def __extract_energy(output: str) -> float:

        match = re.search(r"Affinity:\s*([-+]?\d*\.\d+|\d+)", output)
        if match:
            energy = float(match.group(1))
        else:
            raise ValueError("Could not find total energy in smina output.")

        return energy

    @staticmethod
    def __extract_terms(output: str) -> Tuple[List[float], List[float]]:

        # Weights      Terms
        # -0.035579    gauss(o=0,_w=0.5,_c=8)
        # -0.005156    gauss(o=3,_w=2,_c=8)
        # 0.840245     repulsion(o=0,_c=8)
        # -0.035069    hydrophobic(g=0.5,_b=1.5,_c=8)
        # -0.587439    non_dir_h_bond(g=-0.7,_b=0,_c=8)
        # 1.923        num_tors_div

        lines = output.splitlines()
        term_vals = None
        for line in lines:
            if line.startswith("##"):
                xs = line.strip("##").strip().split()
                if any([x.isalpha() for x in xs]):
                    continue
                else:
                    term_vals = [float(x) for x in xs]

        return term_vals


class VinardoScorer(SMINAScorer):

    name = "vinardo"

    def score_rl(
        self, receptor: Structure, ligand: Structure, ligand_fmt: str = "pdbqt"
    ) -> float:

        (receptor_pdbqt_fp, ligand_pdbqt_fp) = self._preprocess(
            receptor, ligand, ligand_fmt=ligand_fmt
        )

        result = subprocess.run(
            [
                "smina", "-r", receptor_pdbqt_fp, "-l", ligand_pdbqt_fp,
                "--scoring", "vinardo", "--score_only"
            ],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        output = result.stdout + result.stderr
        print(output)

        if self.return_terms:
            result = self.__extract_terms(output)
        else:
            result = self.__extract_energy(output)

        return result

    @staticmethod
    def __extract_energy(output: str) -> float:

        match = re.search(r"Affinity:\s*([-+]?\d*\.\d+|\d+)", output)
        if match:
            energy = float(match.group(1))
        else:
            raise ValueError("Could not find total energy in smina output.")

        return energy

    @staticmethod
    def __extract_terms(output: str) -> Tuple[List[float], List[float]]:

        # Weights      Terms
        # -0.045       gauss(o=0,_w=0.8,_c=8)
        # 0.8          repulsion(o=0,_c=8)
        # -0.035       hydrophobic(g=0,_b=2.5,_c=8)
        # -0.6         non_dir_h_bond(g=-0.6,_b=0,_c=8)
        # 0            num_tors_div

        lines = output.splitlines()
        term_vals = None
        for line in lines:
            if line.startswith("##"):
                xs = line.strip("##").strip().split()
                if any([x.isalpha() for x in xs]):
                    continue
                else:
                    term_vals = [float(x) for x in xs]

        return term_vals


class HADDOCK3Scorer(ReceptorLigandScorer):

    name = "haddock3"

    def score_rl(self, receptor: Structure, ligand: Structure) -> float:

        return
