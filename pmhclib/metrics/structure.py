
from pymol import cmd

from ._metric import Metric
from ..components import Structure


class RMSD(Metric):

    def __init__(self, method: str = "align"):

        super(RMSD, self).__init__()

        self.method = method

        return

    def forward(self, struct_1: Structure, struct_2: Structure) -> float:

        cmd.load(struct_1.fp, "struct_1")
        cmd.load(struct_2.fp, "struct_2")

        result = cmd.align("struct_1", "struct_2")

        cmd.delete("struct_1")
        cmd.delete("struct_2")

        return result[0]
