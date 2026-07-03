
from typing import Union

from scipy.spatial.distance import hamming
from torch import nn, Tensor

from Bio import Align

from ._metric import Metric


BLOSUM62 = Align.substitution_matrices.load("BLOSUM62")
BLOSUM62_ALPHABET = BLOSUM62.alphabet


class Tokenizer(nn.Module):

    def __init__(self, alphabet: str = BLOSUM62_ALPHABET):

        super(Tokenizer, self).__init__()

        self.alphabet = alphabet

        return

    def forward(self, input: str) -> Tensor:

        return Tensor([self.alphabet.find(c) for c in input])


class BLOSUM62Score(Metric):

    def __init__(self):

        super(BLOSUM62Score, self).__init__()

        return

    def forward(self, input_1: str, input_2: str) -> float:

        return sum(
            [2 ** (-0.5 * BLOSUM62[input_1[i], input_2[i]])
             for i in range(len(input_1))]
        )


class HammingDistance(Metric):

    def __init__(self):

        super(HammingDistance, self).__init__()

        self.tokenizer = Tokenizer()

        return

    def forward(
        self, input_1: Union[str, Tensor], input_2: Union[str, Tensor]
    ) -> Tensor:

        if isinstance(input_1, str):
            input_1 = self.tokenizer(input_1)
        if isinstance(input_2, str):
            input_2 = self.tokenizer(input_2)

        return hamming(input_1, input_2)
