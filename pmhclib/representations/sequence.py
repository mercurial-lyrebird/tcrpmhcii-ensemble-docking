
from tqdm import tqdm
from typing import List, Union

import torch
from torch import Tensor

import esm

from ._representation import Representation
from .constants import ATCHLEY_FACTORS, KIDERA_FACTORS


class Atchley(Representation):

    def __init__(self, flatten: bool = False):

        super(Atchley, self).__init__()

        self.flatten = flatten

        return

    def forward(self, input: str) -> Tensor:

        X = torch.stack([Tensor(ATCHLEY_FACTORS[c]) for c in input])
        if self.flatten:
            X = X.flatten()

        return X


class Kidera(Representation):

    def __init__(self, flatten: bool = False):

        super(Kidera, self).__init__()

        self.flatten = flatten

        return

    def forward(self, input: str) -> Tensor:

        X = torch.stack([Tensor(KIDERA_FACTORS[c]) for c in input])
        if self.flatten:
            X = X.flatten()

        return X


class Embedding(Representation):

    def __init__(self):

        super(Embedding, self).__init__()

        return

    def forward(self, input: Union[str, List[str]]) -> Tensor:

        raise NotImplementedError


class ESMEmbedding(Representation):

    def __init__(self):

        super(ESMEmbedding, self).__init__()

        self.model, self.alphabet = esm.pretrained.esm2_t6_8M_UR50D()
        self.batch_converter = self.alphabet.get_batch_converter()
        self.model.eval()

        self.batch_size = 8

        return

    def forward(self, input: Union[str, List[str]]) -> Tensor:

        if isinstance(input, str):
            input = [input]

        data = [(i, seq) for i, seq in enumerate(input)]
        for j in tqdm(range(len(input) // self.batch_size + 1)):
            _, _, batch_tokens = self.batch_converter(
                data[j:(j + self.batch_size)]
            )
            with torch.no_grad():
                results = self.model(
                    batch_tokens, repr_layers=[6], return_contacts=True
                )
            x = results["representations"][6].mean(axis=1)
            if j == 0:
                X = x
            else:
                X = torch.cat((X, x), dim=0)

        return X[:len(input)]
