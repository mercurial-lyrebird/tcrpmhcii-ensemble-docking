
from __future__ import annotations

from typing import List, Union


class Sequence:

    def __init__(self, seq: str):

        self.sequence = seq

        return

    def __len__(self) -> int:

        return len(self.sequence)

    def __str__(self) -> str:

        return self.sequence

    def __getitem__(self, idx: Union[int, slice]) -> Union[str, List[str]]:

        if isinstance(idx, slice):
            idx = slice(idx.start, idx.stop + 1, idx.step)

        return self.sequence[idx]

    def find(self, s: Union[str, Sequence]) -> int:

        if isinstance(s, Sequence):
            s = s.sequence

        idx = self.sequence.find(s)

        return idx, idx + len(s) - 1
