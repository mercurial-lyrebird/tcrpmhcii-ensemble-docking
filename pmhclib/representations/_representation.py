
import abc
from typing import Any

from torch import nn, Tensor


class Representation(nn.Module):

    def __init__(self):

        super(Representation, self).__init__()

        return

    @abc.abstractmethod
    def forward(self, input: Any) -> Tensor:

        raise NotImplementedError
