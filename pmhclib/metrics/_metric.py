
import abc
from typing import Any

from torch import nn


class Metric(nn.Module):

    def __init__(self):

        super(Metric, self).__init__()

        return

    @abc.abstractmethod
    def forward(self, input_1: Any, input_2: Any) -> float:

        raise NotImplementedError
