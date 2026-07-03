
import numpy as np
from scipy.spatial.distance import euclidean
from torch import Tensor

from ._metric import Metric


class Euclidean(Metric):

    def __init__(self):

        super(Euclidean, self).__init__()

        return

    def forward(self, input_1: Tensor, input_2: Tensor) -> float:

        return euclidean(input_1, input_2)


class Norm(Metric):

    def __init__(self, ord: int = 2):

        super(Norm, self).__init__()

        self.ord = ord

        return

    def forward(self, input_1: Tensor, input_2: Tensor) -> float:

        return np.linalg.norm(input_1 - input_2, ord=self.ord)
