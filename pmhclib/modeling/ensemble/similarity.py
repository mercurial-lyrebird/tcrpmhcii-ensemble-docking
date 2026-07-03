"""
    Adapted from:

        A. Aina, S. C. C. Hsueh, and S. S. Plotkin, “PROTHON: A Local
        Order Parameter-Based Method for Efficient Comparison of Protein
        Ensembles,” J. Chem. Inf. Model., vol. 63, no. 11, pp. 3453-3461, Jun.
        2023, doi: 10.1021/acs.jcim.3c00145.

        https://github.com/PlotkinLab/Prothon
"""

from itertools import product
from typing import List

import numpy as np
from scipy.spatial.distance import jensenshannon
from sklearn.neighbors import KernelDensity


def get_estimators(X: np.ndarray) -> List[KernelDensity]:
    """
        X \in R^{m * n * d}

        returns n independent estimators of d-ary residue-wise feature
        distributions fitted over m samples
    """

    n = X.shape[1]

    models = list()
    for i in range(n):
        model = KernelDensity(bandwidth="silverman", kernel="gaussian")
        model = model.fit(X[:, i, :])
        models.append(model)

    return models


def sample_distributions(
    X: np.ndarray, models: List[KernelDensity]
) -> np.ndarray:
    """
        X \in R^{m * d}

        p \in R^{m * n}
    """

    [m, d] = X.shape
    n = len(models)
    p = np.zeros((m, n))

    for i in range(n):
        p[:, i] = models[i].score_samples(X)

    return p


def estimate_and_sample(
    X: np.ndarray, min: float, max: float, v: int = 100
) -> np.ndarray:

    [m, n] = X.shape[:2]
    if len(X.shape) == 2:
        X = np.reshape(X, (m, n, 1))
    d = X.shape[2]

    models = get_estimators(X)
    vals = np.linspace(min, max, v)
    samples = np.array(list(product(vals, repeat=d)))
    p = sample_distributions(samples, models)

    return p


def local_jsd(
    X: np.ndarray, Y: np.ndarray, min: float = None, max: float = None,
    v: int = 100
) -> np.ndarray:

    if min is None:
        min = np.min([X.min(), Y.min()])
    if max is None:
        max = np.max([X.max(), Y.max()])

    p_x = estimate_and_sample(X, min, max, v=v)
    p_y = estimate_and_sample(Y, min, max, v=v)

    dists = jensenshannon(p_x, p_y, base=2)
    dists[np.isinf(dists)] = 0
    dists[np.isnan(dists)] = 0

    return dists
