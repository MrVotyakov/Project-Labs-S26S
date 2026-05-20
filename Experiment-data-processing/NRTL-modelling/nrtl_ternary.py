import numpy as np
import pandas as pd
from scipy.optimize import minimize

from nrtl import nrtl, bubbleP_nrtl


def build_tau(params, T):
    tb = np.zeros((3, 3))
    tT = np.zeros((3, 3))
    pairs = [(0,1),(0,2),(1,0),(1,2),(2,0),(2,1)]
    for k, (i, j) in enumerate(pairs):
        tb[i, j] = params[k]
        tT[i, j] = params[k + 6]
    return tb + tT / T


def _objective(params, alpha, data):
    total = 0.0
    for T, x, gamma_exp in data:
        tau = build_tau(params, T)
        try:
            gamma_calc = nrtl(alpha, tau, T, x)
            total += np.sum((np.log(gamma_calc) - np.log(gamma_exp)) ** 2)
        except Exception:
            total += 1e6
    return total


def fit_nrtl_ternary(df, alpha=None, x0=None):
    if alpha is None:
        alpha = np.array([[0,   0.3, 0.3],
                          [0.3, 0,   0.3],
                          [0.3, 0.3, 0  ]])

    data = [
        (row['T'],
         np.array([row['x1'], row['x2'], row['x3']]),
         np.array([row['gamma1'], row['gamma2'], row['gamma3']]))
        for _, row in df.iterrows()
    ]

    if x0 is None:
        x0 = np.array([3.0, 4.0, -2.0, 2.0, -3.0, -2.0,
                       -700., -900., 500., -600., 800., 600.])

    result = minimize(
        _objective, x0,
        args=(alpha, data),
        method='Nelder-Mead',
        options={'maxiter': 50000, 'xatol': 1e-8, 'fatol': 1e-8}
    )

    return result.x, alpha


class NRTLTernary:
    def __init__(self, params, alpha):
        self.params = params
        self.alpha  = alpha

    @classmethod
    def from_csv(cls, path, alpha=None, x0=None):
        df = pd.read_csv(path, skipinitialspace=True)
        df.columns = df.columns.str.strip()
        params, alpha = fit_nrtl_ternary(df, alpha=alpha, x0=x0)
        return cls(params, alpha)

    def get_gamma(self, T, x1, x2, x3):
        x = np.array([x1, x2, x3], dtype=float)
        x /= x.sum()
        tau = build_tau(self.params, T)
        return nrtl(self.alpha, tau, T, x)