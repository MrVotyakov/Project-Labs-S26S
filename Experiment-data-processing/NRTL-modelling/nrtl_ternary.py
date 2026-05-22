import numpy as np
import pandas as pd
from scipy.optimize import minimize, differential_evolution
from functools import partial
from nrtl import nrtl, bubbleP_nrtl


def build_tau(params, T):
    pairs = [(0,1),(0,2),(1,0),(1,2),(2,0),(2,1)]
    tau = np.zeros((3, 3))
    for k, (i, j) in enumerate(pairs):
        tau[i, j] = params[k] + params[k + 6] / T
    return tau


def _objective(params, data):
    tau_base  = params[:6]
    tau_T     = params[6:12]
    alpha_val = params[12]

    alpha = np.array([[0,         alpha_val, alpha_val],
                      [alpha_val, 0,         alpha_val],
                      [alpha_val, alpha_val, 0        ]])

    pairs = [(0,1),(0,2),(1,0),(1,2),(2,0),(2,1)]
    total = 0.0
    for T, x, gamma_exp in data:
        tau = np.zeros((3, 3))
        for k, (i, j) in enumerate(pairs):
            tau[i, j] = tau_base[k] + tau_T[k] / T
        try:
            gamma_calc = nrtl(alpha, tau, T, x)
            gamma_calc = np.clip(gamma_calc, 1e-10, None)
            total += np.sum((np.log(gamma_calc) - np.log(gamma_exp)) ** 2)
        except Exception:
            total += 1e6
    return total


def fit_nrtl_ternary(df, x0=None):
    data = [
        (row['T'],
         np.array([row['x1'], row['x2'], row['x3']]),
         np.array([row['gamma1'], row['gamma2'], row['gamma3']]))
        for _, row in df.iterrows()
    ]

    bounds = [(-15, 15)] * 6 + [(-5000, 5000)] * 6 + [(0.1, 0.6)]
    obj = partial(_objective, data=data)

    best_result = None
    for seed in [42, 123, 777, 1337, 9999]:
        r = differential_evolution(
            obj, bounds,
            maxiter=2000,
            popsize=10,
            mutation=(0.5, 1.9), recombination=0.9,
            seed=seed, workers=-1, tol=1e-8, disp=False
        )
        r2 = minimize(obj, r.x, method='L-BFGS-B', bounds=bounds,
                      options={'maxiter': 10000, 'ftol': 1e-14})
        best = r2 if r2.fun < r.fun else r
        print(f"seed={seed:6d}  MSE={best.fun:.4e}")
        if best_result is None or best.fun < best_result.fun:
            best_result = best

    print(f"\nЛучший MSE: {best_result.fun:.4e}")

    best_x     = best_result.x
    tau_params = best_x[:12]
    alpha_val  = best_x[12]
    alpha = np.array([[0,         alpha_val, alpha_val],
                      [alpha_val, 0,         alpha_val],
                      [alpha_val, alpha_val, 0        ]])

    return tau_params, alpha


class NRTLTernary:
    def __init__(self, params, alpha):
        self.params = params
        self.alpha  = alpha

    @classmethod
    def from_csv(cls, path, x0=None):
        df = pd.read_csv(path, skipinitialspace=True)
        df.columns = df.columns.str.strip()
        df['T'] = df['T'] + 273.15
        df = df[df['gamma3'] > 0].reset_index(drop=True)
        params, alpha = fit_nrtl_ternary(df, x0=x0)
        return cls(params, alpha)

    @classmethod
    def from_df(cls, df):
        params, alpha = fit_nrtl_ternary(df)
        return cls(params, alpha)

    def get_gamma(self, T, x1, x2, x3):
        x = np.array([x1, x2, x3], dtype=float)
        x /= x.sum()
        tau = build_tau(self.params, T)
        return nrtl(self.alpha, tau, T, x)