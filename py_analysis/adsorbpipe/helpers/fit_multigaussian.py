import pandas as pd
import numpy as np
from typing import Sequence
from scipy.optimize import curve_fit
import math
from py_analysis.config.custom_types import FitResult

def _multi_gaussian(x: np.ndarray, *pars: float) -> np.ndarray:
    """Sum of N Gaussians; *pars* = [A0, mu0, sigma0,  A1, mu1, sigma1, ...]."""
    y = np.zeros_like(x, dtype=float)
    n = len(pars) // 3
    for i in range(n):
        A, mu, sigma = pars[3*i : 3*i+3]
        y += A * np.exp(-(x - mu)**2 / (2.0 * sigma**2))
    return y


def _initial_guess_bounds(centres: Sequence[int],
                            y: np.ndarray,
                            bp: np.ndarray,
                            *,
                            sigma0: float = 1.0,
                            sigma_min: float = 1.0,
                            sigma_max: float = 5.0) -> tuple[list[float], list[float], list[float]]:
    """Return p0, lower, upper for *curve_fit*."""
    p0: list[float] = []
    lo: list[float] = []
    hi: list[float] = []
    for mu in centres:
        A0 = float(y[bp == mu][0])  # use raw peak height as start
        p0.extend([A0, mu, sigma0])
        lo.extend([0.0, mu - 1.0, sigma_min])
        hi.extend([math.inf, mu + 1.0, sigma_max])
    return p0, lo, hi



def fit_multi_gaussian(dwelltime: np.ndarray | pd.Series,
                        bp: np.ndarray,
                        peakcentres: Sequence[int],
                        *,
                        sigma0: float = 1.0,
                        sigma_min: float = 1.0,
                        sigma_max: float = 5.0) -> FitResult:
    """Fit len(peak centres) Gaussians simultaneously; returns a 'FitResult' """
    p0, lo, hi = _initial_guess_bounds(peakcentres, np.asarray(dwelltime), 
                                       bp, 
                                       sigma0=sigma0, sigma_min=sigma_min, sigma_max=sigma_max)
    pars, pcov = curve_fit(_multi_gaussian,
                            xdata=bp,
                            ydata=dwelltime,
                            p0=p0,
                            bounds=(lo, hi))
    
    return FitResult(np.asarray(peakcentres), pars, pcov)
