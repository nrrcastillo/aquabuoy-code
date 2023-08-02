from testdata import T, F, a, b
from typing import Callable

import numpy as np
import scipy as sp


def exciting_force(
    F: Callable[[float], float], t: np.array, A: float, T: float
) -> float:
    """
    Summary:
        Exciting force for regular waves
    Args:
        F (Callable[[float], float]): Interpolated function of given wave data
        t (np.array): time
        A (float): amplitude, same as H/2
        T (float): wave period
    """
    omega = 2 * np.pi / T
    gamma = 0
    return A * np.cos(omega * t + gamma)


def hw(b: Callable[[float], float], omega: Callable[[float], float], t: float) -> float:
    """
    NOTE: input t depends on calculation outside of this function.

    Summary:
        Retardation function for irregular waves
    Args:
        b (Callable[[float], float]): Interpolated function of given wave data
        omega (Callable[[float], float]): Interpolated function of given wave data
        t (float): time
    """
    return np.trapz(b(t) * omega(t), t)
