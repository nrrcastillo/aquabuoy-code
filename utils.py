import numpy as np
import scipy as sp


def relative_displacement(z1: np.array, z2: np.array) -> np.array:
    """
    Summary:
        Computes relative displacement between float and piston
    Args:
        z1 (np.array): float displacement
        z2 (np.array): piston displacement
    """
    return z1 - z2


def fpto():
    pass
