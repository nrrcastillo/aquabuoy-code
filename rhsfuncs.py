from testdata import T, F, a, b
from typing import Callable

import numpy as np
import scipy as sp

"""
For tspan(-250:dtau:endtime) -> counterpart in sp.integrate is t_eval
- starts 
- t_eval parameter specifies times to store computed solution

Note:
- to pass these parameters into sp.integrate, need to pass as tuple
"""


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


def hw_sum(
    new_time: np.array,
    new_velocity: np.array,
    b_i: np.array,
    omega_i: np.array,
    t: float,
) -> float:
    """
    NOTE: Could likely  optimise this to not have to go through if statement every time
    - t array is current time
    """
    dt = 0.5
    n = np.arange(1, 14)
    t_n = t - n * dt
    if len(new_time) == 1:
        vel0 = new_velocity[0]  # for 0th term in sum
        vel = new_velocity[0]
    else:
        vel0 = np.interp(
            t, new_time, new_velocity
        )  # interpolated velocity for 0th term in sum
        vel = np.interp(t_n, new_time, new_velocity)
    res0 = 0.5 * hw(new_time, b_i, omega_i, t) * vel0 * dt
    res = np.sum(
        hw(b_i, omega_i, t_n) * vel * dt
    )  # uses vectorised version of hw (replaces for loop)


def hw(b_i: np.array, omega_i: np.array, t: float) -> float:
    """
    NOTE: input t depends on calculation outside of this function.

    Summary:
        Retardation function for irregular waves
        Cannot be calculated if time array is length == 1
    Args:
        b_i (Callable[[float], float]): Interpolated function of given wave data
        omega_i (Callable[[float], float]): Interpolated function of given wave data
        t (float): time
    """

    temp = (
        (1 / 2)
        * (b_i[1:] + b_i[:-1])
        * np.cos(0.5 * (omega[1:] + omega_i[:-1]) * t)
        * (omega_i[:-1] - omega_i[1:])
    )
    return (2 / np.pi) * np.sum(temp)


def func(
    t: float,
    y: np.array,
    b_i: np.array,
    omega_i: np.array,
    spring_stiffness: float, # Same as c_in, given by 0.5*S
    b2: float, # friction between float and piston, same as Ffric in paper (using nonlinear)
    hydrodynamic_stiffness: float,
    Fff: float, # might be better written directly into the function
    M1: float,
    M2: float,
    am: float, # TODO: Should I use the integral of the wave data instead?
) -> np.array:  # for irregular and nonlinear waves
    """
    Parameters:
        t (float): time
        y (np.array): current state vector [v1, v2, z1, z2]
        b_i (np.array): hydrodynamic damping data
        omega_i (np.array): wave frequency data
        c_in (float): spring stiffness
    """
    pass
