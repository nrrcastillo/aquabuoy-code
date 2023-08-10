from typing import Callable
from fpg import FloatParameters
from waveparams import WaveParameters, HD

import numpy as np
import scipy as sp
from scipy.interpolate import interp1d

"""
For tspan(-250:dtau:endtime) -> counterpart in sp.integrate is t_eval
- starts 
- t_eval parameter specifies times to store computed solution

Note:
- to pass these parameters into sp.integrate, need to pass as tuple
"""


class RHS:
    def __init__(self, buoy: FloatParameters, state: WaveParameters, data: HD) -> None:
        self.buoy = buoy
        self.state = state
        self.data = data
        self.F_i = interp1d(self.data.T, self.data.F, kind="cubic")  # F_i = F(T)
        self.a_i = interp1d(self.data.T, self.data.a, kind="cubic")  # a_i = a(T)
        self.b_i = interp1d(self.data.T, self.data.b, kind="cubic")  # b_i = b(T)
        self._jspars = np.array(
            [[0, 0, 1, 0], [0, 0, 0, 1], [1, 1, 1, 1], [1, 1, 1, 1]]
        )

    # Fw only needs to be called once
    def f(self, t: float, y: np.array):
        v1p = (1 / self.buoy.total_mass) * (
            self.Fw
            - self.hw_sum()
            - self.coulomb_damping(v1, v2)
            - self.buoy.hydro_stiffness * z1
            - self.Fff
        )
        v2p = (1 / self.buoy.tube_water_mass) * (self.Fpto(v1, v2))

    def hw_sum(
        self,
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
        res0 = 0.5 * self.hw(b_i, omega_i, 0) * vel0 * dt
        res = np.sum(
            self.hw(b_i, omega_i, t_n) * vel * dt
        )  # uses vectorised version of hw (replaces for loop)
        return res + res0

    def coulomb_damping(self, v1: float, v2: float) -> float:
        return self.Ffric * np.tanh(100 * (v1 - v2))

    def hw(self, b_i: np.array, omega_i: np.array, t: float) -> float:
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
            * np.cos(0.5 * (omega_i[1:] + omega_i[:-1]) * t)
            * (omega_i[:-1] - omega_i[1:])
        )
        return (2 / np.pi) * np.sum(temp)

    def Fpto(self, v1: float, v2: float) -> float:
        """
        Coulomb damping force for nonlinear damping

        Args:
            Ffric (float): Coulomb friction force
            v1 (float): velocity of body 1
            v2 (float): velocity of body 2

        Returns:
            float: Coulomb damping force
        """
        return (
            self.state.fluid_fric_coef
            * self.buoy.water_plane_area
            * self.buoy.g
            * np.tanh(100 * (v1 - v2))
        )

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


# NOTE: (1) Can move this variable elsewhere
def func(
    t: float,
    y: np.array,
    b_i: np.array,  # (1)
    omega_i: np.array,  # (2)``
    b2: float,  # friction between float and piston, same as Ffric in paper (using nonlinear)
    Fff: float,  # might be better written directly into the function
    am: float,  # TODO: Should I use the integral of the wave data instead?
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
