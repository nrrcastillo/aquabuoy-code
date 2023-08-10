from _collections_abc import Callable, Iterable
from abc import ABC, abstractmethod

import numpy as np
import scipy as sp

"""
TODO:
- Should I make a class for the wave data? Or should I just use a dictionary?
- Should I make a separate class for regular and irregular waves?
  - Likely don't need to do this for linear and nonlinear waves, can use @property.setter with boolean
- [ ] add docstrings
"""


class HD:
    """
    NOTE: Planning to remove this class once NEMOH code is integrated
    """

    def __init__(
        self, T: np.ndarray, F: np.ndarray, a: np.ndarray, b: np.ndarray
    ) -> None:
        self.T = T
        self.F = F
        self.a = a
        self.b = b


"""
Likely do need to inherit from HD class
    - PM spectrum needs interpolated data from HD
"""


class WaveParameters:
    def __init__(
        self,
        regular: bool,
        linear: bool,
        avg_wave_height: float,  #
        avg_wave_period: float,
        fluid_fric_coef: float,
    ) -> None:
        """
        Args:
            regular (bool): True if regular waves, False if irregular waves
            linear (bool): True if linear waves, False if nonlinear waves
            avg_wave_height (float): average wave height
            avg_wave_period (float): average wave period
            fluid_fric_coef (float): fluid friction coefficient. used to compute Fff parameter in DE
        """
        self.regular = regular
        self.linear = linear
        self.avg_wave_height = avg_wave_height  # if regular, this is the wave height. If irregular, this is the significant wave height
        self.avg_wave_period = avg_wave_period
        self.fluid_fric_coef = fluid_fric_coef  # Cd

        ### Non-initalised parameters ###
        self.peak_wave_period = self.avg_wave_period * 1.4  # Tp
        self.peak_wave_freq = 1 / self.peak_wave_period
        ### PM spectrum parameters ###
        self.A = (5 / 16) * (self.avg_wave_height**2) * (self.peak_wave_freq**4)
        self.B = (5 / 4) * (self.peak_wave_freq**2)
        self.fi = (
            0.5 * np.random.rand(40) + np.arange(40)
        ) * 0.01  # random wave phase for each component. 0.01 is df
        self.Ti = 1 / self.fi
        self.Sf = (self.A / (self.fi**5)) * np.exp(
            -self.B / (self.fi**4)
        )  # PM spectrum
        self.ai = np.sqrt(
            2 * self.Sf * self.fi * 0.01
        )  # wave amplitude of each component
        # self.Fwt = np.sum(self.ai * np.sin(2 * np.pi * self.fi * t))


"""
Note: 
- For WAMIT prototype, dry_mass was predefined rather than calculated in-program
- Same is true for other buoy parameters. 
"""

"""
For Reference: (Notation for unknown variables)
t0 = 0       Initial time of simulation
y(1) = z1    Displacement of float mass center of gravity
y(2) = z2    Displacement of piston center of gravity
y(3) = v1    velocity of float mass center of gravity
y(4) = v2    velocity of piston center of gravity
Initial positions and initial velocities: y0 = [z1(t0) z2(t0) v1(t0) v2(t0)] = [0 0 0 0] --> zero dirichlet initial conditions
"""
