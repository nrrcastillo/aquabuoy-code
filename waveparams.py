from _collections_abc import Callable, Iterable
from abc import ABC, abstractmethod
from testdata import (
    T_data as wave_period,
    F_data as wave_ef,
    a_data as added_mass,
    b_data as hydro_damping,
)

import numpy as np
import scipy as sp

"""
TODO:
- Should I make a class for the wave data? Or should I just use a dictionary?
- Should I make a separate class for regular and irregular waves?
  - Likely don't need to do this for linear and nonlinear waves, can use @property.setter with boolean
- [ ] add docstrings
"""


class WaveData:
    def __init__(
        self, T: np.ndarray, F: np.ndarray, a: np.ndarray, b: np.ndarray
    ) -> None:
        self.T = T
        self.F = F
        self.a = a
        self.b = b


class DEParams(WaveData):
    def __init__(
        self,
        regular: bool,
        linear: bool,
        avg_wave_height: float,
        avg_wave_period: float,
    ) -> None:
        """
        Args:
            regular (bool): True if regular waves, False if irregular waves
            linear (bool): True if linear waves, False if nonlinear waves
            avg_wave_height (float): average wave height
            avg_wave_period (float): average wave period
        """
        self.regular = regular
        self.linear = linear
        self.avg_wave_height = avg_wave_height  # if regular, this is the wave height. If irregular, this is the significant wave height
        self.avg_wave_period = avg_wave_period
        self.peak_wave_period = self.avg_wave_period * 1.4
        self.peak_wave_freq = 1 / self.peak_wave_period
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
        self.T_data = None
        self.F_data = None
        self.a_data = None
        self.b_data = None


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
