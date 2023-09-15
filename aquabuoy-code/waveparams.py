from _collections_abc import Callable, Iterable
from abc import ABC, abstractmethod

import numpy as np
import scipy as sp


class WaveParameters:
    def __init__(
        self,
        regular: bool,
        linear: bool,
        avg_wave_height: float,
        avg_wave_period: float,
        drag_coef: float,
        damping_coef: float,
    ) -> None:
        """
        Summary:
            Initialises wave parameters and computes values for irregular waves. PM spectrum is computed here, to be used in irregular exciting force.
            
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
        self.drag_coef = drag_coef  # 
        self.damping_coef = damping_coef  # Ffric for nonlinear

        ### Non-initalised parameters ###
        self.peak_wave_period = self.avg_wave_period * 1.4  # Tp
        self.peak_wave_freq = 1 / self.peak_wave_period

        ### PM spectrum parameters ###
        self.A = (5 / 16) * (self.avg_wave_height**2) * (self.peak_wave_freq**4)
        self.B = (5 / 4) * (self.peak_wave_freq**4)
        self.fi = (
            0.5 * np.random.rand(40) + np.arange(40)
        ) * 0.01  # random wave phase for each component. 0.01 is df
        self.random_phase = np.random.rand(40) * 2 * np.pi
        self.Ti = 1 / self.fi
        self.Sf = (self.A * (self.Ti**5)) * np.exp(
            -self.B * (self.Ti**4)
        )  # PM spectrum
        self.ai = np.sqrt(2 * self.Sf * 0.01)  # wave amplitude of each component

    def _show_params(self):
        print("Wave Parameters:")
        print(f"Regular Wave: {self.regular}")
        print(f"Linear Wave: {self.linear}")
        print(f"Average Wave Height: {self.avg_wave_height}")
        print(f"Average Wave Period: {self.avg_wave_period}")
        print(f"Fluid Friction Coefficient: {self.fluid_fric_coef}")
        print(f"Damping Coefficient: {self.damping_coef}")
        print(f"Peak Wave Period: {self.peak_wave_period}")
        print(f"Peak Wave Frequency: {self.peak_wave_freq}")
        print(f"A: {self.A}")
        print(f"B: {self.B}")
        print(f"fi: {self.fi}")
        print(f"Ti: {self.Ti}")
        print(f"Sf: {self.Sf}")
        print(f"ai: {self.ai}")
