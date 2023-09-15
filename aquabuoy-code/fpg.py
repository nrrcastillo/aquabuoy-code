import os
from dataclasses import dataclass, field
from typing import List, Tuple, Dict, Any, Union
from collections import namedtuple

import numpy as np
import scipy as sp


class HD:
    """
    Summary:
        Class to store hydrodynamic data. Interpolates data using cubic spline interpolation.
    
    Args:
        T (np.array): period data
        F (np.array): force data
        a (np.array): added mass data
        b (np.array): fluid damping data -- NOTE: Check if this is correct
    """
    def __init__(self, T: np.array, F: np.array, a: np.array, b: np.array) -> None:
        self.T = T
        self.F = F
        self.a = a
        self.b = b
        self.F_interp = sp.interpolate.CubicSpline(self.T, self.F, bc_type="natural")
        self.a_interp = sp.interpolate.CubicSpline(self.T, self.a, bc_type="natural")
        self.b_interp = sp.interpolate.CubicSpline(self.T, self.b, bc_type="natural")


@dataclass(frozen=False)
class FloatParameters:
    """
    Summary:
        Stores float specifications and physical parameters. 

    Raises:
        FloatParameterError: if any values are negative.
        FloatParameterError: if outer parameters are smaller than inner parameters.
    Args: (Only will write about more confusing ones first)
        rho (float): denoted as rho in paper. default = 1000 \n
        g (float): denoted as g in paper. default = 9.81 \n
    Returns:
        all other arguments specified above \n
        water_plane_area (float): computed using float_diameter \n
        dry_mass (float): computed using float_draft and rho \n
        extra_mass (float): get back to this later
    """

    float_diameter: float
    float_draft: float
    tube_length: float
    tube_inner_diameter: float  # treat like diameter, same formulas so could just change variable names later
    tube_outer_diameter: float
    rho: float = 1000.0  # fluid density
    g: float = 9.81  # g
    spring_stiffness_ratio: float = 0.5  # c/c_in value with default c_in = 0.5*S
    water_plane_area: float = field(init=False)
    dry_mass: float = field(init=False)
    extra_mass: float = field(init=False)  # a, or
    total_mass: float = field(init=False)  ## M1
    tube_water_mass: float = field(init=False)  # M2
    hydro_stiffness: float = field(init=False)  # S in paper
    spring_stiffness: float = field(
        init=False
    )  # c or c_in in matlab (PTO spring stiffness)

    def __post_init__(self):
        self._inspect_init()
        self.water_plane_area = (np.pi / 4) * (self.float_diameter**2)
        """self.dry_mass = self.float_draft * self.rho * self.water_plane_area
        self.extra_mass = (
            (self.tube_outer_diameter**2 - self.tube_inner_diameter**2)
            * (np.pi / 4)
            * self.tube_length
            * self.rho
        )"""
        self.dry_mass = 237000  # M1
        self.added_mass = 149000
        self.total_mass = self.dry_mass + self.added_mass
        """self.tube_water_mass = (
            (self.tube_inner_diameter**2) * (np.pi / 4) * self.tube_length * self.rho
        )"""
        self.tube_water_mass = 363000
        self.hydro_stiffness = self.water_plane_area * self.rho * self.g
        self.spring_stiffness = self.spring_stiffness_ratio * self.hydro_stiffness
        # object.__setattr__(self, "_frozen", True)

    def _inspect_init(self):
        if self.float_diameter <= 0:
            raise FloatParameterError("Float diameter must be positive")
        if self.float_draft <= 0:
            raise FloatParameterError("Float draft must be positive")
        if self.tube_length <= 0:
            raise FloatParameterError("Tube length must be positive")
        if self.tube_inner_diameter <= 0:
            raise FloatParameterError("Tube inner diameter must be positive")
        if self.tube_outer_diameter <= 0:
            raise FloatParameterError("Tube outer diameter must be positive")
        if self.rho <= 0:
            raise FloatParameterError("Fluid density must be positive")
        if self.g <= 0:
            raise FloatParameterError("Gravity must be positive")
        if self.spring_stiffness_ratio < 0:
            raise FloatParameterError("Spring stiffness ratio must be non-negative")

        if self.tube_inner_diameter >= self.tube_outer_diameter:
            raise FloatParameterError(
                "Tube inner diameter must be smaller than tube outer diameter"
            )
        if self.tube_inner_diameter >= self.float_diameter:
            raise FloatParameterError(
                "Tube inner diameter must be smaller than float diameter"
            )

    def _show_params(self):
        print("Parameters:")
        for key, value in self.__dict__.items():
            print(key, ":", value)


class FloatParameterError(Exception):
    pass
