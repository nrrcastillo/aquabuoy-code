import os
from dataclasses import dataclass, field
from typing import List, Tuple, Dict, Any, Union
from collections import namedtuple

import numpy as np

import cProfile


class HD:
    def __init__(self, T: np.array, F: np.array, a: np.array, b: np.array) -> None:
        self.T = T
        self.F = F
        self.a = a
        self.b = b


@dataclass(frozen=False)
class FloatParameters:
    """
    Getting parameters for aquabuoy and calculating deducible values, operates similar to __dict__

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
    pipe_inner_radius: float
    pipe_outer_radius: float
    tube_length: float
    tube_inner_radius: float
    tube_outer_radius: float
    rho: float = 1000.0  # fluid density
    g: float = 9.81  # g
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
        """if self.pipe_inner_radius >= self.pipe_outer_radius:
            raise FloatParameterError(
                "Pipe inner radius must be less than pipe outer radius"
            )
        if self.tube_inner_radius >= self.tube_outer_radius:
            raise FloatParameterError(
                "Tube inner radius must be less than tube outer radius"
            )"""
        self.water_plane_area = (np.pi / 4) * self.float_diameter**2
        self.dry_mass = self.float_draft * self.rho
        self.extra_mass = (
            (self.tube_outer_radius**2 - self.tube_inner_radius**2)
            * (np.pi / 4)
            * self.tube_length
            * self.rho
        )
        self.total_mass = self.dry_mass + self.extra_mass
        self.tube_water_mass = (
            self.tube_inner_radius**2 * (np.pi / 4) * self.tube_length * self.rho
        )
        self.hydro_stiffness = self.water_plane_area * self.rho * self.g
        self.spring_stiffness = 0.5 * self.hydro_stiffness
        object.__setattr__(self, "_frozen", True)

    def __setattr__(self, name, value):
        if getattr(self, "_frozen", False):
            raise TypeError(f"Cannot assign to field '{name}'")
        super().__setattr__(name, value)

    def __repr__(self):
        return f"FloatParameters(float_diameter={self.float_diameter}, float_draft={self.float_draft}, pipe_inner_radius={self.pipe_inner_radius}, pipe_outer_radius={self.pipe_outer_radius}, tube_length={self.tube_length}, tube_inner_radius={self.tube_inner_radius}, tube_outer_radius={self.tube_outer_radius}, rho={self.rho}, g={self.g})"


class FloatParameterError(Exception):
    pass
