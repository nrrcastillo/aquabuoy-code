import os
from dataclasses import dataclass, field
from typing import List, Tuple, Dict, Any, Union
from collections import defaultdict

import numpy as np


@dataclass
class FloatParameters:
    """
    Getting parameters for aquabuoy and calculating deducible values, operates similar to __dict__

    Raises:
        FloatParameterError: if any values are negative.
        FloatParameterError: if outer parameters are smaller than inner parameters.
    Args: (Only will write about more confusing ones first)
        fluid_density (float): denoted as rho in paper. default = 1000 \n
        acceleration_due_to_gravity (float): denoted as g in paper. default = 9.81 \n
    Returns:
        all other arguments specified above \n
        water_plane_area (float): computed using float_diameter \n
        dry_mass (float): computed using float_draft and fluid_density \n
        extra_mass (float): get back to this later
    """

    float_diameter: float
    float_draft: float
    pipe_inner_radius: float
    pipe_outer_radius: float
    tube_length: float
    tube_inner_radius: float
    tube_outer_radius: float
    fluid_density: float = 1000.0  # rho
    acceleration_due_to_gravity: float = 9.81  # g
    water_plane_area: float = field(init=False)
    dry_mass: float = field(init=False)
    extra_mass: float = field(init=False)  # a, or
    total_mass: float = field(init=False)  ## M1
    hydro_stiffness: float = field(init=False)

    def __post_init__(self):
        if self.pipe_inner_radius >= self.pipe_outer_radius:
            raise FloatParameterError(
                "Pipe inner radius must be less than pipe outer radius"
            )
        if self.tube_inner_radius >= self.tube_outer_radius:
            raise FloatParameterError(
                "Tube inner radius must be less than tube outer radius"
            )
        self.water_plane_area = (np.pi / 4) * self.float_diameter**2
        self.dry_mass = self.float_draft * self.fluid_density
        self.extra_mass = (
            (self.tube_outer_radius**2 - self.tube_inner_radius**2)
            * (np.pi / 4)
            * self.tube_length
            * self.fluid_density
        )
        self.total_mass = self.dry_mass + self.extra_mass
        self.hydrodynamic_stiffness = (
            self.water_plane_area
            * self.fluid_density
            * self.acceleration_due_to_gravity
        )


class FloatParameterError(Exception):
    pass
