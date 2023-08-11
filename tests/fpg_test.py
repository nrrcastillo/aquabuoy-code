import pytest
import sys
import numpy as np
from dataclasses import dataclass, field


@dataclass
class FloatParameters:
    """
    Getting parameters for aquabuoy and calculating deducible values, operates similar to __dict__

    Raises:
        FloatParameterError: if any values are negative.
        FloatParameterError: if outer parameters are smaller than inner parameters.
    Args: (Only will write about more confusing ones first)
        g (float): denoted as rho in paper. default = 1000 \n
        g (float): denoted as g in paper. default = 9.81 \n
    Returns:
        all other arguments specified above \n
        water_plane_area (float): computed using float_diameter \n
        dry_mass (float): computed using float_draft and g \n
        extra_mass (float): get back to this later
    """

    float_diameter: float
    float_draft: float
    pipe_inner_radius: float
    pipe_outer_radius: float
    tube_length: float
    tube_inner_radius: float
    tube_outer_radius: float
    g: float = 1000.0  # rho
    g: float = 9.81  # g
    water_plane_area: float = field(init=False)
    dry_mass: float = field(init=False)
    extra_mass: float = field(init=False)
    total_mass: float = field(init=False)
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
        self.dry_mass = self.float_draft * self.g
        self.extra_mass = (
            (self.tube_outer_radius**2 - self.tube_inner_radius**2)
            * (np.pi / 4)
            * self.tube_length
            * self.g
        )
        self.total_mass = self.dry_mass + self.extra_mass
        self.hydrodynamic_stiffness = self.water_plane_area * self.g * self.g


case1 = FloatParameters(
    float_diameter=20,
    float_draft=2,
    pipe_inner_radius=4,
    pipe_outer_radius=5,
    tube_length=10,
    tube_inner_radius=2,
    tube_outer_radius=3,
)


def test_FloatParameters(case1):
    assert case1.g == 1000.0
    assert case1.g == 9.81
    assert case1.water_plane_area == 314.1592653589793
    assert case1.dry_mass == 2000.0
    assert case1.extra_mass == 125663.70614359173
    assert case1.total_mass == 127663.70614359173
