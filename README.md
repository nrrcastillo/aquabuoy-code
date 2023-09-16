# aquabuoy-code

Python port of original MATLAB code for numerical simulations of the AquaBuOY wave energy converter.

## Description

This repository contains the Python port of the original MATLAB code for numerical simulations of the AquaBuOY wave energy converter. Models regular and irregular waves and linear and nonlinear damping. The model is described in [Mathematical and numerical modeling of the AquaBuOY wave energy converter.](https://www.researchgate.net/publication/48776715_Mathematical_and_numerical_modeling_of_the_AquaBuOY_wave_energy_converter)

## Getting Started

### Requirements

Download `fpg.py`, `rhsfuncs.py`, `waveparams.py`, `rhs_utils.py`, `testdata.py` from the `Aquabuoy` folder in this repository. If one wishes to use their own hydrodynamic data, then `testdata.py` is not needed.

### Executing program

* Import the module in Python script or jupyter notebook.

``` python

import fpg
import waveparams as wp
import rhsfuncs as rhs
import rhs_utils as _rhs
from testdata import *

```

* Input float specifications for device.
By default, if rho not specified, rho = 1000. If g not specified, g = 9.81. If spring_stiffness_ratio not specificed, spring_stiffness_ratio = 0.5.

``` python

# Device specifications
test_device = fpg.FloatParameters(
    float_diameter=7.0,
    float_draft=4.5,
    tube_length=30.0,
    tube_inner_diameter=4.7,
    tube_outer_diameter=5.0,
    rho=1025.0,
    g=9.81,
    spring_stiffness_ratio=1.0,
)
    
```

* Input corresponding hydrodynamic data for device. If using own hydrodynamic data, input as `np.array`.

``` python

test_data = fpg.HD(
    T=T_data,
    F=F_data,
    a=a_data,
    b=b_data
)

```

* Input wave parameters.
Set boolean `regular = True` if modeling with regular waves. This sets the `Fw(t)` function to use the regular wave function. If using irregular waves, set `regular = False` and `Fw(t)` will point to the irregular wave exciting force function. `avg_wave_height` and `avg_wave_period` will be treated as the significant wave height and significant wave period for irregular waves, respectively.

Set boolean `linear = True` if modeling with linear damping. This indicates whether to use the linear or nonlinear coulomb damping. Set `linear = False` to use nonlinear damping.

Regular waves with linear damping example:

``` python

state_test_rl = wp.WaveParameters(
        regular=True,
        linear=True,
        avg_wave_height=2.5,
        avg_wave_period=7,
        drag_coef=3.0, # used in Fff(t)
        damping_coef=250000,
    )
```

Irregular waves with nonlinear damping example:

``` python
# treats avg_wave_height and avg_wave_period as significant wave height and period
state_test_in = wp.WaveParameters(
    regular=False,
    linear=False,
    avg_wave_height=2.0,
    avg_wave_period=7.0,
    drag_coef=3.0, # used in Fff(t)
    damping_coef=250000,
)
```

* Input buoy specifications.
By default, `rho=1000.0`, `g=9.81`, and `spring_stiffness_ratio=0.5` if not specificed. `spring_stiffness_ratio` should be between 0 and 1. Currently, `dry_mass`, `added_mass`, and `tube_water_mass` are hard-coded to the values used in the original MATLAB code. If you wish to change these values, modify the code in class `FloatParameters` in `fpg.py`. Formulas for calculating these values are also provided as block comments in the class.

Example using full-scale device specifications:

``` python
buoy_test = fpg.FloatParameters(
    float_diameter=7,
    float_draft=4.5,
    tube_length=30,
    tube_inner_diameter=4.7,
    tube_outer_diameter=5.0,
    rho=1025.0,
    g=9.81,
    spring_stiffness_ratio=1.0,
)
```

* Running simulations.

Examples: (using full-scale device specifications)

Note that `end_time=200` for regular waves with linear damping and `end_time=1200` for irregular waves with nonlinear damping.

Solver calls `scipy.integrate.solve_ivp` using `LSODA` method by default. To use other methods such as `RK45` or `BDF`, see below.

``` python

# regular waves with linear damping.
test_rhs_rl = rhs.RHS(
    buoy=buoy_test,
    state=state_test_rl,
    data=test_data,
    end_time=200
)

test_results = test_rhs_rl.solve() # uses LSODA method by default

# irregular waves with nonlinear damping.
test_rhs_in = rhs.RHS(
    buoy=buoy_test,
    state=state_test_in,
    data=test_data,
    end_time=1200
)

test_results = test_rhs_in.solve(use_method=`BDF`) # if using BDF method

```

Note that other ODE solver functions can be used as well. The system of equations can be accessed through `test_rhs_rl.fun` or `test_rhs_in.fun`.

* Plotting results.

To get a summary of the simulation parameters and plots of the results, use the following functions:

``` python
_rhs.display_parameters(test_rhs_nl) # displays parameters used in simulation

_rhs.show_results(test_rhs_nl, test_results) # outputs data summary and plots

```

## Developer Notes

Some potential optimisations to be considered can be found in the block comments of certain functions, such as `hw_sum` in `rhsfuncs.py`.

### Current Issues and Limitations

ODE solver occasionally fails to converge with nonlinear damping. Using same specifications, it may need to be run again to converge. This is currently being investigated.

### Contact and Support

For any questions or issues, please email me at [ncastillo@seattleu.edu](ncastillo@seattleu.edu)
