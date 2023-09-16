import numpy as np
import scipy as sp

import rhs_utils as _rhs
from fpg import HD, FloatParameters
from waveparams import WaveParameters


class RHS:
    def __init__(
        self, buoy: FloatParameters, state: WaveParameters, data: HD, end_time: float
    ) -> None:
        """
        Summary:
            Generates right hand side of system of equations.

        Args:
            buoy (FloatParameters): instance of FloatParameters class
            state (WaveParameters): instance of WaveParameters class
            data (HD): instance of HD class
            end_time (float): end time for simulation
        """
        self.buoy = buoy
        self.state = state
        self.data = data
        self.end_time = end_time

        self._jspars = np.array(
            [[0, 0, 1, 0], [0, 0, 0, 1], [1, 1, 1, 1], [1, 1, 1, 1]]
        )

        self.T_i = np.linspace(self.data.T[0], self.data.T[-1], 85)
        self.omega_i = 2 * np.pi / self.T_i
        self.b_i = self.data.b_interp(self.T_i)

        self.t_Fa = np.arange(-250, self.end_time + 1)  # given dtau = 1. This is tspan

        if self.state.regular == True:
            self.Fw = self.EF
        else:
            self.Fw = _rhs.set_Fw_nl(self.state, self.data, self.t_Fa)
            self.buoy.added_mass += 46000

        if self.state.linear == True:
            self.coulomb_damping = self.linear_damping
            self.hydro_damping = self.bv
        else:
            self.coulomb_damping = self.nonlinear_damping
            self.hydro_damping = self.hw_sum

        self.t0 = np.round(1 - self.t_Fa[0])
        self._time_data = []
        self._velocity_data = []

    def EF(self, t: float) -> float:
        """
        Summary:
            Calculates exciting force for regular waves.

        Args:
            t (float): current time

        Returns:
            EF(t) (float): exciting force at time t
        """
        amplitude = self.state.avg_wave_height / 2
        cycliq_freq = 2 * np.pi / self.state.avg_wave_period  # omega
        phase = 0
        return (
            amplitude
            * self.F_interp(self.state.avg_wave_period)
            * np.sin(cycliq_freq * t + phase)
        )

    def update_data(self, t: float, y: np.array) -> np.array:
        """
        Summary:
            Updates the time and velocity data arrays and verifies that unique time values are stored. Velocity is used to calculate retardation integral by pchip interpolation. Only uses most recent 200 values.
            Not called if linear==True

        Args:
            t (float): current time
            y (np.array): current velocity

        Returns:
            _time_data (np.array): updated array of time values
            _velocity_data (np.array): updated array of velocity values

        NOTE: Consider allocating memory for these arrays beforehand to avoid copying every time
            - could use np.empty() to allocate memory and then np.put() to update values
            - although need to have an estimate for how many values will be stored
            - If I proceed with this, probably allocate array of size 2^n
        """
        self._time_data.append(t)
        self._velocity_data.append(y[2])
        new_time, I = np.unique(self._time_data, return_index=True)
        new_velocity = [self._velocity_data[i] for i in I]

        self._time_data = new_time.tolist()
        self._velocity_data = new_velocity

        if len(self._time_data) > 200:
            return np.array(self._time_data[-200:]), np.array(
                self._velocity_data[-200:]
            )
        else:
            return np.array(self._time_data), np.array(self._velocity_data)

    # y = [z1, z2, v1, v2]
    def fun(self, t: float, y: np.array):
        """
        Summary:
            Right hand side of the system of equations.

        Args:
            t (float): current time
            y (np.array): current state vector
        Returns:
            dydt (np.array): RHS of system of equations
        """
        v1p = (1 / self.buoy.total_mass) * (
            self.Fw(t)
            - self.hydro_damping(t, y)
            - self.coulomb_damping(y[2], y[3])
            - self.buoy.hydro_stiffness * y[0]
            - self.buoy.spring_stiffness * (y[0] - y[1])
            - self.Fff(y[2])
        )
        v2p = (1 / self.buoy.tube_water_mass) * (
            self._coulomb_damping(y[2], y[3])
            + self.buoy.spring_stiffness * (y[0] - y[1])
        )
        return [y[2], y[3], v1p, v2p]

    def bv(self, t: float, y: np.array) -> float:
        """
        Summary:
            Calculates the hydrodynamic damping in linear damping case.

        Args:
            t (float): current time
            y (np.array): current state vector

        Returns:
            float: hydrodynamic damping value (b*v1 in paper)
        """
        return self.data.b_interp(self.state.avg_wave_period) * y[2]

    def hw_sum(
        self,
        t: float,
        y: np.array,
    ) -> float:
        """
        Summary:
            Summation of retardation integral for irregular waves. Calls hw(t) function to calculate retardation function for each time step.
        Args:
            t (float): current time
            y (np.array): current state vector
        Returns:
            res (float): summation of retardation integral

        TODO: Should try to optimise
            1. if-else statement can be avoided by implementing a flag to determine if len(new_time) > 1.
                - possibly use a class variable to store this flag
                - use different function depending on flag
                - could change with @property decorator
            2. Potentially vectorise this function to avoid for loop
                - initially tried with tn = t - np.arange(1, 14) * 0.5
                - Problem is it casts arrays with uneven length inside hw(t) function
            3. Find to avoid recalculation of array operations inside hw(t)
                - could try to store these values in a class variable.
                   ex. self.dw = self.omega_i[:-1] - self.omega_i[1:]
                - For some reason, this seems to throw an error when I try to do this (ODE solver doesn't like it)
        """
        new_time, new_velocity = self.update_data(t, y)
        dt = 0.5
        n = np.arange(1, 14)
        t_n = t - n * dt
        res = 0
        if len(new_time) == 1:
            vel = new_velocity[0]
            for i in n:
                res += self.hw(i * dt) * vel * dt
            return res + 0.5 * self.hw(0) * vel * dt * y[2]
        else:
            vel = sp.interpolate.PchipInterpolator(new_time, new_velocity)
            for i in n:
                res += self.hw(i * dt) * vel(t - i * dt) * dt
            return res + 0.5 * self.hw(0) * vel(t) * dt * y[2]

    def hw(self, t: float) -> float:
        """
        NOTE: input t depends on function calls from hw_sum

        Summary:
            Retardation function for irregular waves
            Cannot be calculated if time array is length == 1

        Args:
            t (float): current time
        Returns:
            res (float): retardation function at time t
        """

        temp = (
            0.5
            * (self.b_i[1:] + self.b_i[:-1])
            * np.cos(0.5 * (self.omega_i[1:] + self.omega_i[:-1]) * t)
            * (self.omega_i[:-1] - self.omega_i[1:])
        )
        res = (2 / np.pi) * np.sum(temp)
        return res

    def linear_damping(self, v1: float, v2: float) -> float:
        """
        Summary:
            Calculates linear coulomb damping.

        Args:
            v1 (float): velocity of float centre of gravity
            v2 (float): velocity of piston centre of gravity

        Returns:
            float: coulomb damping value
        """
        return self.state.damping_coef * (v1 - v2)

    def nonlinear_damping(self, v1: float, v2: float) -> float:
        """
        Summary:
            Calculates nonlinear coulomb damping.

        Args:
            v1 (float): velocity of float centre of gravity
            v2 (float): velocity of piston centre of gravity

        Returns:
            float: coulomb damping value
        """
        return self.state.damping_coef * np.tanh(100 * (v1 - v2))

    def Fff(self, v1: float) -> float:
        """
        Summary:
            Calculates fluid friction for a given velocity

        Args:
            v1 (float): velocity of float centre of gravity

        Returns:
            float: fluid friction value
        """
        return (
            0.5
            * self.buoy.water_plane_area
            * self.state.drag_coef
            * self.buoy.rho
            * v1
            * np.abs(v1)
        )

    def solve(self, use_method: str = "LSODA"):
        """
        NOTE:
            Do not necessarily need to use this function. Could just call sp.integrate.solve_ivp directly.
        """
        solver = sp.integrate.solve_ivp(
            fun=self.fun,
            t_span=(-250, self.end_time),
            y0=[0, 0, 0, 0],
            method=use_method,
            # jac_sparsity=self._jspars,
            rtol=1e-4,
            atol=1e-6,
            t_eval=self.t_Fa[self.t0 :],
        )
        
        return solver
