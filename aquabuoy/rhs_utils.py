import numpy as np
import scipy as sp
import waveparams, fpg
import matplotlib.pyplot as plt


def set_Fw_lin(wave_height: float, wave_period: float):
    amplitude = wave_height / 2


def set_Fw_nl(
    state: waveparams.WaveParameters,
    data: fpg.HD,
    t_Fa: np.array,
):
    Fa = np.zeros(len(t_Fa))
    for i in range(len(t_Fa)):
        Fa[i] = np.sum(
            state.ai
            * np.sin(
                (2 * np.pi * state.fi * i) + state.random_phase
            )  # NOTE: random_phase is a scalar
            * data.F_interp(state.Ti),
        )
    Fw = sp.interpolate.CubicSpline(t_Fa, Fa, bc_type="natural")
    return Fw

def display_parameters(rhs):
    if rhs.state.regular:
        print("Regular Waves")
    else:
        print("Irregular Waves")
    if rhs.state.linear:
        print("Linear Damping", end="\n\n")
    else:
        print("Nonlinear Damping", end="\n\n")

    if rhs.state.regular:
        print(f"Wave Height: {rhs.state.avg_wave_height}")
        print(f"Wave Period: {rhs.state.avg_wave_period}", end="\n\n")
    else:
        print(f"Significant Wave Height: {rhs.state.avg_wave_height}")
        print(f"Average Wave Period: {rhs.state.avg_wave_period}", end="\n\n")

    print(f"Hydro Stiffness (S): {rhs.buoy.hydro_stiffness}")
    print(f"Spring Stiffness (c): {rhs.buoy.spring_stiffness}")
    print(f"Drag Coefficient (Cd): {rhs.state.drag_coef}")
    print(f"Damping Coefficient (b2): {rhs.state.damping_coef}", end="\n\n")

    print(f"rho: {rhs.buoy.rho}")
    print(f"g: {rhs.buoy.g}")
    print(f"Dry Mass (M1): {rhs.buoy.dry_mass}")
    print(f"Added Mass (am): {rhs.buoy.added_mass}")
    print(f"Water Mass in Tube (M2): {rhs.buoy.tube_water_mass}")
    print(f"Tube Length (Lt): {rhs.buoy.tube_length}")
    print(f"Tube Inner Diameter (Dt_in): {rhs.buoy.tube_inner_diameter}")
    print(f"Tube Outer Diameter (Dt_out): {rhs.buoy.tube_outer_diameter}")
    print(f"Float Diameter (D): {rhs.buoy.float_diameter}")
    print(f"Float Draft (hf): {rhs.buoy.float_draft}", end="\n\n")


def show_results(rhs, results):
    ts, ys = results.t, results.y
    if len(ts) == 0:
        print(results.message)
        return
    else:
        print("Generating Results...")
    zr = ys[0] - ys[1]  # relative displacement
    vr = ys[2] - ys[3]  # relative velocity

    if rhs.state.linear:
        force = (rhs.state.damping_coef * vr + rhs.buoy.spring_stiffness * zr) / 1000
    else:
        force = (
            rhs.state.damping_coef * np.tanh(100 * vr) + rhs.buoy.spring_stiffness * zr
        ) / 1000
    avg_force = np.cumsum(np.abs(force)) / np.arange(1, len(force) + 1)
    power = force * vr
    avg_power = np.cumsum(power) / np.arange(1, len(power) + 1)

    print(f"Peak Float Amplitude [m]: {np.max(np.abs(ys[0]))}")
    print(f"Average Float Amplitude [m]: {np.mean(np.abs(ys[0]))}")
    print(
        f"Root Mean Square Float Amplitude [m]: {np.sqrt(np.mean(np.power(ys[0], 2)))}",
        end="\n\n",
    )

    print(f"Peak Relative Amplitude [m]: {np.max(np.abs(zr))}")
    print(f"Average Relative Amplitude [m]: {np.mean(np.abs(zr))}")
    print(
        f"Root Mean Square Relative Amplitude [m]: {np.sqrt(np.mean(np.power(zr, 2)))}",
        end="\n\n",
    )

    print(f"Peak PTO Force [kN]: {np.max(force)}")
    print(f"Average PTO Force [kN]: {np.mean(np.abs(force))}")
    print(
        f"Root Mean Square PTO Force [kN]: {np.sqrt(np.mean(np.power(force, 2)))}",
        end="\n\n",
    )

    print(f"Peak Power [kW]: {np.max(power)}")
    print(f"Average Power [kW]: {np.mean(power)}")
    print(
        f"Root Mean Square Power [kW]: {np.sqrt(np.mean(np.power(power, 2)))}",
        end="\n\n",
    )

    print(f"Peak Velocity of Float [m/s]: {np.max(ys[2])}")
    print(f"Average Velocity of Float [m/s]: {np.mean(np.abs(ys[2]))}")
    print(
        f"Root Mean Square Velocity of Float [m/s]: {np.sqrt(np.mean(np.power(ys[2], 2)))}",
        end="\n\n",
    )

    print(f"Peak Velocity of Piston [m/s]: {np.max(ys[3])}")
    print(f"Average Velocity of Piston [m/s]: {np.mean(np.abs(ys[3]))}")
    print(
        f"Root Mean Square Velocity of Piston [m/s]: {np.sqrt(np.mean(np.power(ys[3], 2)))}",
        end="\n\n",
    )

    # plt.hist(ts, bins=100)
    ## this histogram is to show density of time steps

    fig, ax = plt.subplots(2, 1, sharex=True)
    fig.set_size_inches(10, 6)

    ax[0].set_title("Displacement")
    (line_z1,) = ax[0].plot(
        ts, ys[0], color="blue", label="z1", alpha=0.85, linewidth=0.65
    )
    (line_z2,) = ax[0].plot(
        ts, ys[1], color="red", label="z2", alpha=0.85, linewidth=0.65
    )
    (line_z_diff,) = ax[0].plot(
        ts, ys[0] - ys[1], color="green", label="z1 - z2", alpha=0.85, linewidth=0.65
    )
    ax[0].legend(
        handles=[line_z1, line_z2, line_z_diff],
        loc="upper right",
        bbox_to_anchor=(1, 1),
    )

    ax[1].set_title("Velocity")
    (line_v1,) = ax[1].plot(
        ts, ys[2], color="blue", label="v1", alpha=0.85, linewidth=0.65
    )
    (line_v2,) = ax[1].plot(
        ts, ys[3], color="red", label="v2", alpha=0.85, linewidth=0.65
    )
    ax[1].set_xlabel("Time (s)")
    ax[1].legend(handles=[line_v1, line_v2], loc="upper right", bbox_to_anchor=(1, 1))

    fig, ax = plt.subplots(2, 1, sharex=True)  # force graphs
    fig.set_size_inches(10, 6)
    ax[0].plot(ts, force, color="blue", label="Fpto", alpha=0.85, linewidth=0.65)
    ax[0].set_title("Force of Power Take Off [kN]")
    ax[1].plot(
        ts, avg_force, color="red", label="Average Fpto", alpha=0.85, linewidth=0.65
    )
    ax[1].set_title("Average Force of Power Take Off [kN]")
    ax[1].set_xlabel("Time (s)")
    ax[0].legend(loc="upper right", bbox_to_anchor=(1, 1))
    ax[1].legend(loc="upper right", bbox_to_anchor=(1, 1))

    fig = plt.figure()
    fig.set_size_inches(10, 6)
    plt.scatter(zr, force, color="green", alpha=1.0, linewidth=0.65, marker="*")
    plt.xlabel("Extension [m]")
    if rhs.state.linear:
        plt.ylabel("Force F_{pto} = b_{2}(v_1 - v_2) [kN]")
    else:
        plt.ylabel("Force F_{pto} = F_{fric}tanh(100(v_1 - v_2)) + c(z_1 - z_2) [kN]")

    fig = plt.figure()
    fig.set_size_inches(10, 6)
    plt.plot(ts, power, color="blue", alpha=0.85, linewidth=0.65)
    plt.xlabel("Time [s]")
    plt.ylabel("P_{pto} = F_{pto}(v_1 - v_2) [kW]")
