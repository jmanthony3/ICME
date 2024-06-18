from engineering_notation import EngNumber as engr
import joby_m_anthony_iii.numerical_methods as nm
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import sympy as sp
import sys
import warnings
warnings.filterwarnings(action="ignore", module="numpy/*", message="Empty input file", category=UserWarning)

filepath = os.path.dirname(os.path.abspath(__file__))

print("")

dt, width = 0.500, 40 * 5   # lattice distance * unit distance (approx)
offset      = 40.0
L_buff      = 0.0 + offset
R_buff      = width - offset
midpoint    = width / 2

lattice_parameter = 2.781e-10 # m
slip_direction = np.array([1, 1, 1])
burger_vec = lattice_parameter*(nm.Norm(slip_direction).l_two()/2)
print("Magnitude of Burger's vector, ||b|| = " + f"{burger_vec} m", end="\n\n")

# must be same examined in `rescale_commands.sh`
TEMP = np.array([300]) # np.arange(150, 550, 50) # K
SIGMA = np.arange(25, 325, 25) # MPa
# for s in np.array([700, 800, 900, 1000, 1100, 1200]):
#     SIGMA = np.append(SIGMA, s)

skip = 1
columns = ["Particle Identifier", "X", "Y", "Z", "Centrosymmetry", "R", "G", "B", "nan", "Time"]

DRAG = {}

old_stdout = sys.stdout
sys.stdout = open(os.devnull, "w")
try:
    for temp in TEMP:
        fig_temp, ax_temp = plt.subplots(1, 1)
        ax_temp.set(
            title="Dislocation Velocity for Stress",
            xlabel="Stress (" + r"$\tau$) [$MPa$]",
            ylabel="Velocity (" + r"$v_{x}$) [$\frac{m}{s}$]",
        )
        with open(f"dislocation_velocity-{temp}", "w") as f:
            fig_sigma, ax_sigma = plt.subplots(1, 1)
            ax_sigma.set(
                title="Dislocation Position for Time",
                xlabel="Time (t) [" + r"$ps$]",
                ylabel="Displacement (x) [" + r"$nm$]",
            )
            f.write(f"Stress Velocity\n")
            for sigma in SIGMA:
                print(f"Looking at {temp:3} K, {sigma:5.3f} MPa...", end="\r")
                x_ave, lap_count = 0, 0
                data_dir = f"{filepath}/PositionFrameData/{temp}/{sigma:.3f}"

                velocities, offsets = [], []
                i, time, position = 0, [], []
                for frame in range(0, len(next(os.walk(data_dir))[2]), skip):
                    with open(f"{data_dir}/lammps.{frame}.xyz") as g:
                        try:
                            py_array = np.genfromtxt(g.readlines()[2:])
                            if len(py_array) > 0:
                                data = pd.DataFrame(py_array, columns=columns)
                                time.append(frame*dt)
                                x_store = data[["X"]].values
                                # check to see if the dsl cluster is on wraping both sides
                                if min(x_store) < L_buff and R_buff < max(x_store):
                                    for j in range(len(x_store)):
                                        #values on left side (wrapped around)
                                        if x_store[j] < midpoint:
                                            x_store[j] = x_store[j] + width
                                x_ave_prev = x_ave
                                x_ave = np.average(x_store) + width*lap_count
                                # check if a lap just happened
                                if x_ave < (x_ave_prev - midpoint):
                                    lap_count += 1
                                    x_ave += width
                                position.append(x_ave)
                        except ValueError: pass
                        except UserWarning: pass
                    i += 1
                time, position = np.array(time), np.array(position)/10 # ps, nm
                ax_sigma.scatter(time, position)#, label="Raw")

                try:
                    position_expr = nm.least_squares.linear(time, position, 1)[0]
                    other_pos = np.array([position_expr(t) for t in time])
                    # ax_sigma.plot(time, other_pos, "k", label="Homogeneous Fit")

                    if "t" != str(position_expr(sp.Symbol("t")))[-1]:
                        velocity = float(str(position_expr(sp.Symbol("t"))).split(" ")[0].split("*")[0])
                        offset = float(str(position_expr(sp.Symbol("t"))).split(" ")[-1])
                    else:
                        velocity = float(str(position_expr(sp.Symbol("t"))).split(" ")[-1].split("*")[0])
                        offset = float(str(position_expr(sp.Symbol("t"))).split(" ")[0])
                    velocities.append(velocity)
                    offsets.append(offset)
                except np.linalg.LinAlgError: pass

                velocity = np.average(velocities)
                f.write(f"{sigma:.9f} {velocity*1e3:.9f}\n") # MPa m/s
                offset = np.average(offsets)

                try:
                    pft_string = f"{round(velocity, 3)}*t + {round(offset, 3)}"
                except ValueError:
                    pft_string = f"Could not round..."
                time = np.arange(0, 100, dt*skip)
                other_pos = velocity*time + offset
                ax_sigma.plot(time, other_pos, label=f"{int(sigma)}".zfill(4) + f" MPa | {pft_string}".ljust(15))
                print(f"Closing {temp} K, {sigma:.3f} MPa...", end="\r")
            ax_sigma.text(
                x=ax_sigma.get_xlim()[1]*0.175,
                y=ax_sigma.get_ylim()[1]*-0.025,
                s=r"The $v$'s in legend equations are dislocation velocities."
            )
            ax_sigma.legend(
                title=r"$x(t) = v[\frac{nm}{ps}]*t[ps] + [nm]$"
            )
            fig_sigma.savefig(f"position_for_time-{temp}.png")
        with open(f"dislocation_velocity-{temp}", "r") as f:
            data = pd.DataFrame(np.genfromtxt(f.readlines())[1:], columns=["Stress", "Velocity"])
            drag = burger_vec/float(str(nm.least_squares.linear(data['Stress'][:12], data['Velocity'][:12], 1)[0](sp.Symbol('t'))).split('*')[0])*1e6 # Pa-s
            DRAG[f"{temp}"] = (drag, 1/drag)
            ax_temp.scatter(data["Stress"], data["Velocity"], label=f"{temp} K" + r"$\rightarrow$ B = " + f"{drag:3e}" + r" Pa-s")
        ax_temp.legend()
        fig_temp.savefig(f"velocity_for_stress.png")
finally:
    sys.stdout.close()
    sys.stdout = old_stdout

    # Print the names of the columns.
    print("{:<15} | {:<30} | {:<35}".format("Temperature [K]", "Drag Coefficient (B) [Pa-s]", "Dislocation Mobility (Ms) [1/Pa-s]"))
    print("-"*16 + "|" + "-"*32 + "|" + "-"*36)

    # print each data item.
    for key, value in DRAG.items():
        coeff, mob = value
        print("{:<15} | {:<30} | {:<35}".format(key, coeff, mob))