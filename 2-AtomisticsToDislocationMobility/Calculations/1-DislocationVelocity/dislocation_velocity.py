from engineering_notation import EngNumber as engr
import joby_m_anthony_iii.numerical_methods as nm
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import sympy as sp
import sys
import warnings
# warnings.filterwarnings("ignore", module="matplotlib\..*")
warnings.filterwarnings(action="ignore", module="numpy/*", message="Empty input file", category=UserWarning)

# 194.8
dt, width = 0.500, 40 * 5   # lattice distance * unit distance (approx)
offset      = 40.0
L_buff      = 0.0 + offset
R_buff      = width - offset
midpoint    = width / 2

lattice_parameter = 2.781e-10 # m
slip_direction = np.array([1, 1, 1])
burger_vec = lattice_parameter*(nm.norms(slip_direction).l_two()/2) # (np.sqrt(3)/2)
print(burger_vec)

# must be same examined in `rescale_commands.sh`
# if not a range, then define as blank string: TEMP = ""
TEMP = np.array([300]) # np.arange(150, 500, 50) # K
# if not a range, then define as blank string: SIGMA = ""
SIGMA = np.arange(25, 325, 25) # MPa
for s in np.array([700, 800, 900, 1000, 1100, 1200]):
    SIGMA = np.append(SIGMA, s)

skip = 1
columns = ["Particle Identifier", "X", "Y", "Z", "Centrosymmetry", "R", "G", "B", "nan", "Time"]

do_every = done_every = False
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
                x_ave, lap_count = 0, 0
                data_dir = f"./PositionFrameData/{temp}/{sigma:.3f}"
                if do_every:
                    particles = {}
                    for i in range(0, len(next(os.walk(data_dir))[2]), skip):
                        with open(f"{data_dir}/lammps.{i}.xyz") as g:
                            try:
                                py_array = np.genfromtxt(g.readlines()[2:])
                                data = pd.DataFrame(py_array, columns=columns)
                                particles[i] = np.array(data["Particle Identifier"])
                            except ValueError:
                                particles[i] = np.array([])

                    particle_ids = {}
                    for frame in particles:
                        ids = particles[frame]
                        for id in ids:
                            particle_ids[int(id)] = []
                    for frame in particles:
                        ids = particles[frame]
                        for id in ids:
                            particle_ids[int(id)].append(frame)
                    foo = 0
                    for ids in particle_ids:
                        if len(particle_ids[ids]) == 1:
                            foo += 1
                        else: foo -= 1
                    if foo == len(particle_ids): done_every = False
                    else: done_every = True

                velocities, offsets = [], []
                if done_every:
                    for id in particle_ids:
                        i, time, position = 0, [], []
                        for frame in particle_ids[int(id)]:
                            with open(f"{data_dir}/lammps.{frame}.xyz") as g:
                                py_array = np.genfromtxt(g.readlines()[2:])
                                data = pd.DataFrame(py_array, columns=columns)
                                time.append(frame*dt)#np.average(data["Time"]))
                                # position.append(np.average(data[["X"]]))
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
                        time, position = np.array(time), np.array(position)/10 # ps, nm
                        # i, position2 = 0, position - position[0]
                        # for i in range(1, len(position)):
                        #     position2[i] = position2[i-1] + np.abs(position[i] - position[i-1])
                        # position = position2
                        ax_sigma.scatter(time, position)

                        try:
                            position_expr = nm.least_squares.linear(time, position, 1)[0]
                            other_pos = np.array(position_expr(t) for t in time)
                            # ax_sigma.plot(time, other_pos)

                            if "t" != str(position_expr(sp.Symbol("t")))[-1]:
                                velocity = float(str(position_expr(sp.Symbol("t"))).split(" ")[0].split("*")[0])
                                offset = float(str(position_expr(sp.Symbol("t"))).split(" ")[-1])
                            else :
                                velocity = float(str(position_expr(sp.Symbol("t"))).split(" ")[-1].split("*")[0])
                                offset = float(str(position_expr(sp.Symbol("t"))).split(" ")[0])
                            velocities.append(velocity)
                            offsets.append(offset)
                        except np.linalg.LinAlgError: pass
                else:
                    i, time, position = 0, [], []
                    for frame in range(0, len(next(os.walk(data_dir))[2]), skip):
                        with open(f"{data_dir}/lammps.{frame}.xyz") as g:
                            try:
                                py_array = np.genfromtxt(g.readlines()[2:])
                                if len(py_array) > 0:
                                    data = pd.DataFrame(py_array, columns=columns)
                                    time.append(frame*dt)#np.average(data["Time"]))
                                    # position.append(np.average(data[["X"]]))
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
                    # i, position2 = 0, position - position[0]
                    # for i in range(1, len(position2)):
                    #     position2[i] = position2[i-1] + np.abs(position[i] - position[i-1])
                    # position = position2
                    ax_sigma.scatter(time, position)#, label="Raw")

                    try:
                        position_expr = nm.least_squares.linear(time, position, 1)[0]
                        other_pos = np.array([position_expr(t) for t in time])
                        # ax_sigma.plot(time, other_pos, "k", label="Homogeneous Fit")

                        if "t" != str(position_expr(sp.Symbol("t")))[-1]:
                            velocity = float(str(position_expr(sp.Symbol("t"))).split(" ")[0].split("*")[0])
                            offset = float(str(position_expr(sp.Symbol("t"))).split(" ")[-1])
                        else :
                            velocity = float(str(position_expr(sp.Symbol("t"))).split(" ")[-1].split("*")[0])
                            offset = float(str(position_expr(sp.Symbol("t"))).split(" ")[0])
                        velocities.append(velocity)
                        offsets.append(offset)
                    except np.linalg.LinAlgError: pass

                velocity = np.average(velocities)
                f.write(f"{sigma:.9f} {velocity*1e3:.9f}\n") # MPa m/s
                offset = np.average(offsets)

                try:
                    pft_string = f"{engr(velocity, 3)}*t + {engr(offset, 3)}"
                    # ax_sigma.text(
                    #     x=ax_sigma.get_xlim()[1]*0.65,
                    #     y=ax_sigma.get_ylim()[1]*0.05,
                    #     s=pft_string
                    # )
                except ValueError:
                    pft_string = f"Could not round..."
                    # ax_sigma.text(
                    #     x=ax_sigma.get_xlim()[1]*0.65,
                    #     y=ax_sigma.get_ylim()[1]*0.05,
                    #     s=pft_string
                    # )
                time = np.arange(0, 100, dt*skip)
                other_pos = velocity*time + offset
                ax_sigma.plot(time, other_pos, label=f"{sigma} MPa: {pft_string}")
            ax_sigma.legend()
            fig_sigma.savefig(f"position_for_time-{temp}.pdf")
        with open(f"dislocation_velocity-{temp}", "r") as f:
            data = pd.DataFrame(np.genfromtxt(f.readlines())[1:], columns=["Stress", "Velocity"])
            drag = burger_vec/float(str(nm.least_squares.linear(data['Stress'][:12], data['Velocity'][:12], 1)[0](sp.Symbol('t'))).split('*')[0])*1e6 # Pa-s
            DRAG[f"{temp}"] = 1/drag
            ax_temp.scatter(data["Stress"], data["Velocity"], label=f"{temp} K" + r"$\rightarrow$ B = " + f"{drag:3e}" + r" Pa-s")
        ax_temp.legend()
        fig_temp.savefig(f"velocity_for_stress.pdf")
finally:
    sys.stdout.close()
    sys.stdout = old_stdout
    print(DRAG)