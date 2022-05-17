from engineering_notation import EngNumber as engr
import joby_m_anthony_iii.numerical_methods as nm
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import sympy as sp

dt = 500
# must be same examined in `rescale_commands.sh`
# if not a range, then define as blank string: TEMP = ""
TEMP = np.array([300]) # np.arange(150, 550, 50)
# if not a range, then define as blank string: SIGMA = ""
SIGMA = np.arange(25, 325, 25)

skip = 3
columns = ["Particle Identifier", "X", "Y", "Z", "Centrosymmetry", "R", "G", "B", "nan", "Time"]

do_every = False

for temp in TEMP:
    with open(f"dislocation_velocity-{int(temp)}", "w") as f:
        f.write(f"# AppliedShearStress(MPa) DislocationVelocity(m/s)")
        for sigma in SIGMA:
            fig_sigma, ax_sigma = plt.subplots(1, 1)
            ax_sigma.set(
                title="Dislocation Position for Time",
                xlabel="Time (t) [" + r"$ps$]",
                ylabel="Displacement (x) [" + r"$\AA$]",
            )
            data_dir = f"./PositionFrameData/{int(temp)}/{int(sigma)}"
            if do_every:
                particles = {}
                for i in range(0, len(next(os.walk(data_dir))[2]), skip):
                    with open(f"{data_dir}/lammps.{i}.xyz") as f:
                        try:
                            py_array = np.genfromtxt(f.readlines()[2:])
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

            velocities, offsets = [], []
            if do_every:
                for id in particle_ids:
                    i, time, position = 0, [], []
                    for frame in particle_ids[int(id)]:
                        with open(f"{data_dir}/lammps.{frame}.xyz") as f:
                            py_array = np.genfromtxt(f.readlines()[2:])
                            data = pd.DataFrame(py_array, columns=columns)
                            time.append(np.average(data["Time"]))
                            position.append(np.average(data[["X"]]))
                            i += 1
                    time, position = np.array(time), np.array(position)
                    i, position2 = 0, position - position[0]
                    for i in range(1, len(position)):
                        position2[i] = position[i-1] + np.abs(position[i] - position[i-1])
                    position = position2
                    ax_sigma.scatter(time, position)

                    try:
                        position_expr = nm.least_squares.linear(time, position, 1)[0]
                        other_pos = []
                        for t in time: other_pos.append(position_expr(t))
                        # ax_sigma.plot(time, other_pos)

                        velocity = float(str(position_expr(sp.Symbol("t"))).split(" ")[0].split("*")[0])
                        velocities.append(velocity)
                        offset = float(str(position_expr(sp.Symbol("t"))).split(" ")[-1])
                        offsets.append(offset)
                    except np.linalg.LinAlgError: pass
            else:
                i, time, position = 0, [], []
                for frame in range(0, 201, skip):
                    with open(f"{data_dir}/lammps.{frame}.xyz") as f:
                        py_array = np.genfromtxt(f.readlines()[2:])
                        data = pd.DataFrame(py_array, columns=columns)
                        time.append(np.average(data["Time"]))
                        position.append(np.average(data[["X"]]))
                        i += 1
                time, position = np.array(time), np.array(position)
                i, position2 = 0, position - position[0]
                for i in range(1, len(position)):
                    position2[i] = position[i-1] + np.abs(position[i] - position[i-1])
                position = position2
                ax_sigma.scatter(time, position, label="Raw")

                try:
                    position_expr = nm.least_squares.linear(time, position, 1)[0]
                    other_pos = []
                    for t in time: other_pos.append(position_expr(t))
                    ax_sigma.plot(time, other_pos, label="Fit")

                    velocity = float(str(position_expr(sp.Symbol("t"))).split(" ")[0].split("*")[0])
                    velocities.append(velocity)
                    offset = float(str(position_expr(sp.Symbol("t"))).split(" ")[-1])
                    offsets.append(offset)
                except np.linalg.LinAlgError: pass

            velocity = np.average(velocities)
            f.write(f"{sigma} {velocity}")
            offset = np.average(offsets)

            print(velocity)

            if do_every:
                time = np.linspace(0, 100e3, dt)
                other_pos = velocity*time + offset
                ax_sigma.plot(time, other_pos, "k", label="Fit")

            # plt.show()
            ax_sigma.text(
                x=40e3,
                y=150,
                s=f"{engr(velocity)}*t + {engr(offset)}"
            )
            ax_sigma.legend()
            fig_sigma.savefig(f"position_for_time-{int(temp)}/{int(sigma)}.pdf")