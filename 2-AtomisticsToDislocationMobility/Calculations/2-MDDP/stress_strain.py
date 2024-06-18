import joby_m_anthony_iii.numerical_methods as nm
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import warnings
warnings.filterwarnings(action="ignore", module="numpy/*", message="Empty input file", category=UserWarning)

filepath = os.path.dirname(os.path.abspath(__file__))

print("")

dt, width = 500, 40 * 5   # lattice distance * unit distance (approx)
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
STRAIN = np.array([3], dtype=int) # np.arange(1, 7, 1) # 1/s
# for s in np.array([-3, -2, -1]):
#     STRAIN = np.append(STRAIN, s)
STRAIN = 10**STRAIN

skip = int(1.5e2)
columns = ["timenow", "disDensity", "Stress", "Strain", "S1", "S2", "S3", "S23", "S31", "S12", "p1", "p2", "p3", "p23", "p31",  "p12", "jogs", "junctions", "CrossSlip"]

do_monitor = True

if do_monitor:
    fig, ax = plt.subplots(1, 1)
    ax.set(
        title="Stress-Strain Curve",
        # xscale="symlog",
        xlabel="True Strain (" + r"$\epsilon$)",
        # yscale="symlog",
        ylabel="True Stress (" + r"$\sigma$) [$MPa$]",
    )
    with open(f"{filepath}/DDtimeResults.out", "r") as f:
        data = pd.DataFrame(np.genfromtxt(f.readlines()[4:]), columns=columns)
        ax.scatter(
            np.array([data["Strain"][i] for i in range(0, len(data), skip)]),
            np.array([data["Stress"][i] for i in range(0, len(data), skip)])/1e6, # MPa
            # label=f"{strain} " + r"$\frac{1}{s}$")
        )
    # ax.legend()
    fig.savefig(f"stress_strain-monitor.png")
else:
    for temp in TEMP:
        fig, ax = plt.subplots(1, 1)
        ax.set(
            title=f"Stress-Strain @ {temp}" + r"[$K$]",
            # xscale="symlog",
            xlabel="True Strain (" + r"$\epsilon$)",
            # yscale="symlog",
            ylabel="True Stress (" + r"$\sigma$) [$MPa$]",
        )
        for strain in STRAIN:
            with open(f"{filepath}/{temp}/{strain}/DDtimeResults.out", "r") as f:
                data = pd.DataFrame(np.genfromtxt(f.readlines()[4:]), columns=columns)
                ax.scatter(
                    np.array([data["Strain"][i] for i in range(0, len(data), skip)]),#/data["Strain"].values[-1],
                    np.array([data["Stress"][i] for i in range(0, len(data), skip)])/1e6, # MPa
                    label=r"$\epsilon_{f}$ = " + f"{data['Strain'].values[-1]} @ " + r"$\dot{\epsilon}$ = " + f"{strain} " + r"$\frac{1}{s}$")
        ax.legend()
        fig.savefig(f"{filepath}/stress_strain-{temp}.png")