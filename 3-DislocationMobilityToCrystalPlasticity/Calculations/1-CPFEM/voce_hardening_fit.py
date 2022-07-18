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
# warnings.filterwarnings(action="ignore", module="numpy/*", message="Empty input file", category=UserWarning)

lattice_parameter = 2.781e-10 # m
slip_direction = np.array([1, 1, 1])
burger_vec = lattice_parameter*(nm.Norm(slip_direction).l_two()/2) # (np.sqrt(3)/2)
# print(burger_vec)
mu = 195.343e3 # 117e3 # MPa

# must be same examined in `rescale_commands.sh`
# if not a range, then define as blank string: TEMP = ""
TEMP = np.array([300]) # np.arange(150, 500, 50) # K
# if not a range, then define as blank string: SIGMA = ""
STRAIN = np.array([1000]) # np.arange(25, 325, 25) # MPa
# for s in np.array([700, 800, 900, 1000, 1100, 1200]):
#     SIGMA = np.append(SIGMA, s)

skip = int(1.5e2)
columns = ["timenow", "disDensity", "Stress", "Strain", "S1", "S2", "S3", "S23", "S31", "S12", "p1", "p2", "p3", "p23", "p31",  "p12", "jogs", "junctions", "CrossSlip"]

do_monitor = False

H0 = np.arange(0.4, 0.42, 0.0001)*1e6 # Pa
# H0 = np.arange(0, 101, 1)*1e6 # Pa
C = np.array([1000]) # np.arange(73.5, 73.61, 0.01) # 1/s
# C = np.arange(0, 101, 1) # 1/s
Kappa = lambda ks, k0, h0, c, t: ks - (ks - k0)*np.exp(-h0/(ks - k0)*c*t)

# old_stdout = sys.stdout
# sys.stdout = open(os.devnull, "w")
# try:
fig, (ax_kappa, ax_alpha) = plt.subplots(1, 2)
# fig, ax_kappa = plt.subplots(1, 1)
ax_kappa.set(
    title="Stress-Strain @ 300 K",
    # xscale="symlog",
    xlabel=r"Normalized True Strain ($\epsilon$)",
    # yscale="symlog",
    ylabel=r"True Stress ($\sigma$) [$MPa$]",
)
ax_alpha.set(
    title="Average Value of\nJunction Strength @ 300 K",
    # xlim=(0, 0.0002),
    # xscale="symlog",
    xlabel=r"Time ($t$) [$\mu s$]",
    # yscale="symlog",
    ylabel=r"Junction Strength ($\alpha$)",
)
ax_alpha.yaxis.set_label_position("right")
if do_monitor:
    # MDDP_BCC_HW2/MDDP_Zip/Windows/
    with open(f"../Homework 2/MDDP_BCC_HW2/MDDP_Zip/Windows/DDtimeResults.out", "r") as f:
    # with open(f"./input", "r") as f:
        data = pd.DataFrame(np.genfromtxt(f.readlines()[4:]), columns=columns)
        # data.plot(
        #     kind="scatter",
        #     x="Strain",
        #     y="Stress",
        #     color="blue",
        #     ax_kappa=plt.gca())
        ax_kappa.scatter(
            np.array([data["Strain"][i] for i in range(0, len(data), skip)]),#[:int(1e1)],
            np.array([data["Stress"][i] for i in range(0, len(data), skip)])/1e6,#[:int(1e1)]/1e6
            ) # label=f"{strain} " + r"$\frac{1}{s}$")
    # ax_kappa.legend()
    fig.savefig(f"./ENGR851_JobyAnthonyIII_Homework3/stress_strain-monitor.png")
else:
    for temp in TEMP:
        for strain in STRAIN:
            # # MDDP_BCC_HW2/MDDP_Zip/Windows/
            # with open(f"../Homework 2/MDDP_BCC_HW2/MDDP_Zip/Windows/{temp}/{strain}/ResultPerSlip", "r") as f:
            #     data = f.readlines()
            #     search = " variables= \"SlipSystem\""
            #     line, lines = 0, []
            #     for dat in data:
            #         if search in dat:
            #             lines.append(line)
            #         line += 1
            # f.close()
            # lines = np.array(lines)
            with open(f"../Homework 2/MDDP_BCC_HW2/MDDP_Zip/Windows/{temp}/{strain}/DDtimeResults.out", "r") as f:
            # with open(f"./input", "r") as f:
                data = pd.DataFrame(np.genfromtxt(f.readlines()[4:]), columns=columns)
                # print(len(lines), len(data))
                # lines = lines[:len(data)+1]
                # with open(f"../Homework 2/MDDP_BCC_HW2/MDDP_Zip/Windows/{temp}/{strain}/ResultPerSlip", "r") as g:
                #     slip = g.readlines()
                #     i, rho_f = 0, np.zeros(len(lines)) # int, 1/m2
                #     for line in lines[:-2]:
                #         rho_f[i] = np.average(np.genfromtxt(slip[line+2:lines[i+1]-1])[:,1])
                #         i += 1
                #     rho_f[i] = np.average(np.genfromtxt(slip[lines[-2]+2:lines[-1]-1])[:,1])
                # g.close()
                time = np.array([data["timenow"][i] for i in range(0, len(data), skip)])
                eps = np.array([data["Strain"][i] for i in range(0, len(data), skip)])#/data["Strain"].values[-1],#[:int(1e1)],
                sigma = np.array([data["Stress"][i] for i in range(0, len(data), skip)])/1e6#[:int(1e1)]/1e6
                # rho_fp = np.array([rho_f[i] for i in range(0, len(data), skip)])
                # rho_f = rho_fp
                rho_f = np.array([data["disDensity"][i] for i in range(0, len(data), skip)])
                # data.plot(
                #     kind="scatter",
                #     x="Strain",
                #     y="Stress",
                #     color="blue",
                #     ax_kappa=plt.gca())
                ax_kappa.scatter(
                    eps,
                    sigma,
                    label=r"$\epsilon_{f}$ = " + f"{eps[-1]} @ " + r"$\dot{\epsilon}$ = " + f"{strain} " + r"$\frac{1}{s}$"
                )
            ERRORS = {}
            for j in range(0, 15, 1):
                ERRORS[f"{j}"] = {}
                for h0 in H0:
                    ERRORS[f"{j}"][f"{float(h0/1e6)}"] = {}
                    for c in C:
                        i, kappa = 0, np.zeros(len(time))
                        error = np.zeros_like(kappa)
                        for t in time:
                            kappa[i] = Kappa(np.average(sigma[-j:]), sigma[0], h0, c, t)
                            if i > 0:
                                error[i] = np.abs((kappa[i] - sigma[i])/sigma[i])
                            i += 1
                        ERRORS[f"{j}"][f"{float(h0/1e6)}"][f"{c}"] = np.sum(error)
            min_error = ["", "", "", 1e6]
            for j in range(0, 15, 1):
                for h0 in ERRORS[f"{j}"]:
                    for c in ERRORS[f"{j}"][h0]:
                        if ERRORS[f"{j}"][h0][c] < min_error[-1]:
                            min_error = [j, h0, c, ERRORS[f"{j}"][h0][c]]
            print(min_error)
            j, h0, c, error = min_error
            i, kappa = 0, np.zeros(len(time))
            alpha = np.zeros_like(kappa)
            for t in time:
                kappa[i] = Kappa(np.average(sigma[-j:]), sigma[0], float(h0)*1e6, float(c), t)
                alpha[i] = kappa[i]/(mu*burger_vec*np.sqrt(rho_f[i]))
                i += 1
            print(kappa[0])
            # foo = nm.MultiVariableIteration(Kappa, np.array([h0, C]), kappa)
            ax_kappa.scatter(
                eps,
                kappa,
                label=r"Voce Equation, $\kappa = \kappa_{s} - (\kappa_{s} - \kappa_{0})\exp(-\frac{h_{0}}{\kappa_{s} - \kappa_{0}}Ct)$"
            )
            ax_kappa.text(
                x=ax_kappa.get_xlim()[1]*0.65,
                y=ax_kappa.get_ylim()[1]*0.4,
                s=r"$\kappa_{s}$ = " + f"{round(np.average(sigma[-j:]), 4)} MPa\n" + r"$\kappa_{0}$ = " + f"{round(sigma[0], 4)} MPa\n" + r"$C$ = " + f"{int(c)} " + r"$\frac{1}{s}$" + f"\n" + r"$h_{0}$ = " + f"{round(float(h0), 4)} MPa\n" + f"E = {round(error, 4)} " + r"$\frac{MPa}{MPa}$"
            )
            ax_alpha.scatter(
                time*1e3,
                alpha
            )
        ax_kappa.legend()
        fig.savefig(f"./ENGR851_JobyAnthonyIII_Homework3/stress_strain-{temp}.png")
        # finally:
        #     sys.stdout.close()
        #     sys.stdout = old_stdout
        #     print(DRAG)