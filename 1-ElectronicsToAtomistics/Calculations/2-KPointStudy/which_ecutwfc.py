import matplotlib.pyplot as plt
import numpy as np
import os

plot_first = True
cutoff_energies = np.arange(30, 135, 15)
kpoint = 8
legend_labels = [str(eng) for eng in cutoff_energies]
filepath = os.path.dirname(os.path.abspath(__file__))

### equilibrium energy
fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.set(
    xlabel=("K-Points"),
    ylabel=("Equilibrium Energy per Atom [" + r"$eV$]")
)
ax2.yaxis.set_ticks_position("right")
ax2.yaxis.set_label_position("right")
# ax2.set_yscale("symlog")
ax2.set(
    xlabel=("K-Points"),
    ylabel=("Percent Error [%]")
)
esub_normal = float(open(f"{filepath}/{cutoff_energies[-1]}/EvsK.dat").readlines()[-1].split(" ")[-1])
for ecutwfc in cutoff_energies:
    with open(f"{filepath}/{ecutwfc}/EvsK.dat", "r") as f:
        lines = f.readlines()
        K = np.zeros((len(lines), 1))
        E = np.zeros((len(lines), 1))
        N = np.zeros((len(lines)-1, 1))
        i = 0
        for line in lines:
            line_split = line.split(" ")
            if "\n" in line_split[1]:
                line_split[1] = line_split[1][:-1]
            K[i] = int(line_split[0])
            E[i] = float(line_split[1])
            if ecutwfc=="60":
                if K[i]==kpoint: print(E[i])
            if i > 0:
                N[i-1] = np.abs((E[i] - esub_normal)/((E[i] + esub_normal)/2))*100
            i += 1
    ax1.plot(K if plot_first else K[1:], E if plot_first else E[1:])
    ax2.plot(K[1:] if plot_first else K[2:], N if plot_first else N[1:])

    if ecutwfc==cutoff_energies[-1]: print(E[-1])

# Create the legend
fig.legend([ax1, ax2], # The line objects
    labels=legend_labels, # The labels for each line
    loc="lower center", # Position of legend
    ncol=len(legend_labels),
    borderaxespad=0.1, # Small spacing around legend box
    title="Cutoff Energies" # Title for the legend
)
plt.subplots_adjust(bottom=0.225)
plt.savefig("./energy_versus_kpoint.png")



### lattice parameter
fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.set(
    xlabel=("K-Points"),
    ylabel=("Equilibrium Lattice Parameter [" + r"$\AA$]")
)
ax2.yaxis.set_ticks_position("right")
ax2.yaxis.set_label_position("right")
# ax2.set_yscale("symlog")
ax2.set(
    xlabel=("K-Points"),
    ylabel=("Percent Error [%]")
)
alat_normal = float(open(f"{filepath}/{cutoff_energies[-1]}/AvsK.dat").readlines()[-1].split(" ")[-1])
for ecutwfc in cutoff_energies:
    with open(f"{filepath}/{ecutwfc}/AvsK.dat", "r") as f:
        lines = f.readlines()
        K = np.zeros((len(lines), 1))
        A = np.zeros((len(lines), 1))
        N = np.zeros((len(lines)-1, 1))
        i = 0
        for line in lines:
            line_split = line.split(" ")
            if "\n" in line_split[1]:
                line_split[1] = line_split[1][:-1]
            K[i] = int(line_split[0])
            A[i] = float(line_split[1])
            if ecutwfc=="60":
                if K[i]==kpoint: print(A[i])
            if i > 0:
                N[i-1] = np.abs((A[i] - alat_normal)/((A[i] + alat_normal)/2))*100
            i += 1
    ax1.plot(K if plot_first else K[1:], A if plot_first else A[1:])
    ax2.plot(K[1:] if plot_first else K[2:], N if plot_first else N[1:])

    if ecutwfc==cutoff_energies[-1]: print(A[-1])

# Create the legend
fig.legend([ax1, ax2], # The line objects
    labels=legend_labels, # The labels for each line
    loc="lower center", # Position of legend
    ncol=len(legend_labels),
    borderaxespad=0.1, # Small spacing around legend box
    title="Cutoff Energies" # Title for the legend
)
plt.subplots_adjust(bottom=0.225)
plt.savefig("./alat_versus_kpoint.png")



### bulk modulus
fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.set(
    xlabel=("K-Points"),
    ylabel=("Bulk Modulus (" + r"$G$) [" + r"$GPa$]")
)
ax2.yaxis.set_ticks_position("right")
ax2.yaxis.set_label_position("right")
# ax2.set_yscale("symlog")
ax2.set(
    xlabel=("K-Points"),
    ylabel=("Percent Error [%]")
)
bulk_normal = float(open(f"{filepath}/{cutoff_energies[-1]}/GvsK.dat").readlines()[-1].split(" ")[-1])/10
for ecutwfc in cutoff_energies:
    with open(f"{filepath}/{ecutwfc}/GvsK.dat", "r") as f:
        lines = f.readlines()
        K = np.zeros((len(lines), 1))
        G = np.zeros((len(lines), 1))
        N = np.zeros((len(lines)-1, 1))
        i = 0
        for line in lines:
            line_split = line.split(" ")
            if "\n" in line_split[1]:
                line_split[1] = line_split[1][:-1]
            K[i] = int(line_split[0])
            G[i] = float(line_split[1])/10
            if i > 0:
                N[i-1] = np.abs((G[i] - bulk_normal)/((G[i] + bulk_normal)/2))*100
            i += 1
    ax1.plot(K if plot_first else K[1:], G if plot_first else G[1:])
    ax2.plot(K[1:] if plot_first else K[2:], N if plot_first else N[1:])

    if ecutwfc==cutoff_energies[-1]: print(G[-1])

# Create the legend
fig.legend([ax1, ax2], # The line objects
    labels=legend_labels, # The labels for each line
    loc="lower center", # Position of legend
    ncol=len(legend_labels),
    borderaxespad=0.1, # Small spacing around legend box
    title="Cutoff Energies" # Title for the legend
)
plt.subplots_adjust(bottom=0.225)
plt.savefig("./bulk_versus_kpoint.png")



### convergence time
fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.set(
    xlabel=("K-Points"),
    ylabel=("Time to Convergence [s]")
)
ax2.yaxis.set_ticks_position("right")
ax2.yaxis.set_label_position("right")
# ax2.set_yscale("symlog")
ax2.set(
    xlabel=("K-Points"),
    ylabel=("Iterations")
)
time_normal = int(open(f"{filepath}/{cutoff_energies[-1]}/TvsK.dat").readlines()[-1].split(" ")[1])
for ecutwfc in cutoff_energies:
    with open(f"{filepath}/{ecutwfc}/TvsK.dat", "r") as f:
        lines = f.readlines()
        K = np.zeros((len(lines), 1))
        T = np.zeros((len(lines), 1))
        I = np.zeros((len(lines), 1))
        N = np.zeros((len(lines)-1, 1))
        i = 0
        for line in lines:
            line_split = line.split(" ")
            if "\n" in line_split[-1]:
                line_split[-1] = line_split[-1][:-1]
            K[i] = int(line_split[0])
            T[i] = float(line_split[1])
            I[i] = int(line_split[-1])
            if i > 0:
                N[i-1] = np.abs((T[i] - time_normal)/((T[i] + time_normal)/2))*100
            i += 1
    ax1.plot(K if plot_first else K[1:], T if plot_first else T[1:])
    ax2.plot(K if plot_first else K[1:], I if plot_first else I[1:])

# Create the legend
fig.legend([ax1, ax2], # The line objects
    labels=legend_labels, # The labels for each line
    loc="lower center", # Position of legend
    ncol=len(legend_labels),
    borderaxespad=0.1, # Small spacing around legend box
    title="Cutoff Energies" # Title for the legend
)
plt.subplots_adjust(bottom=0.225)
plt.savefig("./convergence_versus_kpoint.png")



# ### EvsA
# fig, (ax1, ax2) = plt.subplots(1, 2)
# ax1.set(
#     xlabel=("Lattice Parameter [" + r"$\AA$]"),
#     ylabel=("Equilibrium Energy per Atom [" + r"$eV$]")
# )
# ax2.yaxis.set_ticks_position("right")
# ax2.yaxis.set_label_position("right")
# # ax2.set_yscale("symlog")
# ax2.set(
#     xlabel=("Lattice Parameter Volume [" + r"$\AA^{3}$]"),
#     ylabel=("Equilibrium Energy per Atom [" + r"$eV$]")
# )
# do_full, do_bcc = False, True
# for ecutwfc in cutoff_energies:
#     with open(f"{filepath}/{ecutwfc}/EvsA.{ecutwfc}.{kpoint}" if do_full else f"../EvsA-{'bcc' if do_bcc else 'fcc'}-Calibration", "r") as f:
#         lines = f.readlines()
#         A = np.zeros((len(lines), 1))
#         E = np.zeros((len(lines), 1))
#         i = 0
#         for line in lines:
#             line_split = line.split(" ")
#             if "\n" in line_split[1]:
#                 line_split[1] = line_split[1][:-1]
#             A[i] = float(line_split[0])
#             E[i] = float(line_split[1])
#             i += 1
#     ax1.plot(A if plot_first else A[1:], E if plot_first else E[1:])
# for ecutwfc in cutoff_energies:
#     with open(f"{filepath}/{ecutwfc}/EvsA.{ecutwfc}.{kpoint}" if do_full else f"../EvsA-{'bcc' if do_bcc else 'fcc'}-Calibration", "r") as f:
#         lines = f.readlines()
#         V = np.zeros((len(lines), 1))
#         E = np.zeros((len(lines), 1))
#         i = 0
#         for line in lines:
#             line_split = line.split(" ")
#             if "\n" in line_split[1]:
#                 line_split[1] = line_split[1][:-1]
#             V[i] = float(line_split[0])
#             E[i] = float(line_split[1])
#             i += 1
#     ax2.plot(V if plot_first else V[1:], E if plot_first else E[1:])

# # Create the legend
# fig.legend([ax1, ax2], # The line objects
#     labels=legend_labels, # The labels for each line
#     loc="lower center", # Position of legend
#     ncol=len(legend_labels),
#     borderaxespad=0.1, # Small spacing around legend box
#     title="Cutoff Energies" # Title for the legend
# )
# plt.subplots_adjust(bottom=0.225)
# plt.savefig("{filepath}/EvsA_for_energy.png")