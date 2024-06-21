#####################################################################################
# which_ecutwfc.jl                                                                  #
#                                                                                   #
# CHANGELOG                                                                         #
# 2024-06-17: Converted from Python3 to Julia and added this title block.           #
#                                                                                   #
# 2022-XX-XX: Original import from Dr. Heechen Cho.                                 #
#####################################################################################

using CSV
using DataFrames
using Plots; gr()

default(grid=false, framestyle=:box)

plot_first = true
cutoff_energies = range(30, 120; step=15)
kpoint = 8
legend_labels = ["$eng" for eng in cutoff_energies]
l = @layout [grid(1, 2); a{1.0w, 0.2h}]
# l = @layout [grid(1, 2); _ a{0.95w, 0.2h} _]
# l = @layout [grid(1, 2); grid(1, 3; widths=[1/9, 7/9, 1/9])]
# l = @layout [grid(1, 2); grid(1, 1)]
filepath = pwd()



ax2                 = plot( # percent difference (base)
    xlabel="K-Points",
    xlims=(0, 12),
    ylabel="Percent Error [%]",
    yguide_position=:right,
    # yscale=log10,
    ymirror=true,
    bottommargin=(1.0, :mm))
plot_EvsK_ax1       = plot( # equilibrium energy
    xlabel="K-Points",
    ylabel="Equilibrium Energy per Atom [\$eV\$]",
    xlims=(0, 12),
    bottommargin=(1.0, :mm))
plot_EvsK_ax2       = deepcopy(ax2)
plot_AvsK_ax1       = plot( # lattice parameter
    xlabel="K-Points",
    ylabel="Equilibrium Lattice Parameter [\$\\AA\$]",
    xlims=(0, 12),
    bottommargin=(1.0, :mm))
plot_AvsK_ax2       = deepcopy(ax2)
plot_GvsK_ax1       = plot( # bulk modulus
    xlabel="K-Points",
    ylabel="Bulk Modulus (\$G\$) [\$GPa\$]",
    xlims=(0, 12),
    bottommargin=(1.0, :mm))
plot_GvsK_ax2       = deepcopy(ax2)
plot_TvsK_ax1       = plot( # bulk modulus
    xlabel="K-Points",
    ylabel="Time to Convergence [\$s\$]",
    xlims=(0, 12),
    bottommargin=(1.0, :mm))
plot_TvsK_ax2       = deepcopy(ax2)
pb1 = plot(legend=false,grid=false,foreground_color_subplot=:white);
pb2 = plot(legend=false,grid=false,foreground_color_subplot=:white);
legend              = plot((1:length(legend_labels))',
    labels=reshape(legend_labels, (1, length(legend_labels))),
    fillcolors=palette(:default)[1:length(legend_labels)]',
    legendcolumns=div(length(legend_labels) + 1, 2),
    legendforeground_color=:black,
    legendposition=:bottom,
    legendtitle="Cutoff Energies",
    framestyle=:none,
    # leftmargin=(-2.5, :mm),
    # rightmargin=(-2.5, :mm),
    # topmargin=(-4, :mm),
    bottommargin=(-6, :mm),
)


esub_normal = last(CSV.read("$filepath/$(last(cutoff_energies))/EvsK.dat", DataFrame)[!, 2])
alat_normal = last(CSV.read("$filepath/$(last(cutoff_energies))/AvsK.dat", DataFrame)[!, 2])
bulk_normal = last(CSV.read("$filepath/$(last(cutoff_energies))/GvsK.dat", DataFrame)[!, 2])
time_normal = last(CSV.read("$filepath/$(last(cutoff_energies))/TvsK.dat", DataFrame)[!, 2])


for ecutwfc in collect(cutoff_energies)
    df_EvsK = CSV.read("$filepath/$ecutwfc/EvsK.dat", DataFrame)
    df_AvsK = CSV.read("$filepath/$ecutwfc/AvsK.dat", DataFrame)
    df_GvsK = CSV.read("$filepath/$ecutwfc/GvsK.dat", DataFrame)
    df_TvsK = CSV.read("$filepath/$ecutwfc/TvsK.dat", DataFrame)

    K = df_EvsK[!, 1]
    E = df_EvsK[!, 2]
    A = df_AvsK[!, 2]
    G = df_GvsK[!, 2]
    T = df_TvsK[!, 2]
    I = df_TvsK[!, 3]

    N_E = zeros(nrow(df_EvsK)-1)
    N_A = zeros(nrow(df_AvsK)-1)
    N_G = zeros(nrow(df_GvsK)-1)
    N_T = zeros(nrow(df_TvsK)-1)
    for (i, kp, e, a, g, t) in zip(1:1:nrow(df_EvsK), K, E, A, G, T)
        if ecutwfc == "60" && kp == kpoint
            println(e)
            println(a)
            println(g)
        end
        if i > 1
            # percent difference
            N_E[i-1] = 100abs((e - esub_normal) / ((e + esub_normal) / 2.))
            N_A[i-1] = 100abs((a - alat_normal) / ((a + alat_normal) / 2.))
            N_G[i-1] = 100abs((g - bulk_normal) / ((g + bulk_normal) / 2.))
            N_T[i-1] = 100abs((t - time_normal) / ((t + time_normal) / 2.))
        end
    end

    # equilibrium energy ~ k-points
    plot!(plot_EvsK_ax1, plot_first ? K : K[2:end], plot_first ? E : E[2:end], label="")
    plot!(plot_EvsK_ax2, plot_first ? K[2:end] : K[3:end], plot_first ? N_E : N_E[2:end], label="")
    # lattice parameter ~ k-points
    plot!(plot_AvsK_ax1, plot_first ? K : K[2:end], plot_first ? A : A[2:end], label="")
    plot!(plot_AvsK_ax2, plot_first ? K[2:end] : K[3:end], plot_first ? N_A : N_A[2:end], label="")
    # bulk modulus ~ k-points
    plot!(plot_GvsK_ax1, plot_first ? K : K[2:end], plot_first ? G : G[2:end], label="")
    plot!(plot_GvsK_ax2, plot_first ? K[2:end] : K[3:end], plot_first ? N_G : N_G[2:end], label="")
    # convergence time ~ iterations
    plot!(plot_TvsK_ax1, plot_first ? K : K[2:end], plot_first ? T : T[2:end], label="")
    plot!(plot_TvsK_ax2, plot_first ? K : K[2:end], plot_first ? I : I[2:end], label="")

    if ecutwfc == last(cutoff_energies)
        println(last(E)); println(last(A)); println(last(G))
    end
end


# plot everything together
p_EvsK = plot(plot_EvsK_ax1, plot_EvsK_ax2, legend, layout=@layout [grid(1, 2); a{1.0w, 0.2h}]) #, bottom_margin=(-5, :mm))
p_AvsK = plot(plot_AvsK_ax1, plot_AvsK_ax2, legend, layout=@layout [grid(1, 2); a{1.0w, 0.2h}]) #, bottom_margin=(-5, :mm))
p_GvsK = plot(plot_GvsK_ax1, plot_GvsK_ax2, legend, layout=@layout [grid(1, 2); a{1.0w, 0.2h}]) #, bottom_margin=(-5, :mm))
p_TvsK = plot(plot_TvsK_ax1, plot_TvsK_ax2, legend, layout=@layout [grid(1, 2); a{1.0w, 0.2h}]) #, bottom_margin=(-5, :mm))


# save the plots
savefig(p_EvsK, "./energy_versus_kpoint.png")
savefig(p_AvsK, "./alat_versus_kpoint.png")
savefig(p_GvsK, "./bulk_versus_kpoint.png")
savefig(p_TvsK, "./convergence_versus_kpoint.png")



# # ### EvsA
# # fig, (ax1, ax2) = plt.subplots(1, 2)
# # ax1.set(
# #     xlabel=("Lattice Parameter [" + r"$\AA$]"),
# #     ylabel=("Equilibrium Energy per Atom [" + r"$eV$]")
# # )
# # ax2.yaxis.set_ticks_position("right")
# # ax2.yaxis.set_label_position("right")
# # # ax2.set_yscale("symlog")
# # ax2.set(
# #     xlabel=("Lattice Parameter Volume [" + r"$\AA^{3}$]"),
# #     ylabel=("Equilibrium Energy per Atom [" + r"$eV$]")
# # )
# # do_full, do_bcc = False, True
# # for ecutwfc in cutoff_energies:
# #     with open(f"{filepath}/{ecutwfc}/EvsA.{ecutwfc}.{kpoint}" if do_full else f"../EvsA-{'bcc' if do_bcc else 'fcc'}-Calibration", "r") as f:
# #         lines = f.readlines()
# #         A = np.zeros((len(lines), 1))
# #         E = np.zeros((len(lines), 1))
# #         i = 0
# #         for line in lines:
# #             line_split = line.split(" ")
# #             if "\n" in line_split[1]:
# #                 line_split[1] = line_split[1][:-1]
# #             A[i] = float(line_split[0])
# #             E[i] = float(line_split[1])
# #             i += 1
# #     ax1.plot(A if plot_first else A[1:], E if plot_first else E[1:])
# # for ecutwfc in cutoff_energies:
# #     with open(f"{filepath}/{ecutwfc}/EvsA.{ecutwfc}.{kpoint}" if do_full else f"../EvsA-{'bcc' if do_bcc else 'fcc'}-Calibration", "r") as f:
# #         lines = f.readlines()
# #         V = np.zeros((len(lines), 1))
# #         E = np.zeros((len(lines), 1))
# #         i = 0
# #         for line in lines:
# #             line_split = line.split(" ")
# #             if "\n" in line_split[1]:
# #                 line_split[1] = line_split[1][:-1]
# #             V[i] = float(line_split[0])
# #             E[i] = float(line_split[1])
# #             i += 1
# #     ax2.plot(V if plot_first else V[1:], E if plot_first else E[1:])

# # # Create the legend
# # fig.legend([ax1, ax2], # The line objects
# #     labels=legend_labels, # The labels for each line
# #     loc="lower center", # Position of legend
# #     ncol=len(legend_labels),
# #     borderaxespad=0.1, # Small spacing around legend box
# #     title="Cutoff Energies" # Title for the legend
# # )
# # plt.subplots_adjust(bottom=0.225)
# # plt.savefig("{filepath}/EvsA_for_energy.png")