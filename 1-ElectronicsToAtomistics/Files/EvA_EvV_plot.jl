#####################################################################################
# EvA_EvV_plot.jl                                                                   #
#                                                                                   #
# The purpose of this code is to create plots for any set of EvsA and EvsV data.    #
# It should be run from the same directory that contains EvsA and EvsV data files.  #
# This code can print out figures both to the screen and to external files          #
# depending on user settings.                                                       #
#                                                                                   #
# Data files are sorted in 2 columns by the following format:                       #
# EvsA | Lattice Parameter [Angstrom]       energy[eV]                              #
# EvsV | Cell Volume [Angstrom^3]           energy[eV]                              #
#                                                                                   #
#                                                                                   #
# CHANGELOG                                                                         #
# 2024-06-14: Converted from Python3 to Julia and modified this title block.        #
#                                                                                   #
# 2022-XX-XX: Original import from Dr. Heechen Cho.                                 #
#####################################################################################

# import packages
using CSV
using DataFrames
using LaTeXStrings
using Plots; gr()



#######################################################
###                  User Settings                  ###
#######################################################

default(grid=false, framestyle=:box)

# display or print out
write_to_screen = false
printout        = true

energy_offset = 4479.41619611 # Set by running simulation with very large lattice parameter

# color options for plots (https://xkcd.com/color/rgb/)
color1 = "#b66325"
color2 = "#8b2e16"
color3 = "#00035b"

marker_shape    = :circle   # type of marker used for plotting google matplotlib for more options
marker_size     = 4         # markersize for plots, larger is bigger
ax_label_size   = 18        # font size for axis labels
title_font_size = 24        # font size for table title

EvsA_fn         = "Name_of_EvA.png"
EvsV_fn         = "Name_of_EvV.png"
EvsA_EvsV_fn    = "Name_of_Combined.png"

# path to data files
EvsA_df = CSV.read("EvsA", DataFrame; header=false, delim=' ', types=Float64)
EvsV_df = CSV.read("EvsV", DataFrame; header=false, delim=' ', types=Float64)


#######################################################
###            End of  User Settings                ###
#######################################################


# initialize arrays
latt_param  = EvsA_df[!, 1] # array of Lattice Parameter in Angstrom
volume      = EvsV_df[!, 1] # array of Volume in Angstrom^3
energy      = EvsA_df[!, 2] # array of Energy in eV
energy    .+= energy_offset # shifting energy values

# creating plots
p               = scatter(              # basic structure/layout of plots
    markershape=marker_shape,
    markersize=marker_size,
    grid=true,
    gridstyle=:dashdot,
    guidelinecolor=:black,
    guidelinewidth=1,
    guidefontsize=ax_label_size,
    titlefontsize=title_font_size,
    framestyle=:box)
EvsA_plot       = scatter!(deepcopy(p), # bond energy vs lattice parameter
    latt_param, energy,
    label=:none,
    markercolor=color1,
    xlabel="Lattice Parameter (\$\\AA\$)",
    ylabel="Energy (\$eV\$)",
    title="Energy vs Lattice Parameter\nCurve for FCC copper\n")
EvsV_plot       = scatter!(deepcopy(p), # bond energy vs cell volume
    volume, energy,
    label=:none,
    markercolor=color1,
    xlabel="Volume (\$\\AA^{3}\$)",
    ylabel="Energy (\$eV\$)",
    title="Energy vs Volume\nCurve for FCC copper\n")
EvsA_EvsV_plot  = plot(EvsA_plot, plot(EvsV_plot, ylabel=""),
    markercolor=color2,
    title="",
    layout=(1, 2),
    link=:y,
    plot_title="Energy vs Lattice Parameter and Energy vs Cell Volume\n",
    plot_titlefontsize=ax_label_size+4,
    size=(1000, 600)
)

# print to screen
if write_to_screen
    display(EvsA_plot); display(EvsV_plot); display(EvsA_EvsV_plot)
end

# print out plots
if printout
    savefig(EvsA_plot, EvsA_fn)
    savefig(EvsV_plot, EvsV_fn)
    savefig(EvsA_EvsV_plot, EvsA_EvsV_fn)
end
