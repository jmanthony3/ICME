using CSV
using DataFrames
using LinearAlgebra: norm
using NumericalMethods
using Plots

filepath = pwd()

print("")

dt, width = 500, 40 * 5   # lattice distance * unit distance (approx)
offset      = 40.0
L_buff      = 0.0 + offset
R_buff      = width - offset
midpoint    = width / 2

lattice_parameter = 2.7810811e-10 # m
slip_direction = [1, 1, 1]
burger_vec = lattice_parameter * (norm(slip_direction) / 2.)
println("Magnitude of Burger's vector, ||b|| = $burger_vec m\n")

# must be same examined in `rescale_commands.sh`
TEMP = [300] # range(150, 500, 50) # K
STRAIN = [3] # range(1, 6, 1) # 1/s
# for s in [-3, -2, -1]:
#     push!(STRAIN, s)
# end
STRAIN .^= 10.

skip = 150
columns = ["timenow", "disDensity", "Stress", "Strain", "S1", "S2", "S3", "S23", "S31", "S12", "p1", "p2", "p3", "p23", "p31",  "p12", "jogs", "junctions", "CrossSlip"]

do_monitor = true

if do_monitor
    data = CSV.read("$filepath/DDtimeResults.out", DataFrame; header=columns, skipto=5)
    ax = scatter(
        data[!, "Strain"][1:skip:end],
        data[!, "Stress"][1:skip:end] ./ 1e6, # MPa
        # label="$strain \$\\frac{1}{s}\$")
        # xscale=:log10,
        xlabel="True Strain (\$\\epsilon\$)",
        # yscale=:log10,
        ylabel="True Stress (\$\\sigma\$) [\$MPa\$]",
        title="Stress-Strain Curve",
    )
    savefig(ax, "stress_strain-monitor.png")
else
    for temp in TEMP
        ax = plot(
            # xscale=:log10,
            xlabel="True Strain (\$\\epsilon\$)",
            # yscale:log10,
            ylabel="True Stress (\$\\sigma\$) [\$MPa\$]",
            title="Stress-Strain @ $temp [\$K\$]",)
        for strain in STRAIN
            data = CSV.read("$filepath/$temp/$strain/DDtimeResults.out", DataFrame; header=columns, skipto=5)
            scatter!(ax,
                data["Strain"][1:skip:end], # ./ last(data[!, "Strain"]),
                data["Stress"][1:skip:end] ./ 1e6, # MPa
                label="\$\\epsilon_{f}\$ = $(last(data[!, "Strain"])) @ \$\\dot{\\epsilon}\$ = $strain \$\\frac{1}{s}\$")
        end
        savefig(ax, "$filepath/stress_strain-$temp.png")
    end
end