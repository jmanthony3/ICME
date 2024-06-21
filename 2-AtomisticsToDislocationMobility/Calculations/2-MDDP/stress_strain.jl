using CSV
using DataFrames
using LinearAlgebra: norm
using NumericalMethods
using Plots; gr()

default(grid=false, framestyle=:box)

filepath = pwd()

print("")

dt, width = 500., 40. * 5.   # lattice distance * unit distance (approx)
offset      = 40.0
L_buff      = 0.0 + offset
R_buff      = width - offset
midpoint    = width / 2.

lattice_parameter = 2.7810811e-10 # m
slip_direction = [1., 1., 1.]
burger_vec = lattice_parameter * (norm(slip_direction) / 2.)
println("Magnitude of Burger's vector, ||b|| = $burger_vec m\n")

# must be same examined in `rescale_commands.sh`
TEMP = [300.] # float.(range(150, 500; step=50)) # K
STRAIN = [3.] # float.(range(1, 6; step=1)) # 1/s
# append!(STRAIN, float.([-3, -2, -1]))
STRAIN .^= 10.

skip = 150
columns = ["timenow", "disDensity", "Stress", "Strain", "S1", "S2", "S3", "S23", "S31", "S12", "p1", "p2", "p3", "p23", "p31",  "p12", "jogs", "junctions", "CrossSlip"]

do_monitor = true

if do_monitor
    data = CSV.read("$filepath/DDtimeResults.out", DataFrame; header=columns, skipto=4, delim=' ', ignorerepeated=true, types=Float64)
    ax = scatter(
        data[!, "Strain"][begin:skip:end],
        data[!, "Stress"][begin:skip:end] ./ 1e6, # MPa
        label="", # "$strain \$\\frac{1}{s}\$",
        # xscale=:log10,
        xlabel="True Strain (\$\\epsilon\$)",
        # yscale=:log10,
        ylabel="True Stress (\$\\sigma\$) [\$MPa\$]",
        title="Stress-Strain Curve",
    )
    savefig(ax, "stress_strain-monitor.png")
else
    for temp in TEMP
        local ax = plot(
            # xscale=:log10,
            xlabel="True Strain (\$\\epsilon\$)",
            # yscale:log10,
            ylabel="True Stress (\$\\sigma\$) [\$MPa\$]",
            title="Stress-Strain @ $temp [\$K\$]",)
        for strain in STRAIN
            local data = CSV.read("$filepath/$temp/$strain/DDtimeResults.out", DataFrame; header=columns, skipto=5)
            scatter!(ax,
                data["Strain"][1:skip:end], # ./ last(data[!, "Strain"]),
                data["Stress"][1:skip:end] ./ 1e6, # MPa
                label="\$\\epsilon_{f}\$ = $(last(data[!, "Strain"])) @ \$\\dot{\\epsilon}\$ = $strain \$\\frac{1}{s}\$")
        end
        savefig(ax, "$filepath/stress_strain-$temp.png")
    end
end