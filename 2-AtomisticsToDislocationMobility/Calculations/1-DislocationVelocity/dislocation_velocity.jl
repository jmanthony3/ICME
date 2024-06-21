using CSV
using DataFrames
using LinearAlgebra: norm
using NumericalMethods
using Plots; gr()
using Printf
using Symbolics: Num, @variables, simplify, value, degree

default(grid=false, framestyle=:box)

average(array::AbstractVector)::Real = sum(array)/length(array)

filepath = pwd()

println("")

dt, width = 0.500, 40. * 5. # lattice distance * unit distance (approx)
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
SIGMA = float.(range(25, 300; step=25)) # MPa
# append!(SIGMA, float.([700, 800, 900, 1000, 1100, 1200]))

@variables t
skip = 1
columns = ["Particle Identifier", "X", "Y", "Z", "Centrosymmetry", "R", "G", "B", "nan", "Time"]
DRAG = Dict{String, NTuple{2, Float64}}()
for temp in TEMP
    plot_temp = plot(
        xlabel="Stress (\$\\tau\$) [\$MPa\$]",
        ylabel="Velocity (\$v_{x}\$) [\$\\frac{m}{s}\$]",
        title="Dislocation Velocity for Stress")
    open(@sprintf("dislocation_velocity-%.3f.dat", temp), "w") do f
        plottemp_sigma = plot(
            xlabel="Time (\$t\$) [\$ps\$]",
            ylabel="Displacement (\$x\$) [\$nm\$]",
            legendtitle="\$x(t) = v[\\frac{nm}{ps}]*t[ps] + [nm]\$",
            title="Dislocation Position for Time",
            annotation=((0.175, -0.025),
                "The \$v\$'s in legend equations are dislocation velocities.")
        )
        println(f, "Stress Velocity")
        for sigma in SIGMA
            @printf("Looking at %.3f K, %.3f MPa...\r", temp, sigma)
            x_ave, lap_count = 0., 0
            dir_data = @sprintf("%s/PositionFrameData/%.3f/%.3f", filepath, temp, sigma)

            time, position, velocities, offsets = Float64[], Float64[], Float64[], Float64[]
            for (i, frame) in enumerate(range(0, length([_ for _ in walkdir(dir_data)...][3])-1; step=skip))
                data = CSV.read("$dir_data/lammps.$frame.xyz", DataFrame; header=columns, skipto=3)
                if nrow(data) > 0
                    push!(time, frame*dt)
                    x_store = data[!, "X"]
                    # check to see if the dsl cluster is on wraping both sides
                    if minimum(x_store) < L_buff && R_buff < maximum(x_store)
                        for (j, x) in enumerate(x_store)
                            # values on left side (wrapped around)
                            if x < midpoint
                                x_store[j] = x + width
                            end
                        end
                    end
                    x_ave_prev = x_ave
                    x_ave = average(x_store) + width*lap_count
                    # check if a lap just happened
                    if x_ave < (x_ave_prev - midpoint)
                        lap_count += 1
                        x_ave += width
                    end
                    push!(position, x_ave)
                end
            end
            # time        = time  # ps
            position  ./= 10.   # nm
            scatter!(plottemp_sigma, time[begin:5:end], position[begin:5:end], label="")

            pft_f, pft_rmse = linearleastsquares(time, position, 1)
            # other_pos = pft_f.(time)
            # plot!(plottemp_sigma, time, other_pos, "k", label="Homogeneous Fit")

            pft_g_value = value(simplify(pft_f(t); expand=true))
            v           = if hasproperty(pft_g_value, :dict)
                coeffs_map = pft_g_value.dict
                collect(values(coeffs_map))[sortperm(degree.(keys(coeffs_map)))][1]
            end
            pft_offset      = pft_f(0.)
            push!(velocities, v); push!(offsets, pft_offset)

            velocity = average(velocities)
            @printf(f, "%.9f %.9f\n", sigma, 1e3*velocity) # MPa m/s
            pft_offset = average(offsets)

            pft_string = @sprintf("%.3f*t + %.3f ± %.3f", velocity, pft_offset, pft_rmse)
            time = range(0., 100.; step=dt*skip)
            other_pos = velocity .* time .+ pft_offset
            plot!(plottemp_sigma, time, other_pos, label="$sigma MPa | $pft_string")
            print("\e[2K"); @printf("Closing %.3f K, %.3f MPa...\r", temp, sigma)
        end
        savefig(plottemp_sigma, @sprintf("position_for_time-%.3f.png", temp))
    end
    data    = CSV.read(@sprintf("dislocation_velocity-%.3f.dat", temp), DataFrame; header=["Stress", "Velocity"], skipto=2, delim=' ', types=Float64)
    f       = linearleastsquares(data[!, "Stress"][begin:12], data[!, "Velocity"][begin:12], 1)[1]
    g_value = value(simplify(f(t); expand=true))
    drag    = burger_vec / if hasproperty(g_value, :dict)
        coeffs_map = g_value.dict
        collect(values(coeffs_map))[sortperm(degree.(keys(coeffs_map)))][1]
    end * 1e6 # Pa-s
    DRAG["$temp"] = (drag, 1. / drag)
    scatter!(plot_temp, data[!, "Stress"], data[!, "Velocity"], label="$temp K \$\\rightarrow\$ B = $drag Pa-s")
    savefig(plot_temp, "velocity_for_stress.png")
end

# print the names of the columns
print("\e[2K")
@printf("%15s | %30s | %35s\n", "Temperature [K]", "Drag Coefficient (B) [Pa-s]", "Dislocation Mobility (Ms) [1/Pa-s]")
println('-'^16 * '|' * '-'^32 * '|' * '-'^36)

# print each data item
for (key, entry) in sort(collect(pairs(DRAG)), by=x->parse(Int64, x[1]), rev=false)
    coeff, mob = entry
    @printf("%15s | %30s | %35s\n", key, "$coeff", "$mob")
end