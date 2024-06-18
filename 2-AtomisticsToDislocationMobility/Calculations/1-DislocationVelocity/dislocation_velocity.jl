using CSV
using DataFrames
using LinearAlgebra: norm
using NumericalMethods
using Plots

average(array::AbstractVector)::Real = sum(array)/length(array)

filepath = pwd()

print("")

dt, width = 0.500, 40 * 5   # lattice distance * unit distance (approx)
offset      = 40.0
L_buff      = 0.0 + offset
R_buff      = width - offset
midpoint    = width / 2

lattice_parameter = 2.7810811e-10 # m
slip_direction = [1, 1, 1]
burger_vec = lattice_parameter * (norm(slip_direction) / 2.)
println("Magnitude of Burger's vector, ||b|| = $burger_vec m\n")

# must be same examined in `rescale_commands.sh`
TEMP = [300] # range(150, 500; step=50) # K
SIGMA = range(25, 300; step=25) # MPa
# for s in [700, 800, 900, 1000, 1100, 1200]
#     push!(SIGMA, s)
# end

skip = 1
columns = ["Particle Identifier", "X", "Y", "Z", "Centrosymmetry", "R", "G", "B", "nan", "Time"]

DRAG = Dict{String, NTuple{2, Float64}}()

old_stdout = sys.stdout
sys.stdout = open(os.devnull, "w")
try
    for temp in TEMP
        plot_temp = plot(
            xlabel="Stress (" + r"$\tau$) [$MPa$]",
            ylabel="Velocity (" + r"$v_{x}$) [$\frac{m}{s}$]",
            title="Dislocation Velocity for Stress")
        open("dislocation_velocity-$temp", "w") do f
            plottemp_sigma = plot(
                xlabel="Time (t) [" + r"$ps$]",
                ylabel="Displacement (x) [" + r"$nm$]",
                title="Dislocation Position for Time",
            )
            println(f, "Stress Velocity")
            for sigma in SIGMA
                @printf("Looking at %3d K, %5.3f MPa...\r", temp, sigma)
                x_ave, lap_count = 0, 0
                dir_data = @sprintf("%s/PositionFrameData/%d/%.3f", filepath, temp, sigma)

                velocities, offsets = [], []
                time, position = [], []
                for (i, frame) in enumerate(range(1, length([_ for _ in walkdir(dir_data)...][3]), skip))
                    data = CSV.read("$data_dir/lammps.$frame.xyz", DataFrame)
                    try
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
                    catch exc
                        if isa(exc, DomainError)
                            continue
                        end
                    end
                end
                # time        = time  # ps
                position  ./= 10.   # nm
                scatter!(plottemp_sigma, time, position) # , label="Raw")

                try
                    position_expr = linearleastsquares(time, position, 1)[1]
                    other_pos = position_expr.(time)
                    # plot!(plottemp_sigma, time, other_pos, "k", label="Homogeneous Fit")

                    if String(position_expr(:t))[end] != "t"
                        velocity = float(str(position_expr(sp.Symbol("t"))).split(" ")[0].split("*")[0])
                        offset = float(str(position_expr(sp.Symbol("t"))).split(" ")[-1])
                    else
                        velocity = float(str(position_expr(sp.Symbol("t"))).split(" ")[-1].split("*")[0])
                        offset = float(str(position_expr(sp.Symbol("t"))).split(" ")[0])
                    end
                    push!(velocities, velocity)
                    push!(offsets, offset)
                catch exc
                    if isa(exc, BoundsError)
                        continue
                    end
                end

                velocity = average(velocities)
                @printf(f, "%.9f %.9f\n", sigma, 1e3*velocity) # MPa m/s
                offset = average(offsets)

                pft_string = try
                    @sprintf("%.3f*t + %.3f", velocity, offset)
                catch exc
                    if isa(exc, MethodError)
                        "Could not round..."
                    end
                end
                time = range(0, 100, dt*skip)
                other_pos = velocity*time .+ offset
                plot!(plottemp_sigma, time, other_pos, label="$sigma MPa | $pft_string")
                @printf(f"Closing %d K, %.3f MPa...\r", temp, sigma)
            end
            plottemp_sigma.text(
                x=plottemp_sigma.get_xlim()[1]*0.175,
                y=plottemp_sigma.get_ylim()[1]*-0.025,
                s=r"The $v$'s in legend equations are dislocation velocities."
            )
            plottemp_sigma.legend(
                title=r"$x(t) = v[\frac{nm}{ps}]*t[ps] + [nm]$"
            )
            plottemp_sigma.savefig(f"position_for_time-{temp}.png")
        end
        data = CSV.read("dislocation_velocity-$temp", DataFrame; header=["Stress", "Velocity"])
        drag = burger_vec / parse(Float64, split(String(linearleastsquares(data[!, "Stress"][begin:11], data[!, "Velocity"][begin:11], 1)[0](:t)))[1], '*')*1e6 # Pa-s
        DRAG["$temp"] = (drag, 1/drag)
        scatter!(plot_temp, data[!, "Stress"], data[!, "Velocity"], label="$temp K \$\\rightarrow\$ B = $drag Pa-s")
        savefig(plot_temp, "velocity_for_stress.png")
    end
finally
    sys.stdout.close()
    sys.stdout = old_stdout

    # print the names of the columns
    @printf("%15s | %30s | %35s\n", "Temperature [K]", "Drag Coefficient (B) [Pa-s]", "Dislocation Mobility (Ms) [1/Pa-s]")
    println('-'^16 * '|' * '-'^32 * '|' * '-'^36)

    # print each data item
    for (key, entry) in sort(collect(pairs(DRAG)), by=x->parse(Int64, x[1]), rev=false)
        coeff, mob = entry
        @printf("%15s | %30s | %35s\n", key, "$coeff", "$mob")
    end
end