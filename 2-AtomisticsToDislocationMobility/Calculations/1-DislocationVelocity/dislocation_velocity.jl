using CSV
using DataFrames
using LinearAlgebra: norm
using Plots

filepath = pwd()

print("")

dt, width = 0.500, 40 * 5   # lattice distance * unit distance (approx)
offset      = 40.0
L_buff      = 0.0 + offset
R_buff      = width - offset
midpoint    = width / 2

lattice_parameter = 2.781e-10 # m
slip_direction = [1, 1, 1]
burger_vec = lattice_parameter * (norm(slip_direction) / 2.)
println("Magnitude of Burger's vector, ||b|| = $burger_vec m\n")

# must be same examined in `rescale_commands.sh`
TEMP = [300] # range(150, 550; step=50) # K
SIGMA = range(25, 325; step=25) # MPa
# for s in [700, 800, 900, 1000, 1100, 1200]
#     push!(SIGMA, s)
# end

skip = 1
columns = ["Particle Identifier", "X", "Y", "Z", "Centrosymmetry", "R", "G", "B", "nan", "Time"]

DRAG = {}

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
                    open("$data_dir/lammps.$frame.xyz") do g
                        try
                            py_array = np.genfromtxt(g.readlines()[2:])
                            if len(py_array) > 0
                                data = pd.DataFrame(py_array, columns=columns)
                                time.append(frame*dt)
                                x_store = data[["X"]].values
                                # check to see if the dsl cluster is on wraping both sides
                                if min(x_store) < L_buff and R_buff < max(x_store)
                                    for j in range(len(x_store))
                                        #values on left side (wrapped around)
                                        if x_store[j] < midpoint
                                            x_store[j] = x_store[j] + width
                                        end
                                    end
                                end
                                x_ave_prev = x_ave
                                x_ave = np.average(x_store) + width*lap_count
                                # check if a lap just happened
                                if x_ave < (x_ave_prev - midpoint)
                                    lap_count += 1
                                    x_ave += width
                                end
                                position.append(x_ave)
                            end
                        except ValueError: pass
                        except UserWarning: pass
                        end
                    end
                    i += 1
                end
                time, position = np.array(time), np.array(position)/10 # ps, nm
                plottemp_sigma.scatter(time, position)#, label="Raw")

                try
                    position_expr = nm.least_squares.linear(time, position, 1)[0]
                    other_pos = np.array([position_expr(t) for t in time])
                    # plottemp_sigma.plot(time, other_pos, "k", label="Homogeneous Fit")

                    if "t" != str(position_expr(sp.Symbol("t")))[-1]
                        velocity = float(str(position_expr(sp.Symbol("t"))).split(" ")[0].split("*")[0])
                        offset = float(str(position_expr(sp.Symbol("t"))).split(" ")[-1])
                    else
                        velocity = float(str(position_expr(sp.Symbol("t"))).split(" ")[-1].split("*")[0])
                        offset = float(str(position_expr(sp.Symbol("t"))).split(" ")[0])
                    end
                    velocities.append(velocity)
                    offsets.append(offset)
                except np.linalg.LinAlgError: pass
                end

                velocity = np.average(velocities)
                f.write(f"{sigma:.9f} {velocity*1e3:.9f}\n") # MPa m/s
                offset = np.average(offsets)

                try
                    pft_string = f"{round(velocity, 3)}*t + {round(offset, 3)}"
                except ValueError
                    pft_string = f"Could not round..."
                end
                time = np.arange(0, 100, dt*skip)
                other_pos = velocity*time + offset
                plottemp_sigma.plot(time, other_pos, label=f"{int(sigma)}".zfill(4) + f" MPa | {pft_string}".ljust(15))
                print(f"Closing {temp} K, {sigma:.3f} MPa...", end="\r")
            end
            plottemp_sigma.text(
                x=plottemp_sigma.get_xlim()[1]*0.175,
                y=plottemp_sigma.get_ylim()[1]*-0.025,
                s=r"The $v$'s in legend equations are dislocation velocities."
            )
            plottemp_sigma.legend(
                title=r"$x(t) = v[\frac{nm}{ps}]*t[ps] + [nm]$"
            )
            plottemp_sigma.savefig(f"position_for_time-{temp}.svg")
        end
        open(f"dislocation_velocity-{temp}", "r") do f
            data = pd.DataFrame(np.genfromtxt(f.readlines())[1:], columns=["Stress", "Velocity"])
            drag = burger_vec/float(str(nm.least_squares.linear(data['Stress'][:12], data['Velocity'][:12], 1)[0](sp.Symbol('t'))).split('*')[0])*1e6 # Pa-s
            DRAG[f"{temp}"] = (drag, 1/drag)
            plot_temp.scatter(data["Stress"], data["Velocity"], label=f"{temp} K" + r"$\rightarrow$ B = " + f"{drag:3e}" + r" Pa-s")
        end
        plot_temp.legend()
        plot_temp.savefig(f"velocity_for_stress.svg")
    end
finally
    sys.stdout.close()
    sys.stdout = old_stdout

    # Print the names of the columns.
    print("{:<15} | {:<30} | {:<35}".format("Temperature [K]", "Drag Coefficient (B) [Pa-s]", "Dislocation Mobility (Ms) [1/Pa-s]"))
    print("-"*16 + "|" + "-"*32 + "|" + "-"*36)

    # print each data item.
    for key, value in DRAG.items()
        coeff, mob = value
        print("{:<15} | {:<30} | {:<35}".format(key, coeff, mob))
    end
end