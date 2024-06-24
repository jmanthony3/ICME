using CSV
using DataFrames
using LinearAlgebra: norm
using NumericalMethods
using Plots; gr()

default(grid=false, framestyle=:box)

lattice_parameter = 2.781e-10 # m
slip_direction = [1., 1., 1.]
burger_vec = lattice_parameter*(norm(slip_direction) / 2.) # (np.sqrt(3)/2)
# print(burger_vec)
mu = 195.343e3 # 117e3 # MPa

# must be same examined in `rescale_commands.sh`
# if not a range, then define as blank string: TEMP = ""
TEMP = [300.] # float.(range(150, 500; step=50)) # K
# if not a range, then define as blank string: SIGMA = ""
STRAIN = [1000.] # float.(range(25, 325; step=25)) # MPa
# append!(SIGMA, float.([700, 800, 900, 1000, 1100, 1200]))

skip = 150
columns = ["timenow", "disDensity", "Stress", "Strain", "S1", "S2", "S3", "S23", "S31", "S12", "p1", "p2", "p3", "p23", "p31",  "p12", "jogs", "junctions", "CrossSlip"]

do_monitor = false

H0 = range(0.4, 0.42; step=0.0001) .* 1e6 # Pa
# H0 = range(0, 101; step=1) .* 1e6 # Pa
C = [1000.] # float.(range(73.5, 73.61; step=0.01)) # 1/s
# C = float.(range(0, 101; step=1)) # 1/s
Kappa(ks, k0, h0, c, t) = ks - ((ks - k0) * exp(-h0 / (ks - k0) * c * t))

# old_stdout = sys.stdout
# sys.stdout = open(os.devnull, "w")
# try:
fig, (ax_kappa, ax_alpha) = plt.subplots(1, 2)
# fig, ax_kappa = plt.subplots(1, 1)
ax_kappa = plot(
    xlabel="Normalized True Strain (\$\\epsilon\$)",
    # xscale=:log10,
    ylabel="True Stress (\$\\sigma\$) [\$MPa\$]",
    # yscale=:log10,
    title="Stress-Strain @ 300 K",
)
ax_alpha = plot(
    # xlim=(0, 0.0002),
    xlabel="Time (\$t\$) [\$\\mu-s\$]",
    # xscale=:log10,
    ylabel="Junction Strength (\$\\alpha\$)",
    yguide_position=:right,
    # yscale=:log10,
    ymirror=true,
    title="Average Value of\nJunction Strength @ 300 K",
)
if do_monitor
    while do_monitor
        # # MDDP_BCC_HW2/MDDP_Zip/Windows/
        # open(f"../Homework 2/MDDP_BCC_HW2/MDDP_Zip/Windows/DDtimeResults.out", "r") do f
        # # with open(f"./input", "r") as f:
        # end
        data = CSV.read("../Homework 2/MDDP_BCC_HW2/MDDP_Zip/Windows/DDtimeResults.out", DataFrame; header=columns, skipto=4)
        # data.plot(
        #     kind="scatter",
        #     x="Strain",
        #     y="Stress",
        #     color="blue",
        #     ax_kappa=plt.gca())
        scatter!(ax_kappa,
            data[!, "Strain"][begin:skip:end], # [:int(1e1)],
            data[!, "Stress"][begin:skip:end] ./ 1e6, # [:int(1e1)]/1e6
            ) # label=f"{strain} " + r"$\frac{1}{s}$")
        # ax_kappa.legend()
        savefig("./ENGR851_JobyAnthonyIII_Homework3/stress_strain-monitor.png")
        sleep(30)
    end
else
    for temp in TEMP
        for strain in STRAIN
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
            # open("../Homework 2/MDDP_BCC_HW2/MDDP_Zip/Windows/$temp/$strain/DDtimeResults.out", "r") do f
            # # with open(f"./input", "r") as f:
            # end
            data = CSV.read("../Homework 2/MDDP_BCC_HW2/MDDP_Zip/Windows/$temp/$strain/DDtimeResults.out", DataFrame; header=columns, skipto=4)
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
            time = data[!, "timenow"][begin:skip:end]
            eps = data[!, "Strain"][begin:skip:end] # /data["Strain"].values[-1],#[:int(1e1)],
            sigma = data[!, "Stress"][begin:skip:end] ./ 1e6 # [:int(1e1)]/1e6
            # rho_fp = np.array([rho_f[i] for i in range(0, len(data), skip)])
            # rho_f = rho_fp
            rho_f = data[!, "disDensity"][begin:skip:end]
            # data.plot(
            #     kind="scatter",
            #     x="Strain",
            #     y="Stress",
            #     color="blue",
            #     ax_kappa=plt.gca())
            scatter!(ax_kappa,
                eps, sigma,
                label="\$\\epsilon_{f}\$ = $(last(eps)) @ \$\\dot{\\epsilon}\$ = $strain \$\\frac{1}{s}\$")
            ERRORS = Dict{String, Any}()
            for j in 0:1:15
                ERRORS["$j"] = Dict{String, Any}()
                for h0 in H0
                    ERRORS["$j"]["$(h0 / 1e6)"] = Dict{String, Any}()
                    for c in C
                        kappa, error = zeros(length(time)), zeros(length(time))
                        for (i, t) in enumerate(time)
                            kappa[i] = Kappa(average(sigma[end-j:end]), first(sigma), h0, c, t)
                            if i > 1
                                error[i] = abs((kappa[i] - sigma[i]) / sigma[i])
                            end
                        end
                        ERRORS["$j"]["$(h0 / 1e6)"]["$c"] = sum(error)
                    end
                end
            end
            min_error = ["", "", "", 1e6]
            for j in 0:1:15, h0 in ERRORS["$j"], c in ERRORS["$j"][h0]
                if ERRORS["$j"][h0][c] < last(min_error)
                    min_error = [j, h0, c, ERRORS["$j"][h0][c]]
                end
            end
            print(min_error)
            j, h0, c, error = min_error
            kappa, alpha = np.zeros(length(time)), zeros(length(time))
            for (i, t) in enumerate(time)
                kappa[i] = Kappa(average(sigma[end-j:end]), first(sigma), parse(Float64, h0) * 1e6, parse(Float64, c), t)
                alpha[i] = kappa[i] / (mu * burger_vec * sqrt(rho_f[i]))
            end
            print(first(kappa))
            # foo = nm.MultiVariableIteration(Kappa, np.array([h0, C]), kappa)
            scatter!(ax_kappa,
                eps, kappa,
                label="Voce Equation, \$\\kappa = \\kappa_{s} - (\\kappa_{s} - \\kappa_{0})\\exp(-\\frac{h_{0}}{\\kappa_{s} - \\kappa_{0}}Ct)\$",
                annotation=(0.65, 0.4,
                    "\$\\kappa_{s}\$ = $(round(average(sigma[end-j:end]); digits=4)) MPa\n\$\\kappa_{0}\$ = $(round(first(sigma); digits=4)) MPa\n\$C\$ = $c \$\\frac{1}{s}\$\n\$h_{0}\$ = $(round(h0; digits=4)) MPa\nE = $(round(error; digits=4)) \$\\frac{MPa}{MPa}\$"))
            scatter!(ax_alpha, time .* 1e3, alpha)
        end
        savefig("./ENGR851_JobyAnthonyIII_Homework3/stress_strain-$temp.png")
        # finally:
        #     sys.stdout.close()
        #     sys.stdout = old_stdout
        #     print(DRAG)
    end
end