#####################################################################################
# OutputFileCreator.jl                                                              #
#                                                                                   #
# CHANGELOG                                                                         #
# 2024-06-17: Converted from Python2 to Julia and added this title block.           #
#                                                                                   #
# 2022-XX-XX: Original import from Sam Scott                                        #
#####################################################################################

using Dates
using LinearAlgebra: norm


#
### EDIT THIS PARAMETER FOR CODE IN USE
# 'qe' - Quantum Espresso
# 'vasp' - VASP
dft = "qe"

# Number of processors
num_proc = 6

##### MATERIAL SETTINGS - USER EDIT REQUIRED #####

# vvv These parameters are for QE only, VASP will use the available POTCAR file vvv
# Directory containing your pseudopotentials
pp_dir = "."
# Element name
el = "Cu"
# Potential file name
potential = "Cu.UPF"
# Element weight
el_weight = 63.546
# ^^^


##### DFT PARAMETERS - USER EDIT IF DESIRED #####
# Energy cutoff value (eV) (qe only)
energy_cutoff = 30.0*13.6057
# The below energy cutoff setting works for most materials.
# Check your pseudopotential file to see the suggested value for your material
energy_cutoff_rho = energy_cutoff * 4.0
# K-points specification
kpoints = 8
# Smearing parameter (QE only)
smear = 0.06
# Relax the ions in the z direction? (vasp only right now)
relax = false


##### STACKING FAULT PARAMETERS - USER EDIT IF NEEDED #####
# Number or location of simulated points. It MUST start with 0.
#fault_points = [0.0, 1./4., 1./2.]
# fault_points = np.linspace(0,1,16)
fault_points = range(0, 16; step=1)
# Size of vacuum in Ang.
vacuum = 20.0
# Number of stacking layers - CAREFUL, DFT scales very poorly with more atoms
stacking_layers = 10


##### Unit conversions #####
ry_to_ev = 13.6056849587
au_to_ang = 0.52917721092

##### Below this line should not require user edits #####

# function that creates the geometry and writes the input file
function create_stackingfault(structure, lp, layers, slip, system=nothing)
    # structure : reference structure
    # lp : lattice parameter\
    # layers : number of layers of atoms to generate
    # slip : normalized displacement in slip direction
    # system : the type of stacking fault to create

    base_atoms  = zeros(3)
    if lowercase(structure) == "fcc"
        basis       = [
            [ √(2.)/4.  √(6.)/4.        0.]
            [-√(2.)/4.  √(6.)/4.        0.]
            [       0.  √(6.)/6.   √(3)/3.]]
        slip_vector = if isnothing(system) || system == "partial" || system == "shockley"
            # Assume the partial dislocation is desired (most common)
            [0., √(6.)/2., 0.] # Partial dislocation line
        elseif system == "full" || system == "burgers"
            basis[1,:] # burgers vector direction
        end
    elseif lowercase(structure) == "bcc"
        basis       = [
            [ (√(2.)/4. + 0.5)  -(0.5 - √(2.)/4.)         0.]
            [-(0.5 - √(2.)/4.)   (√(2.)/4. + 0.5)         0.]
            [        √(2.)/4.     √(2.)/4.          √(2.)/2.]]
        slip_vector = if isnothing(system) || system == "partial"
            # Assume the partial dislocation is desired (most common)
            [√(2.)/2., √(2.)/2., 0.] # Partial dislocation line?
        elseif system == "longpartial"
            # or?
            [0.5, -0.5, 0.]
        elseif system == "full" || system == "burgers"
            basis[1,:] # burgers vector directionn
        end
    end

    # create atoms
    threshold   = last(cross([basis[1,1:2]..., 0.], [basis[2,1:2]..., 0.]))
    f1(x)       = -basis[1,2]*x[1] + basis[1,1]*x[2]
    f2(x)       = -basis[2,2]*x[1] + basis[2,1]*x[2]
    f3(x)       = -basis[2,2]*x[1] + basis[2,1]*x[2]
    f4(x)       = -basis[1,2]*x[1] + basis[1,1]*x[2]
    b1, b2      = basis[1,:], basis[2,:]
    use_other   = false
    if threshold < 0.
        threshold   = last(cross([basis[2,1:2]..., 0.], [basis[1,1:2]..., 0.]))
        b1, b2      = basis[2,:], basis[1,:]
        use_other   = true
    end
    d1(x) = use_other ? f3(x) : f1(x)
    d2(x) = use_other ? f4(x) : f2(x)

    old_atom = zeros(3)
    atoms = [old_atom]; for i in range(2, layers)
        new_atom = old_atom + basis[3,:]

        # check if atom is within cell boundaries
        in_cellx, in_celly = false, false
        while !in_cellx || !in_celly
            if d1(new_atom) < 0.
                new_atom   += b2
            elseif d1(new_atom) > threshold
                new_atom   -= b2
            else
                in_cellx    = true
            end
            if d2(new_atom) > 0.
                new_atom   += b1
            elseif d2(new_atom) < -threshold
                new_atom   -= b1
            else
                in_celly    = true
            end
        end

        old_atom = deepcopy(new_atom)

        if i >= Int(div(layers/2, 1))
            new_atom += slip .* slip_vector

            # check if atom is within cell boundaries
            in_cellx, in_celly = false, false
            while !in_cellx || !in_celly
                if d1(new_atom) < 0.
                    new_atom   += b2
                elseif d1(new_atom) > threshold
                    new_atom   -= b2
                else
                    in_cellx    = true
                end
                if d2(new_atom) > 0.
                    new_atom   += b1
                elseif d2(new_atom) < -threshold
                    new_atom   -= b1
                else
                    in_celly    = true
                end
            end
        end

        push!(atoms, new_atom)
    end

    if dft == "vasp"
        write_vasp_inputs(lp, basis, layers, atoms)
    elseif dft == "qe"
        write_qe_inputs(lp, basis, layers, atoms)
    end

    # returns the stacking fault area
    return [norm(cross(basis[1,:], basis[2,:])), norm(slip_vector), length(atoms)]
end # end create_stackingfault

function write_qe_inputs(lp, basis, layers, atoms)
    open("gsfe.in", "w") do f # write the input file
        # control section
        println(f, " &control")
        println(f, "\tprefix=''")
        println(f, "\toutdir='temp'")
        println(f, "\tpseudo_dir = '$pp_dir',")
        println(f, " /")
        # system section
        println(f, " &system")
        println(f, "\tibrav= 0, nat= $layers, ntyp= 1,")
        println(f, "\tcelldm(1) =$(lp/au_to_ang), ")  # convert lattice parameter to a.u.
        println(f, "\tecutwfc =$(energy_cutoff/ry_to_ev),ecutrho =$(energy_cutoff_rho/ry_to_ev), ")  # convert energy to Rydberg
        println(f, "\toccupations='smearing', smearing='mp', degauss=$smear ")
        println(f, " / ")
        # electrons section
        println(f, " &electrons ")
        println(f, "mixing_mode ='local-TF', electron_maxstep = 250,")
        println(f, "mixing_beta = 0.575, conv_thr = 0.000001,")
        println(f, " / ")
        # atomic species - pseudopotential specification
        println(f, "ATOMIC_SPECIES ")
        println(f, " $el  $el_weight $potential ")
        # atomic positions
        println(f, "ATOMIC_POSITIONS alat ")
        for a in atoms
            println(f, " $el\t$(a[1])\t$(a[2])\t$(a[3]) ")
        end
        # k-points
        println(f, "K_POINTS automatic ")
        println(f, " $kpoints $kpoints 2 0 0 0 ")
        # basis vectors
        println(f, "CELL_PARAMETERS alat ")
        println(f, "$(basis[1,:][1])\t$(basis[1,:][2])\t$(basis[1,:][3]) ")
        println(f, "$(basis[2,:][1])\t$(basis[2,:][2])\t$(basis[2,:][3]) ")
        println(f, "0.0\t0.0\t$(basis[3,3]*10. + vacuum/lp) ") # adds vacuum
    end
    return nothing
end # end write_qe_inputs

function write_vasp_inputs(lp, basis, layers, atoms)
    open("POSCAR", "w") do f # write the POSCAR file
        println(f, el)
        println(f, "$lp ")
        println(f, "$(basis[1,:][1])\t$(basis[1,:][2])\t0.0 ")
        println(f, "$(basis[2,:][1])\t$(basis[2,:][2])\t0.0 ")
        println(f, "0.0\t0.0\t$(basis[3,3]*layers + 20.0/lp) ")
        println(f, "$layers ")
        println(f, "Cartesian ")
        for a in atoms
            println(f, "$(a[1])\t$(a[2])\t$(a[3])\tF\tF\tT ")
        end
    end
    open("KPOINTS", "w") do f # write the KPOINTS file
        println(f, el)
        println(f, "0 ")
        println(f, "Monkhorst-Pack ")
        println(f, "$kpoints $kpoints 1 ")
        println(f, "0 0 0 ")
    end
    return nothing
end # end write_vasp_inputs

function gsfe(structure, lp, slip_system=nothing)
    # creates the summary file
    open("GSFE_SUMMARY", "w") do f
        # pass
        #f.println("Generalized Stacking Fault Energy for {} {}\n".format(structure,el))
        #f.println("=============================================\n")
    end

    # initialize the lists for the energy values
    energy = []; for d in fault_points
        # create the stacking fault structure
        area, fault_length, natoms = create_stackingfault(structure, lp, stacking_layers, d, slip_system)

        # if dft == "qe"
        #     E, walltime = run_qe()
        #     E = E*ry_to_ev
        # elseif dft == "vasp"
        #     if length(energy) == 0
        #         try
        #             `rm IBZKPT CHG CONTCAR DOSCAR EIGENVAL OSZICAR OUTCAR`
        #             `rm PCDAT XDATCAR EIGENVAL vasprun.xml cellvol PROCAR CHGCAR WAVECAR`
        #         catch exc
        #             continue
        #         end
        #     end
        #     start_time = time()
        #     E = run_vasp()
        #     walltime = time() - start_time
        # end

        # push!(energy, E)

        # open("GSFE_SUMMARY", "a") do f
        #     # write the normalized displacement and the relaxed and static energies in mJ/m**2
        #     println(f, "$(d*fault_length)\t$((last(energy) - first(energy)) * 1.60217733e-19 * 1e23 / area)\t$walltime ")
        # end
        `mkdir "./input_gens"`
        `cp gsfe.in ./input_gens/gsfe_$filecount.in`
        open("RE_comm.cm", "a") do f
            println(f, "mpirun -np 16 pw.x -i gsfe_$filecount.in > gsfe_$filecount.out")
        end
        filecount += 1
    end
    return nothing
end # end gsfe

function run_qe()
    # run quantum espresso
    # run(Cmd(["mpirun", "-np", "$num_proc", "pw.x", "-i", "gsfe.in", " > ", "gsfe.out"]))
    cmd_runqe = Cmd(["mpirun", "-np", "$num_proc", "pw.x", "-i", "gsfe.in"])
    run(pipeline(cmd_runqe, stdout="gsfe.out"))

    # get energy
    io = IOBuffer()
    cmd_grep_totalenergy    = `grep '! *[ ] total energy' gsfe.out`
    cmd_awk_print5          = `awk '{print $5}'`
    # pipeline_cmd1 = run(Cmd(["$cmd1", " | ", "$cmd2"]))
    cmd_pl_totalenergy      = run(pipeline(cmd_awk_print5; stdin=cmd_grep_totalenergy, stdout=io))
    # totalenergy_line        = read(cmd_grep_totalenergy, String)
    totalenergy_str         = String(take!(io))

    # get walltime
    io = IOBuffer()
    time = try
        cmd_grep_PWSCF      = `grep 'PWSCF' gsfe.out`
        cmd_tail_lastline   = `tail -n1`
        # pipeline_cmd2 = run(Cmd(["$cmd3", " | ", "$cmd4"]))
        cmd_pl_walltime     = run(pipeline(cmd_tail_lastline; stdin=cmd_grep_PWSCF, stdout=io))
        # walltime_line       = read(cmd_grep_totalenergy, String)
        # walltime_str        = read(cmd_pl_walltime, String)
        walltime_str        = String(take!(io))

        time_str = walltime_str # .communicate()[1]
        time_arr = last(findall(r"(\d+)m[ ]*(\d+).(\d+)s", time_str))
        time_arr = Int.(time_arr)

        if length(time_arr) == 3
            60*time_arr[1] + time_arr[2] + 0.01*time_arr[3]
        else
            time_arr[1] + 0.01*time_arr[2]
        end
    catch
        0.
    end

    return parse(Float64, totalenergy_str), time
end # end run_qe

function run_vasp()
    # run vasp - static run first
    cmd_tail_lastline_OSZICAR   = `tail -n1 OSZICAR`
    cmd_awk_print5              = `awk '{print $5}'`
    tries = 0
    has_run = false
    while tries < 2 && !has_run
        try
            if relax
                # run vasp - first a relaxation, then a static run
                `cp relax.INCAR INCAR`
                msg = `mpirun -np $num_proc ./vasp`
                `cp CONTCAR POSCAR`
            end
            `rm WAVECAR CHGCAR`
            `cp gsfe.INCAR INCAR`
            io = IOBuffer()
            msg = `mpirun -np $num_proc ./vasp`
            # pipeline_cmd1 = `$(cmd1) | $(cmd2)`
            cmd_pl_totalenergy  = run(pipeline(cmd_awk_print5; stdin=cmd_tail_lastline_OSZICAR, stdout=io))
            # totalenergy_line    = read(cmd_tail_lastline_OSZICAR, String)
            # totalenergy_str     = read(cmd_pl_totalenergy, String)
            totalenergy_str     = String(take!(io))
            E = parse(Float64, totalenergy_str)
            has_run = true
        catch exc
            println("VASP run failed, trying again after deleting output files.")
            try
                `rm IBZKPT CHG CONTCAR DOSCAR EIGENVAL OSZICAR OUTCAR`
                `rm PCDAT XDATCAR EIGENVAL vasprun.xml cellvol PROCAR CHGCAR WAVECAR`
            catch e
                continue
            end
            tries += 1
        end
    end
    # get energy
    io = IOBuffer()
    cmd_pl_totalenergy  = run(pipeline(cmd_awk_print5; stdin=cmd_tail_lastline_OSZICAR, stdout=io))
    # totalenergy_line    = read(cmd_tail_lastline_OSZICAR, String)
    # totalenergy_str     = read(cmd_pl_totalenergy, String)
    totalenergy_str     = String(take!(io))

    return parse(Float64, totalenergy_str)


    #if msg != 0:
    #	print("VASP exited with an error. It may only be able to run on the CAVS cluster machines.")
end # end run_vasp



# get inputs from the command line
argc        = length(ARGS)
structure   = argc >= 1 ? ARGS[1]                   : nothing
latp        = argc >= 2 ? parse(Float64, ARGS[2])   : nothing
slip        = argc >= 3 ? ARGS[3]                   : nothing
extent      = argc >= 4 ? parse(Float64, ARGS[4])   : nothing

# if the inputs are not specified, prompt for them
#if not element:
#	print("Enter the element name:")
#	element = raw_input("> ")
if isnothing(structure)
    i_struct = nothing; while isnothing(i_struct)
        println("== Enter the desired structure ==")
        println("= 1) FCC                        =")
        println("= 2) BCC                        =")
        println("=================================")
        print("> "); global i_struct = parse(Int64, readline())
        global structure = if i_struct == 1
            println("Structure type: FCC")
            "fcc"
        elseif i_struct == 2
            println("Structure type: BCC")
            "bcc"
        else
            println(" Only FCC and BCC systems are implemented currently ")
            nothing
        end
    end
end
slip_choices = if lowercase(structure) == "fcc"
    ["(111)[1-10] (full)",
        "(111)[11-2] (partial)"]
elseif lowercase(structure) == "bcc"
    ["(110)[-111] (full)",
        "(110)[001] (partial)",
        "(110)[-110] (longpartial)"]
else
    println("Structures other than FCC and BCC are not currently supported")
    exit(1)
end
if isnothing(latp)
    println("Enter the equilibrium lattice parameter for your element:")
    print("> "); global latp = parse(Float64, readline())
end
if isnothing(slip)
    i_slip = nothing; while isnothing(i_slip)
        println("====== Enter the desired slip system ======")
        for (i, sc) in enumerate(slip_choices)
            println("= $i) $(sc) =")
        end
        println("===========================================")
        print("> "); global i_slip = parse(Int64, readline())
        global slip = if i_slip == 1
            println("Using burgers vector direction")
            "full"
        elseif i_slip == 2
            println("Using partial dislocation direction")
            "partial"
        elseif i_slip == 3 && structure == "bcc"
            println("Using longer partial dislocation direction")
            "longpartial"
        else
            println("\nChoice not recognized")
            nothing
        end
    end
end
extent = isnothing(extent) ? 1.0 : nothing

gsfe(structure, latp, slip)
