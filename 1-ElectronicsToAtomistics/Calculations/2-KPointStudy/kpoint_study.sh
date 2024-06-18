#!/bin/bash -x



######################### INPUT FILE ##########################
### define input file variables
ELEMENT_NAME="Fe" # Periodic Table identifier of element
ELEMENT_AMU=55.845 # [g/mol], literature value for element
PSEUDOPOTENTIAL_FILENAME="Fe.pbe-spn-kjpaw_psl.0.2.1.UPF"
ELEMENT_POS=(
    0.00 # x
    0.00 # y
    0.00 # z
)
REFERENCE_STRUCTURE="bcc" # body-centered cubic
LATTICE_PARAMETER=2.866 # [angstrom], literature value for element
MAX_ITER=500
MIXING_BETA=0.5
# define according to: seq FIRST STEP LAST
declare -a CUTOFF_ENERGIES=($(seq 30 15 120)) # [Ry]
# define according to: seq FIRST STEP LAST
declare -a KPOINTS=($(seq 1 1 12)) # k-points for convergence
# working language to perform calculations and plot results
COMPUTING_LANGUAGE="Julia" # can also be "Python"
NUM_PROC=$(nproc) # grabs all cores available by default



########################## `ev_curve` #########################
# input arguments to script
# also uses `REFERENCE_STRUCTURE` and `LATTICE_PARAMETER`
# 1=Birch(First), 2=Birch(Second), 3=Keane, 4=Murnaghan
EQ_OF_STATE=4



####################### `EvA_EvV_plot.py` #####################
ENERGY_OFFSET=4473.05298206 # [Ry]





# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - END OF HUMAN EDITABLE SECTION - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #





set +x # turn script tracing off
computing_language=$(echo $COMPUTING_LANGUAGE | tr '[:upper:]' '[:lower:]')
if [[ "$computing_language" == "julia" ]]; then
    cl_ext="jl"
elif [[ "$computing_language" == "python" ]]; then
    cl_ext="py"
else
    echo "Variable COMPUTING_LANGUAGE=$COMPUTING_LANGUAGE \
        not understood. Must be either 'Julia' or 'Python'."
    exit
fi



######################## MODIFY SCRIPTS #######################
### automatically define other `.in` file variables from inputs
reference_structure=$(echo $REFERENCE_STRUCTURE | tr '[:upper:]' '[:lower:]')
if [[ "$reference_structure" == "fcc" ]]; then
    ibrav=2
elif [[ "$reference_structure" == "bcc" ]]; then
    ibrav=3
else
    echo "Variable REFERENCE_STRUCTURE=$REFERENCE_STRUCTURE \
        not understood. Must be either 'fcc' or 'bcc'."
    exit
fi
echo "Based on $reference_structure, ibrav=$ibrav"
# '*1.88973' converts [angstrom] to [bohr]
lattice_parameter_bohr=$(echo "$LATTICE_PARAMETER*1.88973" | bc -l)


### automatically define other `ev_curve` file variables from inputs
sed -i "64s%^[[:digit:]]*[^ #]%$EQ_OF_STATE%" "../0-Scripts/ev_curve"


### automatically define other `EvA_EvV_plot` file variables from inputs
sed -i "s%^energy_offset = [[:digit:]]*\.*[[:digit:]]*[^ #]%energy_offset = $ENERGY_OFFSET%" \
    "../0-Scripts/EvA_EvV_plot.$cl_ext"



#################### PERFORM K-POINT STUDY ####################
### examine for every cutoff energy
for cutoff_energy in "${CUTOFF_ENERGIES[@]}"; do
    echo "Examining 'ecutwfc=$cutoff_energy' Ry..."
    mkdir "$cutoff_energy" # to store results in directory
    # move into directory for `cutoff_energy`
    cd "$cutoff_energy"
    # create data files
    rm "./EvsK.dat" "./AvsK.dat" "./GvsK.dat" "./TvsK.dat"
    touch "./EvsK.dat" "./AvsK.dat" "./GvsK.dat" "./TvsK.dat"
    # examine for every k-point
    for K in "${KPOINTS[@]}"; do
        echo "Examining $K kpoints..."
        kpoints=( # number of k-points (think as mesh grid element count)
            $K # x
            $K # y
            $K # z
        )
        kpoints_shift=( # shift of `KPOINTS` in some direction
            0 # x
            0 # y
            0 # z
        )
        input_filename="$ELEMENT_NAME.$cutoff_energy.$K" # name input file
        # define input file
        echo """ &control
    prefix=''
    outdir='temp'
    pseudo_dir = '../../',
 /
 &system
    ibrav=  $ibrav, celldm(1) =$lattice_parameter_bohr, nat=  1, ntyp=  1,
    ecutwfc =$cutoff_energy,
    occupations='smearing', smearing='mp', degauss=0.06
 /
 &electrons
    electron_maxstep =$MAX_ITER,
    mixing_beta =$MIXING_BETA,
 /
ATOMIC_SPECIES
 $ELEMENT_NAME  $ELEMENT_AMU $PSEUDOPOTENTIAL_FILENAME
ATOMIC_POSITIONS (alat)
 $ELEMENT_NAME ${ELEMENT_POS[0]} ${ELEMENT_POS[1]} ${ELEMENT_POS[2]}
K_POINTS (automatic)
 ${kpoints[0]} ${kpoints[1]} ${kpoints[2]} ${kpoints_shift[0]} ${kpoints_shift[1]} ${kpoints_shift[2]}\
        """ > "$input_filename.in"

        ### execute QE with input parameters
        echo "Executing QE according to $input_filename.in..."
        (set -x;
            mpirun -np $NUM_PROC pw.x -in "$input_filename.in" \
                > "$input_filename.out" 2> /dev/null
        )

        ### run Fortran codes on input files
        # create appropriate input file to `ev_curve`
        cp "$input_filename.in" "../../0-Scripts/$REFERENCE_STRUCTURE.ev.in"
        rm -r "temp/" # remove calculations temporary folder
        cd "../../0-Scripts" # move into Scripts folder
        gfortran -std=legacy -O2 "evfit.f" -o "evfit" 2> /dev/null # compiles `evfit.f` outputs `evfit`
        sed -i "s%pseudo_dir = '\.\./\.\./',%pseudo_dir = '\.\./'%" \
            "$REFERENCE_STRUCTURE.ev.in"
        # this outputs `evfit.4`: reference structure, lattice parameter
        ./ev_curve $REFERENCE_STRUCTURE $LATTICE_PARAMETER 2>/dev/null
        if [[ "$computing_language" == "julia" ]]; then
            julia "EvA_EvV_plot.jl" # generate plots
        elif [[ "$computing_language" == "python" ]]; then
            python3 "EvA_EvV_plot.py" # generate plots
        else
            echo "Variable COMPUTING_LANGUAGE=$COMPUTING_LANGUAGE \
                not understood. Must be either 'Julia' or 'Python'."
            exit
        fi
        # move all output files back to `cutoff_energy` directory
        # create appropriate input file to `ev_curve`
        mv "$REFERENCE_STRUCTURE.ev.in" "../2-KPointStudy/$cutoff_energy/"
        mv "evfit" "../2-KPointStudy/$cutoff_energy/"
        mv "EvsA" "../2-KPointStudy/$cutoff_energy/EvsA.$cutoff_energy.$K"
        mv "EvsV" "../2-KPointStudy/$cutoff_energy/EvsV.$cutoff_energy.$K"
        mv "SUMMARY" "../2-KPointStudy/$cutoff_energy/SUMMARY.$cutoff_energy.$K"
        mv "evfit.4" "../2-KPointStudy/$cutoff_energy/"
        mv "pw_ev.out" "../2-KPointStudy/$cutoff_energy/"
        mv "Name_of_EvA.png" "../2-KPointStudy/$cutoff_energy/"
        mv "Name_of_EvV.png" "../2-KPointStudy/$cutoff_energy/"
        mv "Name_of_Combined.png" "../2-KPointStudy/$cutoff_energy/"
        rm -r "temp/" # remove calculations temporary folder

        ### post-process output files
        # move back to `cutoff_energy` directory
        cd "../2-KPointStudy/$cutoff_energy/"
        declare -a inputs=(
            "EvsA.$cutoff_energy.$K"
            "EvsV.$cutoff_energy.$K"
        )
        for input in "${inputs[@]}"; do
            echo -e -n "Opening $input to offset by $ENERGY_OFFSET eV...\r"
            i=1
            readarray file < $input
            for line in "${file[@]}"; do
                IFS=" " read -a elem <<< "$line"
                energy=${elem[1]}
                offset_energy=$(echo "$energy+$ENERGY_OFFSET" | bc)
                sed -i "${i}s% \-[[:digit:]]*\.[[:digit:]]*% $offset_energy%" \
                    "$input"
                i=$(echo "$i+1" | bc) # end line `i`
            done # end of `input` file
            echo -e -n "Closing $input...\r"
        done # end of processing
        ## adjust summary file
        echo -e -n "Adjusting summary files by $ENERGY_OFFSET eV...\r"
        # previous energy value
        energy=$(sed -n "s%^Equilibrium Energy per Atom    = %%p" "SUMMARY.$cutoff_energy.$K")
        # offset energy value
        offset_energy=$(echo "$energy+$ENERGY_OFFSET" | bc)
        # replace energy value
        sed -i "s%^Equilibrium Energy per Atom    = \-[[:digit:]]*\.[[:digit:]]*%Equilibrium Energy per Atom    = $offset_energy%" \
            "SUMMARY.$cutoff_energy.$K"
        # get lattice parameter [angstrom]
        lattice_parameter=$(sed -n "s%^Equilibrium lattice constant   = %%p" "SUMMARY.$cutoff_energy.$K")
        # previous bulk modulus [kbar]
        bulk=$(sed -n "s%^Bulk Modulus (kbar)            = %%p" "SUMMARY.$cutoff_energy.$K")
        # convert to [GPa]
        bulk_gpa=$(echo "scale=9;$bulk/10" | bc)
        # replace bulk modulus
        sed -i "s%^Bulk Modulus (kbar)            = [[:digit:]]*\.[[:digit:]]*%Bulk Modulus (GPa)             = $bulk_gpa%" \
            "SUMMARY.$cutoff_energy.$K"
        ## append data files
        echo "$K $offset_energy" >> "./EvsK.dat"
        echo "$K $lattice_parameter" >> "./AvsK.dat"
        echo "$K $bulk_gpa" >> "./GvsK.dat"
        conv_time=$(sed -n "s%^Total RUN time (sec)           = %%p" "SUMMARY.$cutoff_energy.$K")
        conv_iter=$(sed -n "s%^[[:blank:]]*convergence has been achieved in[[:blank:]]*%%p" \
            "$input_filename.out" | sed "s% iterations$%%")
        echo "$K $conv_time $conv_iter" >> "./TvsK.dat"
        echo -e -n "                                                        \r"
        # end of `K`-point
    done # end of `cutoff_energy`
    cd "../"
done # end of study



### visualize data files
if [[ "$computing_language" == "julia" ]]; then
    julia "which_ecutwfc.jl" # generate plots
elif [[ "$computing_language" == "python" ]]; then
    python3 "which_ecutwfc.py" # generate plots
else
    echo "Variable COMPUTING_LANGUAGE=$COMPUTING_LANGUAGE \
        not understood. Must be either 'Julia' or 'Python'."
    exit
fi





# that's all folks