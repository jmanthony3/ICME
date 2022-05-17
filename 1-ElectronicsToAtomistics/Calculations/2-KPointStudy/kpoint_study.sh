#!/bin/bash -x



######################### INPUT FILE ##########################
### define input file variables
ELEMENT_NAME="Fe"
ELEMENT_AMU=55.845
PSEUDOPOTENTIAL_FILENAME="Fe.pbe-spn-kjpaw_psl.0.2.1.UPF"
ELEMENT_POS=(
    0.00 # x
    0.00 # y
    0.00 # z
)
REFERENCE_STRUCTURE="bcc" # body-centered cubic
LATTICE_PARAMETER=2.866 # angstrom
MAX_ITER=500
MIXING_BETA=0.5
# define according to: seq FIRST STEP LAST
declare -a CUTOFF_ENERGIES=($(seq 30 15 120))
# define according to: seq FIRST STEP LAST
declare -a KPOINTS=($(seq 1 1 12)) # k-points for convergence
NUM_PROC=$(nproc) # grabs all cores available by default



########################## `ev_curve` #########################
# input arguments to script
# also uses `REFERENCE_STRUCTURE` and `LATTICE_PARAMETER`
# 1=Birch(First), 2=Birch(Second), 3=Keane, 4=Murnaghan
EQ_OF_STATE=4



####################### `EvA_EvV_plot.py` #####################
ENERGY_OFFSET=4473.05298206 # Ry




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - END OF HUMAN EDITABLE SECTION - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #




set +x
### automatically define other `.in` file variables from inputs
reference_structure=$(echo $REFERENCE_STRUCTURE | tr '[:upper:]' '[:lower:]')
if [[ "$reference_structure" == "fcc" ]]; then
    ibrav=2
elif [[ "$reference_structure" == "bcc" ]]; then
    ibrav=3
else
    echo "Variable REFERENCE_STRUCTURE=$REFERENCE_STRUCTURE not understood. Must be either 'fcc' or 'bcc'."
    exit
fi
echo "Based on $reference_structure, ibrav=$ibrav"
# converts to bohr and sets to large value for offset
lattice_parameter_bohr=$(echo "$LATTICE_PARAMETER*1.88973" | bc -l) # bohr


### automatically define other `ev_curve` file variables from inputs
sed -i "64s%^[[:digit:]]*[^ #]*%$EQ_OF_STATE%" "../0-Scripts/ev_curve"


### automatically define other `EvA_EvV_plot.py` file variables from inputs
sed -i "s%^energy_offset = [[:digit:]]*\.*[[:digit:]]*[^ #]%energy_offset = $ENERGY_OFFSET%" "../0-Scripts/EvA_EvV_plot.py"



#################### PERFORM K-POINT STUDY ####################
### examine for every cutoff energy
for cutoff_energy in "${CUTOFF_ENERGIES[@]}"; do
    echo "Examining 'ecutwfc=$cutoff_energy' Ry..."
    # to store results in directory
    mkdir "$cutoff_energy"
    cd "$cutoff_energy"
    rm "./EvsK" "./AvsK" "./GvsK" "./TvsK"
    touch "./EvsK" "./AvsK" "./GvsK" "./TvsK"
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
 ${kpoints[0]} ${kpoints[1]} ${kpoints[2]} ${kpoints_shift[0]} ${kpoints_shift[1]} ${kpoints_shift[2]}""" > "$input_filename.in"

        ### execute QE with input parameters
        echo "Executing QE according to $input_filename.in..."
        (set -x; mpirun -np $NUM_PROC pw.x -in "$input_filename.in" > "$input_filename.out")


        ### run Fortran codes on input files
        # compiles `evfit.f` outputs `evfit`
        gfortran -O2 "../../0-Scripts/evfit.f" -o "../../0-Scripts/evfit"
        chmod +x "../../0-Scripts/ev_curve" # makes file executable
        cp "$input_filename.in" "../../0-Scripts/$REFERENCE_STRUCTURE.ev.in" # create appropriate input file to `ev_curve`
        rm -r "temp/"
        cd "../../0-Scripts"
        sed -i "s%pseudo_dir = '\.\./\.\./',%pseudo_dir = '\.\./'%" "$REFERENCE_STRUCTURE.ev.in"
        # this outputs `evfit.4`: reference structure, lattice parameter
        ./ev_curve $REFERENCE_STRUCTURE $LATTICE_PARAMETER
        python3 "EvA_EvV_plot.py" # generate plots
        mv "$REFERENCE_STRUCTURE.ev.in" "../2-KPointStudy/$cutoff_energy/$REFERENCE_STRUCTURE.ev.in" # create appropriate input file to `ev_curve`
        mv "evfit" "../2-KPointStudy/$cutoff_energy/evfit"
        mv "EvsA" "../2-KPointStudy/$cutoff_energy/EvsA.$cutoff_energy.$K"
        mv "EvsV" "../2-KPointStudy/$cutoff_energy/EvsV.$cutoff_energy.$K"
        mv "SUMMARY" "../2-KPointStudy/$cutoff_energy/SUMMARY.$cutoff_energy.$K"
        mv "evfit.4" "../2-KPointStudy/$cutoff_energy/evfit.4"
        mv "pw_ev.out" "../2-KPointStudy/$cutoff_energy/pw_ev.out"
        mv "Name_of_EvA.pdf" "../2-KPointStudy/$cutoff_energy/Name_of_EvA.pdf"
        mv "Name_of_EvV.pdf" "../2-KPointStudy/$cutoff_energy/Name_of_EvV.pdf"
        mv "Name_of_Combined.pdf" "../2-KPointStudy/$cutoff_energy/Name_of_Combined.pdf"
        rm -r "temp/"
        cd "../2-KPointStudy/$cutoff_energy/"
        declare -a inputs=(
            "EvsA.$cutoff_energy.$K"
            "EvsV.$cutoff_energy.$K"
        )
        for input in "${inputs[@]}"; do
            echo "Opening $input to offset by $ENERGY_OFFSET eV..."
            i=1
            readarray file < $input
            for line in "${file[@]}"; do
                IFS=" " read -a elem <<< "$line"
                energy=${elem[1]}
                offset_energy=$(echo "$energy+$ENERGY_OFFSET" | bc)
                sed -i "${i}s% \-[[:digit:]]*\.[[:digit:]]*% $offset_energy%" "$input"
                i=$(echo "$i+1" | bc)
            done # end line `i`
            echo "Closing $input..."
        done # end of `input` file
        echo "Adjusting summary files by $ENERGY_OFFSET eV..."
        energy=$(sed -n "s%^Equilibrium Energy per Atom    = %%p" "SUMMARY.$cutoff_energy.$K")
        offset_energy=$(echo "$energy+$ENERGY_OFFSET" | bc)
        sed -i "s%^Equilibrium Energy per Atom    = \-[[:digit:]]*\.[[:digit:]]*%Equilibrium Energy per Atom    = $offset_energy%" "SUMMARY.$cutoff_energy.$K"
        bulk=$(sed -n "s%^Bulk Modulus (kbar)            = %%p" "SUMMARY.$cutoff_energy.$K")
        bulk_gpa=$(echo "scale=9;$bulk/10" | bc)
        sed -i "s%^Bulk Modulus (kbar)            = [[:digit:]]*\.[[:digit:]]*%Bulk Modulus (GPa)             = $bulk_gpa%" "SUMMARY.$cutoff_energy.$K"
        # touch "./EvsK" "./AvsK" "./GvsK" "./TvsK"
        echo "$K $offset_energy" >> "./EvsK"
        alat=$(sed -n "s%^Equilibrium lattice constant   = %%p" "SUMMARY.$cutoff_energy.$K")
        echo "$K $alat" >> "./AvsK"
        echo "$K $bulk_gpa" >> "./GvsK"
        conv_time=$(sed -n "s%^Total RUN time (sec)           = %%p" "SUMMARY.$cutoff_energy.$K")
        echo "$K $conv_time" >> "./TvsK"
    done # end of `K`-point
    cd "../"
done # end of `cutoff_energy`




# that's all folks