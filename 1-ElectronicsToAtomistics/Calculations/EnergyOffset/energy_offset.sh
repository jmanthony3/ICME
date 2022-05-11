#!/bin/bash -x



######################### INPUT FILE ##########################
### define input file variables
ELEMENT_NAME="Fe" # Periodic Table identifier of element
ELEMENT_AMU=55.845 # literature value for element
PSEUDOPOTENTIAL_FILENAME="Fe.pbe-spn-kjpaw_psl.0.2.1.UPF"
ELEMENT_POS=(
    0.00 # x
    0.00 # y
    0.00 # z
)
REFERENCE_STRUCTURE="bcc" # body-centered cubic
LATTICE_PARAMETER=2.866 # angstrom, literature value for element
CUTOFF_ENERGY=90
MAX_ITER=500
KPOINT=12 # k-points for convergence
KPOINTS=( # number of k-points (think as mesh grid element count)
    $KPOINT # x
    $KPOINT # y
    $KPOINT # z
)
KPOINTS_SHIFT=( # shift of `KPOINTS` in some direction
    0 # x
    0 # y
    0 # z
)
NUM_PROC=4




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - END OF HUMAN EDITABLE SECTION - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #




### automatically define other `.in` file variables from inputs
case $REFERENCE_STRUCTURE in

    "fcc" | "FCC")
        IBRAV=2
        ;;

    "bcc" | "BCC")
        IBRAV=3
        ;;
esac
# converts to bohr and sets to large value for offset
LATTICE_PARAMETER_BOHR=$(echo "$LATTICE_PARAMETER*1.88973*4" | bc -l) # bohr
INPUT_FILENAME="${ELEMENT_NAME}_offset" # name input file


### write input file
# input file
echo """ &control
    prefix=''
    outdir='temp'
    pseudo_dir = '../',
 /
 &system
    ibrav=  $IBRAV, celldm(1) =$LATTICE_PARAMETER_BOHR, nat=  1, ntyp=  1,
    ecutwfc =$CUTOFF_ENERGY,
    occupations='smearing', smearing='mp', degauss=0.06
 /
 &electrons
    electron_maxstep =$MAX_ITER,
    mixing_beta =0.5,
 /
ATOMIC_SPECIES
 $ELEMENT_NAME  $ELEMENT_AMU $PSEUDOPOTENTIAL_FILENAME
ATOMIC_POSITIONS (alat)
 $ELEMENT_NAME ${ELEMENT_POS[0]} ${ELEMENT_POS[1]} ${ELEMENT_POS[2]}
K_POINTS (automatic)
 ${KPOINTS[0]} ${KPOINTS[1]} ${KPOINTS[2]} ${KPOINTS_SHIFT[0]} ${KPOINTS_SHIFT[1]} ${KPOINTS_SHIFT[2]}""" > "$INPUT_FILENAME.in"



################### CALCULATE OFFSET ENERGY ###################
### execute QE with input parameters
mpirun -np $NUM_PROC pw.x -in "$INPUT_FILENAME.in" > "$INPUT_FILENAME.out"


### grab `total energy`, which is in Ry, from `.out` file
energy_offset=$(
    # grabs only the Ry value
    sed -n "s%\![[:space:]]*total energy[[:space:]]* = [[:space:]]*%%p" "$INPUT_FILENAME.out" | sed "s% Ry$%%"
)
energy_offset=$(echo "$energy_offset*-13.6057" | bc) # '*-13.6057' converts Ry to eV
echo "Energy offset found to be $energy_offset eV."


### replace offset energies in companion scripts
sed -i "s%^ENERGY_OFFSET=[[:digit:]]*\.*[[:digit:]]*[^ #]%ENERGY_OFFSET=$energy_offset%" "../KPointStudy/kpoint_study.sh"
sed -i "s%^energy_offset = [[:digit:]]*\.*[[:digit:]]*[^ #]%energy_offset = $energy_offset%" "../Scripts/EvA_EvV_plot.py"
sed -i "s%^ENERGY_OFFSET=[[:digit:]]*\.*[[:digit:]]*[^ #]%ENERGY_OFFSET=$energy_offset%" "../GSFE/gsfe_create.sh"
sed -i "s%^ENERGY_OFFSET=[[:digit:]]*\.*[[:digit:]]*[^ #]%ENERGY_OFFSET=$energy_offset%" "../GSFE/gsfe_process.sh"


### run Fortran codes on input files
gfortran -O2 "../Scripts/evfit.f" -o "../Scripts/evfit" # compiles `evfit.f` outputs `evfit`
chmod +x "../Scripts/ev_curve" # makes file executable
cp "$INPUT_FILENAME.in" "$REFERENCE_STRUCTURE.ev.in" # create appropriate input file to `ev_curve`
cp "$INPUT_FILENAME.in" "../Scripts/$REFERENCE_STRUCTURE.ev.in" # create appropriate input file to `ev_curve`
cd "../Scripts"
# this outputs `evfit.4`: reference structure, lattice parameter
./ev_curve $REFERENCE_STRUCTURE $LATTICE_PARAMETER
cp "EvsA" "../EnergyOffset/EvsA_offset"
cp "EvsV" "../EnergyOffset/EvsV_offset"
cp "SUMMARY" "../EnergyOffset/SUMMARY_offset"
python3 "EvA_EvV_plot.py" # generate plots
cp "Name_of_EvA.pdf" "../EnergyOffset/Name_of_EvA.pdf"
cp "Name_of_EvV.pdf" "../EnergyOffset/Name_of_EvV.pdf"
cp "Name_of_Combined.pdf" "../EnergyOffset/Name_of_Combined.pdf"
rm "EvsA" "EvsV" "SUMMARY"
rm "$REFERENCE_STRUCTURE.ev.in" "evfit.4" "pw_ev.out"
rm "Name_of_EvA.pdf" "Name_of_EvV.pdf" "Name_of_Combined.pdf"
rm -r "temp"
cd "../EnergyOffset"
declare -a inputs=(
    "EvsA_offset"
    "EvsV_offset"
)
for input in "${inputs[@]}"; do
    i=1
    readarray file < $input
    for line in "${file[@]}"; do
        IFS=" " read -a elem <<< "$line"
        energy=${elem[1]}
        offset_energy=$(echo "$energy+$energy_offset" | bc)
        sed -i "${i}s% \-[[:digit:]]*\.[[:digit:]]*% $offset_energy%" "$input"
        i=$(echo "$i+1" | bc)
    done # end line `i`
done # end of `input` file
energy=$(sed -n "s%^Equilibrium Energy per Atom    = %%p" "SUMMARY_offset")
offset_energy=$(echo "$energy+$energy_offset" | bc)
sed -i "s%^Equilibrium Energy per Atom    = \-[[:digit:]]*\.[[:digit:]]*%Equilibrium Energy per Atom    = $offset_energy%" "SUMMARY_offset"
lattice_parameter=$(sed -n "s%^Equilibrium lattice constant   = %%p" "SUMMARY_offset")
bulk=$(sed -n "s%^Bulk Modulus (kbar)            = %%p" "SUMMARY_offset")
bulk_gpa=$(echo "scale=9;$bulk/10" | bc)
sed -i "s%^Bulk Modulus (kbar)            = [[:digit:]]*\.[[:digit:]]*%Bulk Modulus (GPa)             = $bulk_gpa%" "SUMMARY_offset"

cat "SUMMARY_offset"
