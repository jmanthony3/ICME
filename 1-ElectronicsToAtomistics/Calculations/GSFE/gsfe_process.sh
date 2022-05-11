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
CUTOFF_ENERGY=64
MAX_ITER=500
KPOINT=8 # k-points for convergence
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



########################## `ev_curve` #########################
# input arguments to script
# also uses `REFERENCE_STRUCTURE` and `LATTICE_PARAMETER`
# 1=Birch(First), 2=Birch(Second), 3=Keane, 4=Murnaghan
EQ_OF_STATE=4



####################### `EvA_EvV_plot.py` #####################
# ENERGY_OFFSET=4473.05298206 # Ry
ENERGY_OFFSET=4479.41619611 # Ry



##################### `OutputFileCreator.py` ##################
DISLOCATION_GRADE="full"




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
LATTICE_AREA=$(echo "$LATTICE_PARAMETER*$LATTICE_PARAMETER" | bc -l) # 
LATTICE_PARAMETER_BOHR=$(echo "$LATTICE_PARAMETER*1.88973" | bc -l) # bohr
INPUT_FILENAME=$ELEMENT_NAME # name input file


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


### automatically define other `ev_curve` file variables from inputs
sed -i "64s%^[[:digit:]]*[^ #]*%$EQ_OF_STATE%" "../Scripts/ev_curve"


### automatically define other `EvA_EvV_plot.py` file variables from inputs
sed -i "s%^energy_offset = [[:digit:]]*\.*[[:digit:]]*[^ #]%energy_offset = $ENERGY_OFFSET%" "../Scripts/EvA_EvV_plot.py"


### automatically define other `OutputFileCreator.py` file variables from inputs
case $DISLOCATION_GRADE in

    "full" | "FULL")
        DISP_SCALE=$(echo "scale=9; sqrt(3) / 2" | bc)
        ;;

    "longpartial" | "LONGPARTIAL")
        DISP_SCALE=$(echo "scale=9; sqrt(2) / 2" | bc)
        ;;

    "partial" | "PARTIAL")
        DISP_SCALE=$(echo "scale=9; 1" | bc)
        ;;
esac

# modify script
sed -i "s%^num_proc = [[:digit:]]*[^ #]*%num_proc = 16%" "../Scripts/OutputFileCreator.py"
sed -i "s%^el = '[[:print:]]*'[^ #]*%el = '$ELEMENT_NAME'%" "../Scripts/OutputFileCreator.py"
sed -i "s%^potential = '[[:print:]]*'[^ #]*%potential = '$PSEUDOPOTENTIAL_FILENAME'%" "../Scripts/OutputFileCreator.py"
sed -i "s%^el_weight = [[:digit:]]*\.*[[:digit:]]*[^ #]*%el_weight = $ELEMENT_AMU%" "../Scripts/OutputFileCreator.py"
sed -i "s%^energy_cutoff = [[:digit:]]*\.*[[:digit:]]*\*13.6057[^ #]*%energy_cutoff = $CUTOFF_ENERGY\*13.6057%" "../Scripts/OutputFileCreator.py"
sed -i "s%^kpoints = [[:digit:]]*[^ #]*%kpoints = $KPOINT%" "../Scripts/OutputFileCreator.py"

# modify script
sed -i "s%^num_proc = [[:digit:]]*[^ #]*%num_proc = 16%" "../Scripts/OutputSummarizer.py"
sed -i "s%^el = '[[:print:]]*'[^ #]*%el = '$ELEMENT_NAME'%" "../Scripts/OutputSummarizer.py"
sed -i "s%^potential = '[[:print:]]*'[^ #]*%potential = '$PSEUDOPOTENTIAL_FILENAME'%" "../Scripts/OutputSummarizer.py"
sed -i "s%^el_weight = [[:digit:]]*\.*[[:digit:]]*[^ #]*%el_weight = $ELEMENT_AMU%" "../Scripts/OutputSummarizer.py"
sed -i "s%^energy_cutoff = [[:digit:]]*\.*[[:digit:]]*\*13.6057[^ #]*%energy_cutoff = $CUTOFF_ENERGY\*13.6057%" "../Scripts/OutputSummarizer.py"
sed -i "s%^kpoints = [[:digit:]]*[^ #]*%kpoints = $KPOINT%" "../Scripts/OutputSummarizer.py"



###################### GENERATE GSFE DATA #####################
### execute QE with input parameters
mpirun -np $NUM_PROC pw.x -in "$INPUT_FILENAME.in" > "$ELEMENT_NAME.out"


### run Fortran codes on input files
gfortran -O2 "../Scripts/evfit.f" -o "../Scripts/evfit" # compiles `evfit.f` outputs `evfit`
chmod +x "../Scripts/ev_curve" # makes file executable
cp "$INPUT_FILENAME.in" "$REFERENCE_STRUCTURE.ev.in" # create appropriate input file to `ev_curve`
cp "$INPUT_FILENAME.in" "../Scripts/$REFERENCE_STRUCTURE.ev.in" # create appropriate input file to `ev_curve`
cd "../Scripts"
# this outputs `evfit.4`: reference structure, lattice parameter
./ev_curve $REFERENCE_STRUCTURE $LATTICE_PARAMETER
cp "EvsA" "../GSFE/EvsA"
cp "EvsV" "../GSFE/EvsV"
cp "SUMMARY" "../GSFE/SUMMARY"
python3 "EvA_EvV_plot.py" # generate plots
cp "Name_of_EvA.pdf" "../GSFE/Name_of_EvA.pdf"
cp "Name_of_EvV.pdf" "../GSFE/Name_of_EvV.pdf"
cp "Name_of_Combined.pdf" "../GSFE/Name_of_Combined.pdf"
rm "EvsA" "EvsV" "SUMMARY"
rm "$REFERENCE_STRUCTURE.ev.in" "evfit.4" "gsfe.in" "pw_ev.out"
rm "Name_of_EvA.pdf" "Name_of_EvV.pdf" "Name_of_Combined.pdf"
rm -r "temp"
cd "../GSFE"
declare -a inputs=(
    "EvsA"
    "EvsV"
)
for input in "${inputs[@]}"; do
    i=1
    init_energy=$(sed -n "${i}s%^[[:digit:]]\.[[:digit:]]* %%p" "$input")
    init_energy=$(echo "$init_energy+$ENERGY_OFFSET" | bc)
    init_energy=$(sed -n "s%-*%%p" <<< $init_energy)
    readarray file < $input
    for line in "${file[@]}"; do
        IFS=" " read -a elem <<< "$line"
        energy=${elem[1]}
        offset_energy=$(echo "$energy+$init_energy+$ENERGY_OFFSET" | bc)
        sed -i "${i}s% \-[[:digit:]]*\.[[:digit:]]*% $offset_energy%" "$input"
        i=$(echo "$i+1" | bc)
    done # end line `i`
done # end of `input` file
energy=$(sed -n "s%^Equilibrium Energy per Atom    = %%p" "SUMMARY")
offset_energy=$(echo "$energy+$init_energy+$ENERGY_OFFSET" | bc)
sed -i "s%^Equilibrium Energy per Atom    = \-[[:digit:]]*\.[[:digit:]]*%Equilibrium Energy per Atom    = $offset_energy%" "SUMMARY"
lattice_parameter=$(sed -n "s%^Equilibrium lattice constant   = %%p" "SUMMARY")
bulk=$(sed -n "s%^Bulk Modulus (kbar)            = %%p" "SUMMARY")
bulk_gpa=$(echo "scale=9;$bulk/10" | bc)
sed -i "s%^Bulk Modulus (kbar)            = [[:digit:]]*\.[[:digit:]]*%Bulk Modulus (GPa)             = $bulk_gpa%" "SUMMARY"

cp "output_gens/gsfe_"*".out" "../Scripts"

cd "../Scripts"
# reference structure, lattice parameter, and block motion
python2 "OutputSummarizer.py" $REFERENCE_STRUCTURE $LATTICE_PARAMETER $DISLOCATION_GRADE
cp "GSFE_SUMMARY" "../GSFE/GSFE_SUMMARY"
rm "GSFE_SUMMARY" "gsfe_"*".out"
cd "../GSFE"
declare -a inputs=(
    "GSFE_SUMMARY"
)
for input in "${inputs[@]}"; do
    i=1
    readarray file < $input
    for line in "${file[@]}"; do
        read -a elem <<< "$line"
        energy=${elem[1]}
        scaled_energy=$(echo "scale=9; $energy/$LATTICE_AREA" | bc)
        sed -i "${i}s%[[:space:]]\-*[[:digit:]]*\.[[:digit:]]*%\t$scaled_energy%" "$input"
        i=$(echo "$i+1" | bc)
    done
done
cp "GSFE_SUMMARY" "GSFE_SUMMARY-Calibration"
declare -a inputs=(
    "GSFE_SUMMARY-Calibration"
)
for input in "${inputs[@]}"; do
    i=1
    readarray file < $input
    for line in "${file[@]}"; do
        read -a elem <<< "$line"
        displacement=${elem[0]}
        scaled_displacement=$(echo "$displacement*$lattice_parameter" | bc)
        sed -i "${i}s%^[[:digit:]]*\.[[:digit:]]*\t%$scaled_displacement\t%" "$input"
        sed -i "${i}s%\t${elem[2]} %%" "$input"
        i=$(echo "$i+1" | bc)
    done
    sed -i "s%\t% %" "$input"
done