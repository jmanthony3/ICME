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
CUTOFF_ENERGY=64 # Ry
MAX_ITER=500
MIXING_BETA=0.5
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
NUM_PROC=$(nproc) # grabs all cores available by default



########################## `ev_curve` #########################
# input arguments to script
# also uses `REFERENCE_STRUCTURE` and `LATTICE_PARAMETER`
# 1=Birch(First), 2=Birch(Second), 3=Keane, 4=Murnaghan
EQ_OF_STATE=4



####################### `EvA_EvV_plot.py` #####################
ENERGY_OFFSET=4479.41619611 # Ry



##################### `OutputFileCreator.py` ##################
DISLOCATION_GRADE="full"




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
# echo "Based on $reference_structure, ibrav=$ibrav"
# lattice_parameter_bohr=$(echo "$LATTICE_PARAMETER*1.88973" | bc -l) # bohr
# input_filename=$ELEMENT_NAME # name input file


# ### write input file
# # input file
# echo """ &control
#     prefix=''
#     outdir='temp'
#     pseudo_dir = '../',
#  /
#  &system
#     ibrav=  $IBRAV, celldm(1) =$lattice_parameter_bohr, nat=  1, ntyp=  1,
#     ecutwfc =$CUTOFF_ENERGY,
#     occupations='smearing', smearing='mp', degauss=0.06
#  /
#  &electrons
#     electron_maxstep =$MAX_ITER,
#     mixing_beta =$MIXING_BETA,
#  /
# ATOMIC_SPECIES
#  $ELEMENT_NAME  $ELEMENT_AMU $PSEUDOPOTENTIAL_FILENAME
# ATOMIC_POSITIONS (alat)
#  $ELEMENT_NAME ${ELEMENT_POS[0]} ${ELEMENT_POS[1]} ${ELEMENT_POS[2]}
# K_POINTS (automatic)
#  ${KPOINTS[0]} ${KPOINTS[1]} ${KPOINTS[2]} ${KPOINTS_SHIFT[0]} ${KPOINTS_SHIFT[1]} ${KPOINTS_SHIFT[2]}""" > "$input_filename.in"
# echo "========== Contents of $input_filename.in =========="
# cat "$input_filename.in"
# echo "===================================================="


### automatically define other `ev_curve` file variables from inputs
sed -i "64s%^[[:digit:]]*[^ #]*%$EQ_OF_STATE%" "../0-Scripts/ev_curve"


### automatically define other `EvA_EvV_plot.py` file variables from inputs
sed -i "s%^energy_offset = [[:digit:]]*\.*[[:digit:]]*[^ #]%energy_offset = $ENERGY_OFFSET%" "../0-Scripts/EvA_EvV_plot.py"


### automatically define other `OutputFileCreator.py` file variables from inputs
dislocation_grade=$(echo $DISLOCATION_GRADE | tr '[:upper:]' '[:lower:]')
if [[ "$dislocation_grade" == "full" ]]; then
    disp_scale=$(echo "scale=9; sqrt(3) / 2" | bc)
elif [[ "$dislocation_grade" == "longpartial" ]]; then
    disp_scale=$(echo "scale=9; sqrt(2) / 2" | bc)
elif [[ "$dislocation_grade" == "partial" ]]; then
    disp_scale=$(echo "scale=9; 1" | bc)
else
    echo "Variable DISLOCATION_GRADE=$DISLOCATION_GRADE not understood. Must be either 'full', 'longpartial', or 'partial'."
    exit
fi

# modify script
sed -i "s%^num_proc = [[:digit:]]*[^ #]*%num_proc = 16%" "../0-Scripts/OutputFileCreator.py"
sed -i "s%^el = '[[:print:]]*'[^ #]*%el = '$ELEMENT_NAME'%" "../0-Scripts/OutputFileCreator.py"
sed -i "s%^potential = '[[:print:]]*'[^ #]*%potential = '$PSEUDOPOTENTIAL_FILENAME'%" "../0-Scripts/OutputFileCreator.py"
sed -i "s%^el_weight = [[:digit:]]*\.*[[:digit:]]*[^ #]*%el_weight = $ELEMENT_AMU%" "../0-Scripts/OutputFileCreator.py"
sed -i "s%^energy_cutoff = [[:digit:]]*\.*[[:digit:]]*\*13.6057[^ #]*%energy_cutoff = $CUTOFF_ENERGY\*13.6057%" "../0-Scripts/OutputFileCreator.py"
sed -i "s%^kpoints = [[:digit:]]*[^ #]*%kpoints = $KPOINT%" "../0-Scripts/OutputFileCreator.py"
sed -i "s%f.write(\"mixing_mode ='local-TF', electron_maxstep = [[:digit:]]*,\" + os.linesep)%f.write(\"mixing_mode ='local-TF', electron_maxstep = $MAX_ITER,\" + os.linesep)%" "../0-Scripts/OutputFileCreator.py"
sed -i "s%f.write(\"mixing_beta = [[:digit:]]*\.*[[:digit:]]*, conv_thr = 0.000001,\" + os.linesep)%f.write(\"mixing_beta = $MIXING_BETA, conv_thr = 0.000001,\" + os.linesep)%" "../0-Scripts/OutputFileCreator.py"

# modify script
sed -i "s%^num_proc = [[:digit:]]*[^ #]*%num_proc = 16%" "../0-Scripts/OutputSummarizer.py"
sed -i "s%^el = '[[:print:]]*'[^ #]*%el = '$ELEMENT_NAME'%" "../0-Scripts/OutputSummarizer.py"
sed -i "s%^potential = '[[:print:]]*'[^ #]*%potential = '$PSEUDOPOTENTIAL_FILENAME'%" "../0-Scripts/OutputSummarizer.py"
sed -i "s%^el_weight = [[:digit:]]*\.*[[:digit:]]*[^ #]*%el_weight = $ELEMENT_AMU%" "../0-Scripts/OutputSummarizer.py"
sed -i "s%^energy_cutoff = [[:digit:]]*\.*[[:digit:]]*\*13.6057[^ #]*%energy_cutoff = $CUTOFF_ENERGY\*13.6057%" "../0-Scripts/OutputSummarizer.py"
sed -i "s%^kpoints = [[:digit:]]*[^ #]*%kpoints = $KPOINT%" "../0-Scripts/OutputSummarizer.py"
sed -i "s%f.write(\"mixing_mode ='local-TF', electron_maxstep = [[:digit:]]*,\" + os.linesep)%f.write(\"mixing_mode ='local-TF', electron_maxstep = $MAX_ITER,\" + os.linesep)%" "../0-Scripts/OutputSummarizer.py"
sed -i "s%f.write(\"mixing_beta = [[:digit:]]*\.*[[:digit:]]*, conv_thr = 0.000001,\" + os.linesep)%f.write(\"mixing_beta = $MIXING_BETA, conv_thr = 0.000001,\" + os.linesep)%" "../0-Scripts/OutputSummarizer.py"



###################### GENERATE 3-GSFE DATA #####################
# ### execute QE with input parameters
# echo "Executing QE according to $input_filename.in..."
# (set -x; mpirun -np $NUM_PROC pw.x -in "$input_filename.in" > "$ELEMENT_NAME.out")


# ### run Fortran codes on input files
# gfortran -O2 "../0-Scripts/evfit.f" -o "../0-Scripts/evfit" # compiles `evfit.f` outputs `evfit`
# chmod +x "../0-Scripts/ev_curve" # makes file executable
# cp "$input_filename.in" "../0-Scripts/$reference_structure.ev.in" # create appropriate input file to `ev_curve`
# cd "../0-Scripts"
# # this outputs `evfit.4`: reference structure, lattice parameter
# ./ev_curve $reference_structure $LATTICE_PARAMETER
# python3 "EvA_EvV_plot.py" # generate plots
# mv "$reference_structure.ev.in" "../3-GSFE/$reference_structure.ev.in"
# mv "evfit" "../3-GSFE/evfit"
# mv "EvsA" "../3-GSFE/EvsA"
# mv "EvsV" "../3-GSFE/EvsV"
# mv "SUMMARY" "../3-GSFE/SUMMARY"
# mv "evfit.4" "../3-GSFE/evfit.4"
# mv "pw_ev.out" "../3-GSFE/pw_ev.out"
# mv "Name_of_EvA.pdf" "../3-GSFE/Name_of_EvA.pdf"
# mv "Name_of_EvV.pdf" "../3-GSFE/Name_of_EvV.pdf"
# mv "Name_of_Combined.pdf" "../3-GSFE/Name_of_Combined.pdf"
# rm -r "temp/"
# cd "../3-GSFE"
# declare -a inputs=(
#     "EvsA"
#     "EvsV"
# )
# for input in "${inputs[@]}"; do
#     echo "Opening $input to offset by $ENERGY_OFFSET eV..."
#     i=1
#     init_energy=$(sed -n "${i}s%^[[:digit:]]\.[[:digit:]]* %%p" "$input")
#     init_energy=$(echo "$init_energy+$ENERGY_OFFSET" | bc)
#     init_energy=$(sed -n "s%-*%%p" <<< $init_energy)
#     readarray file < $input
#     for line in "${file[@]}"; do
#         IFS=" " read -a elem <<< "$line"
#         energy=${elem[1]}
#         offset_energy=$(echo "$energy+$init_energy+$ENERGY_OFFSET" | bc)
#         sed -i "${i}s% \-[[:digit:]]*\.[[:digit:]]*% $offset_energy%" "$input"
#         i=$(echo "$i+1" | bc)
#     done # end line `i`
#     echo "Closing $input..."
# done # end of `input` file
# echo "Adjusting summary files by $ENERGY_OFFSET eV..."
# energy=$(sed -n "s%^Equilibrium Energy per Atom    = %%p" "SUMMARY")
# offset_energy=$(echo "$energy+$init_energy+$ENERGY_OFFSET" | bc)
# sed -i "s%^Equilibrium Energy per Atom    = \-[[:digit:]]*\.[[:digit:]]*%Equilibrium Energy per Atom    = $offset_energy%" "SUMMARY"
# lattice_parameter=$(sed -n "s%^Equilibrium lattice constant   = %%p" "SUMMARY")
# bulk=$(sed -n "s%^Bulk Modulus (kbar)            = %%p" "SUMMARY")
# bulk_gpa=$(echo "scale=9;$bulk/10" | bc)
# sed -i "s%^Bulk Modulus (kbar)            = [[:digit:]]*\.[[:digit:]]*%Bulk Modulus (GPa)             = $bulk_gpa%" "SUMMARY"

# cat "SUMMARY"

rm -r "RescaleUpload/"
mkdir "RescaleUpload"
cd "../0-Scripts/"
mkdir "input_gens"
# reference structure, lattice parameter, and block motion
python2 "OutputFileCreator.py" $reference_structure $LATTICE_PARAMETER $dislocation_grade
mv "input_gens/"* "../3-GSFE/RescaleUpload/"
rm -r "input_gens/"
rm "RE_comm.cm" "gsfe.in" "GSFE_SUMMARY"
cp "rescale_commands.sh" "../3-GSFE/RescaleUpload/rescale_commands.sh"