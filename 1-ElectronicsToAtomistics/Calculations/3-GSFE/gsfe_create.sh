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
CUTOFF_ENERGY=64 # [Ry]
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
# working language to perform calculations and plot results
COMPUTING_LANGUAGE="Julia" # can also be "Python"
NUM_PROC=$(nproc)   # grabs all cores available by default
UNDER_PROC=2        # number of cores to under-subscribe, if possible



########################## `ev_curve` #########################
# input arguments to script
# also uses `REFERENCE_STRUCTURE` and `LATTICE_PARAMETER`
# 1=Birch(First), 2=Birch(Second), 3=Keane, 4=Murnaghan
EQ_OF_STATE=4



####################### `EvA_EvV_plot.py` #####################
ENERGY_OFFSET=4479.47132505 # [Ry]



###################### `OutputFileCreator` ####################
DISLOCATION_GRADE="full"





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
num_proc=$(echo "$NUM_PROC - $UNDER_PROC" | bc)
if [[ $(nproc) -eq 1 ]]; then
    num_proc=1
elif [[ $num_proc -lt 1 ]]; then
    echo "(NUM_PROC=$NUM_PROC) - (UNDER_PROC=$UNDER_PROC) = \
        $num_proc which must be greater than 1."
    exit
# elif [[ $NUM_PROC - $UNDER_PROC -le $(nproc) ]]; then
    # num_proc=$NUM_PROC - $UNDER_PROC
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
input_filename=$ELEMENT_NAME # name input file


### write input file
# input file
echo """ &control
    prefix=''
    outdir='temp'
    pseudo_dir = '../',
 /
 &system
    ibrav=  $ibrav, celldm(1) =$lattice_parameter_bohr, nat=  1, ntyp=  1,
    ecutwfc =$CUTOFF_ENERGY,
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
 ${KPOINTS[0]} ${KPOINTS[1]} ${KPOINTS[2]} ${KPOINTS_SHIFT[0]} ${KPOINTS_SHIFT[1]} ${KPOINTS_SHIFT[2]}\
""" > "$input_filename.in"
# print input file to terminal
echo "========== Contents of $input_filename.in =========="
cat "$input_filename.in"
echo "========================================"


### automatically define other `ev_curve` file variables from inputs
sed -i "64s%^[[:digit:]]*[^ #]%$EQ_OF_STATE%" "../0-Scripts/ev_curve"


### automatically define other `EvA_EvV_plot` file variables from inputs
sed -i "s%^energy_offset = [[:digit:]]*\.*[[:digit:]]*[^ #]%energy_offset = $ENERGY_OFFSET%" \
    "../0-Scripts/EvA_EvV_plot.$cl_ext"


### automatically define other `OutputFileCreator` file variables from inputs
dislocation_grade=$(echo $DISLOCATION_GRADE | tr '[:upper:]' '[:lower:]')
if [[ "$dislocation_grade" == "full" ]]; then
    disp_scale=$(echo "scale=9; sqrt(3) / 2" | bc)
elif [[ "$dislocation_grade" == "longpartial" ]]; then
    disp_scale=$(echo "scale=9; sqrt(2) / 2" | bc)
elif [[ "$dislocation_grade" == "partial" ]]; then
    disp_scale=$(echo "scale=9; 1" | bc)
else
    echo "Variable DISLOCATION_GRADE=$DISLOCATION_GRADE \
        not understood. Must be either 'full', 'longpartial', or 'partial'."
    exit
fi

# modify `OutputFileCreator` script
sed -i "s%^num_proc = [[:digit:]]*[^ #]%num_proc = 16%" \
    "../0-Scripts/OutputFileCreator.$cl_ext"
sed -i "s%^el = '[[:print:]]*'[^ #]%el = '$ELEMENT_NAME'%" \
    "../0-Scripts/OutputFileCreator.$cl_ext"
sed -i "s%^potential = '[[:print:]]*'[^ #]%potential = '$PSEUDOPOTENTIAL_FILENAME'%" \
    "../0-Scripts/OutputFileCreator.$cl_ext"
sed -i "s%^el_weight = [[:digit:]]*\.*[[:digit:]]*[^ #]%el_weight = $ELEMENT_AMU%" \
    "../0-Scripts/OutputFileCreator.$cl_ext"
sed -i "s%^energy_cutoff = [[:digit:]]*\.*[[:digit:]]*\*13.6057[^ #]%energy_cutoff = $CUTOFF_ENERGY\*13.6057%" \
    "../0-Scripts/OutputFileCreator.$cl_ext"
sed -i "s%^kpoints = [[:digit:]]*[^ #]%kpoints = $KPOINT%" \
    "../0-Scripts/OutputFileCreator.$cl_ext"
sed -i "s%f.write(\"mixing_mode ='local-TF', electron_maxstep = [[:digit:]]*,\" + os.linesep)%f.write(\"mixing_mode ='local-TF', electron_maxstep = $MAX_ITER,\" + os.linesep)%" \
    "../0-Scripts/OutputFileCreator.$cl_ext"
sed -i "s%f.write(\"mixing_beta = [[:digit:]]*\.*[[:digit:]]*, conv_thr = 0.000001,\" + os.linesep)%f.write(\"mixing_beta = $MIXING_BETA, conv_thr = 0.000001,\" + os.linesep)%" \
    "../0-Scripts/OutputFileCreator.$cl_ext"

# modify `OutputFileSummarizer` script
sed -i "s%^num_proc = [[:digit:]]*[^ #]%num_proc = 16%" \
    "../0-Scripts/OutputFileSummarizer.$cl_ext"
sed -i "s%^el = '[[:print:]]*'[^ #]%el = '$ELEMENT_NAME'%" \
    "../0-Scripts/OutputFileSummarizer.$cl_ext"
sed -i "s%^potential = '[[:print:]]*'[^ #]%potential = '$PSEUDOPOTENTIAL_FILENAME'%" \
    "../0-Scripts/OutputFileSummarizer.$cl_ext"
sed -i "s%^el_weight = [[:digit:]]*\.*[[:digit:]]*[^ #]%el_weight = $ELEMENT_AMU%" \
    "../0-Scripts/OutputFileSummarizer.$cl_ext"
sed -i "s%^energy_cutoff = [[:digit:]]*\.*[[:digit:]]*\*13.6057[^ #]%energy_cutoff = $CUTOFF_ENERGY\*13.6057%" \
    "../0-Scripts/OutputFileSummarizer.$cl_ext"
sed -i "s%^kpoints = [[:digit:]]*[^ #]%kpoints = $KPOINT%" \
    "../0-Scripts/OutputFileSummarizer.$cl_ext"
sed -i "s%f.write(\"mixing_mode ='local-TF', electron_maxstep = [[:digit:]]*,\" + os.linesep)%f.write(\"mixing_mode ='local-TF', electron_maxstep = $MAX_ITER,\" + os.linesep)%" \
    "../0-Scripts/OutputFileSummarizer.$cl_ext"
sed -i "s%f.write(\"mixing_beta = [[:digit:]]*\.*[[:digit:]]*, conv_thr = 0.000001,\" + os.linesep)%f.write(\"mixing_beta = $MIXING_BETA, conv_thr = 0.000001,\" + os.linesep)%" \
    "../0-Scripts/OutputFileSummarizer.$cl_ext"



###################### GENERATE GSFE DATA #####################
### execute QE with input parameters
echo "Executing QE according to $input_filename.in..."
(set -x;
    mpirun -np $num_proc --use-hwthread-cpus pw.x -in "$input_filename.in" \
        > "$input_filename.out" 2> /dev/null
)
rm -r "temp/" # remove calculations temporary folder


### run Fortran codes on input files
# create appropriate input file to `ev_curve`
cp "$input_filename.in" "../0-Scripts/$reference_structure.ev.in"

# move to Scripts folder
cd "../0-Scripts"
gfortran -std=legacy -O2 "evfit.f" -o "evfit" 2> /dev/null # compiles `evfit.f` outputs `evfit`
# this outputs `evfit.4`: reference structure, lattice parameter
./ev_curve $reference_structure $LATTICE_PARAMETER 2> /dev/null
if [[ "$computing_language" == "julia" ]]; then
    julia "EvA_EvV_plot.jl" # generate plots
elif [[ "$computing_language" == "python" ]]; then
    python3 "EvA_EvV_plot.py" # generate plots
else
    echo "Variable COMPUTING_LANGUAGE=$COMPUTING_LANGUAGE \
        not understood. Must be either 'Julia' or 'Python'."
    exit
fi

# move all output files back to working directory
mv "$reference_structure.ev.in" "../3-GSFE/"
mv "evfit" "../3-GSFE/"
mv "EvsA" "../3-GSFE/"
mv "EvsV" "../3-GSFE/"
mv "SUMMARY" "../3-GSFE/"
mv "evfit.4" "../3-GSFE/"
mv "pw_ev.out" "../3-GSFE/"
mv "Name_of_EvA.png" "../3-GSFE/"
mv "Name_of_EvV.png" "../3-GSFE/"
mv "Name_of_Combined.png" "../3-GSFE/"
rm -r "temp/" # remove calculations temporary folder


### post-process output files
# move back to working directory
cd "../3-GSFE"
declare -a inputs=(
    "EvsA"
    "EvsV"
)
for input in "${inputs[@]}"; do
    echo -e -n "Opening $input to offset by $ENERGY_OFFSET eV...\r"
    i=1
    # init_energy=$(sed -n "${i}s%^[[:digit:]]\.[[:digit:]]* %%p" "$input")
    # init_energy=$(echo "$init_energy+$ENERGY_OFFSET" | bc)
    # init_energy=$(sed -n "s%-*%%p" <<< $init_energy)
    init_energy=0
    readarray file < $input
    for line in "${file[@]}"; do
        IFS=" " read -a elem <<< "$line"
        energy=${elem[1]}
        offset_energy=$(echo "$energy+$init_energy+$ENERGY_OFFSET" | bc)
        sed -i "${i}s% \-[[:digit:]]*\.[[:digit:]]*% $offset_energy%" "$input"
        i=$(echo "$i+1" | bc) # end of line `i`
    done # end line of `input` file
    echo -e -n "Closing $input...\r"
done # end of processing

# adjust summary file
echo -e -n "Adjusting summary files by $ENERGY_OFFSET eV...\r"
# previous energy value
energy=$(sed -n "s%^Equilibrium Energy per Atom    = %%p" "SUMMARY")
# offset energy value
offset_energy=$(echo "$energy+$init_energy+$ENERGY_OFFSET" | bc)
# replace energy value
sed -i "s%^Equilibrium Energy per Atom    = \-[[:digit:]]*\.[[:digit:]]*%Equilibrium Energy per Atom    = $offset_energy%" \
    "SUMMARY"
# get lattice parameter [angstrom]
lattice_parameter=$(sed -n "s%^Equilibrium lattice constant   = %%p" "SUMMARY")
# replace in `2-AtomisticsToDislocationMobility/1-DislocationVelocity/dislocation_velocity`
# 1 [angstrom] = 1e-10 [m]
sed -i "s%^lattice_parameter = [[:digit:]\.e\-]*[^ #]%lattice_parameter = ${lattice_parameter}e\-10%" \
    "../../../2-AtomisticsToDislocationMobility/Calculations/1-DislocationVelocity/dislocation_velocity.$cl_ext"
# replace in `2-AtomisticsToDislocationMobility/2-MDDP/stress_strain`
# 1 [angstrom] = 1e-10 [m]
sed -i "s%^lattice_parameter = [[:digit:]\.e\-]*[^ #]%lattice_parameter = ${lattice_parameter}e\-10%" \
    "../../../2-AtomisticsToDislocationMobility/Calculations/2-MDDP/stress_strain.$cl_ext"
# previous bulk modulus [kbar]
bulk=$(sed -n "s%^Bulk Modulus (kbar)            = %%p" "SUMMARY")
# convert to [GPa]
bulk_gpa=$(echo "scale=9;$bulk/10" | bc)
# replace bulk modulus value
sed -i "s%^Bulk Modulus (kbar)            = [[:digit:]]*\.[[:digit:]]*%Bulk Modulus (GPa)             = $bulk_gpa%" \
    "SUMMARY"

# print study summary to terminal
echo "                                                        "
cat "SUMMARY"
echo -e -n "\nEnergy offset found to be $ENERGY_OFFSET eV\n"


### create Rescale(Up/Down)load structure
# make `RescaleUpload/` folder from inputs
rm -r "RescaleUpload/" 2> /dev/null ; mkdir "RescaleUpload"
# make `RescaleDownload/` folder for outputs
mkdir "RescaleDownload" 2> /dev/null
# move into Scripts directory
cd "../0-Scripts/"
mkdir "input_gens" # subsequent script places inputs here
echo "                                                        "
echo "                                                        "
if [[ "$computing_language" == "julia" ]]; then
    (set -x;
        # julia OutputFileCreator.jl: reference structure, lattice parameter, block motion
        julia "OutputFileCreator.jl" $reference_structure $LATTICE_PARAMETER $dislocation_grade
    )
elif [[ "$computing_language" == "python" ]]; then
    (set -x;
        # python2 OutputFileCreator.py: reference structure, lattice parameter, block motion
        python2 "OutputFileCreator.py" $reference_structure $LATTICE_PARAMETER $dislocation_grade
    )
else
    echo "Variable COMPUTING_LANGUAGE=$COMPUTING_LANGUAGE \
        not understood. Must be either 'Julia' or 'Python'."
    exit
fi
# move (in/out)put files to `../3-GSFE/RescaleUpload/`
mv "input_gens/"* "../3-GSFE/RescaleUpload/"
rm -r "input_gens/"
rm "RE_comm.cm" "gsfe.in" "GSFE_SUMMARY"
cp "rescale_commands.sh" "../3-GSFE/RescaleUpload/"





# that's all folks