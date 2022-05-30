#!/bin/bash -x



##################### `atom-dislocation` ######################
### define file variables and command line arguments
## reference structure
# 1) FCC
# 2) BCC
# 3) HCP
REFERENCE_STRUCTURE="bcc" # crystal structure

## element name
# 1) FCC
#       1) Ni
#       2) Cu
#       3) Al
# 1) BCC
#       1) Fe
#       2) Nb
#       3) Ta
#       4) Mo
#       5) W (Zhou)
#       6) W (Auckland)
# 1) HCP
#       1) Mg (Sun et al, 2006)
#       1) Mg (Liu et al, 1996)
#       1) Mg (experimental value)
ELEMENT_NAME="Fe" # Periodic Table identifier of element

ELEMENT_NUM="40 20 2" # number of `ELEMENT_NAME` atoms in x, y, z directions

## dislocation type
# edge
# screw
DISLOCATION_TYPE="Edge" # edge or screw

## structure type
# 1) cylinder
# 2) PAD
# 3) Perfect FCC crystal
STRUCTURE_TYPE="PAD"

# define according to: seq FIRST STEP LAST
declare -a TEMP=(300) # ($(seq 150 50 500)) # desired temperature [K]
# define according to: seq FIRST STEP LAST
declare -a SIGMA=($(seq 25 25 300)) # applied stress [MPa]
equilTime=10000 # number of increment to equilibrate the temperature
runTime=100000 # number of increment to calibrate the velocity





# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - END OF HUMAN EDITABLE SECTION - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #





set +x # turn script tracing off



######################### ATOMISTICS ##########################
ATOM_DISLOCATION_ARGUMENTS=( # collect arguments into single tuple
    $REFERENCE_STRUCTURE # 0
    $ELEMENT_NAME # 1
    $ELEMENT_NUM # 2
    $DISLOCATION_TYPE # 3
    $STRUCTURE_TYPE # 4
)
for temp in "${TEMP[@]}"; do
    mkdir "$temp" # make directory for temperature
    for sigma in "${SIGMA[@]}"; do
        # convert [MPa] tp [bar]
        sigma_bar=$(echo "scale=9; $sigma/1000*10" | bc)

        ### modify arguments tuple to replace strings with integer selections
        # reference structure selection
        reference_structure=$(echo $REFERENCE_STRUCTURE | tr '[:upper:]' '[:lower:]')
        if [[ "$reference_structure" == "fcc" ]]; then
            ATOM_DISLOCATION_ARGUMENTS[0]="1"
        elif [[ "$reference_structure" == "bcc" ]]; then
            ATOM_DISLOCATION_ARGUMENTS[0]="2"
        elif [[ "$reference_structure" == "hcp" ]]; then
            ATOM_DISLOCATION_ARGUMENTS[0]="3"
        else
            echo "Variable REFERENCE_STRUCTURE=$REFERENCE_STRUCTURE \
                not understood. Must be either 'fcc', 'bcc', or 'hcp'."
            exit
        fi
        # element name selection
        if [[ $ELEMENT_NAME == "Ni" ]] || [[ $ELEMENT_NAME == "Fe" ]]; then
            ATOM_DISLOCATION_ARGUMENTS[1]="1"
        elif [[ $ELEMENT_NAME == "Cu" ]] || [[ $ELEMENT_NAME == "Nb" ]]; then
            ATOM_DISLOCATION_ARGUMENTS[1]="2"
        elif [[ $ELEMENT_NAME == "Al" ]] || [[ $ELEMENT_NAME == "Ta" ]]; then
            ATOM_DISLOCATION_ARGUMENTS[1]="3"
        elif [[ $ELEMENT_NAME == "Mo" ]]; then
            ATOM_DISLOCATION_ARGUMENTS[1]="4"
        elif [[ $ELEMENT_NAME == "W (Zhou)" ]]; then
            ATOM_DISLOCATION_ARGUMENTS[1]="5"
        elif [[ $ELEMENT_NAME == "W (Auckland)" ]]; then
            ATOM_DISLOCATION_ARGUMENTS[1]="6"
        else
            echo "Variable ELEMENT_NAME=$ELEMENT_NAME not understood."
            exit
        fi
        # dislocation type selection
        dislocation_type=$(echo $DISLOCATION_TYPE | tr '[:upper:]' '[:lower:]')
        if [[ $dislocation_type == "edge" ]]; then
            ATOM_DISLOCATION_ARGUMENTS[3]="1"
        elif [[ $dislocation_type == "screw" ]]; then
            ATOM_DISLOCATION_ARGUMENTS[3]="2"
        else
            echo "Variable DISLOCATION_TYPE=$DISLOCATION_TYPE \
                not understood. Must be either 'edge' or 'screw'."
            exit
        fi
        # structure type selection
        structure_type=$(echo $STRUCTURE_TYPE | tr '[:upper:]' '[:lower:]')
        if [[ $structure_type == "cylinder" ]]; then
            ATOM_DISLOCATION_ARGUMENTS[4]="1"
        elif [[ $structure_type == "pad" ]]; then
            ATOM_DISLOCATION_ARGUMENTS[4]="2"
        elif [[ $structure_type == "perfect fcc crystal" ]]; then
            ATOM_DISLOCATION_ARGUMENTS[4]="3"
        else
            echo "Variable STRUCTURE_TYPE=$STRUCTURE_TYPE \
                not understood. Must be either 'cylinder', 'PAD' or 'Perfect FCC crystal'."
            exit
        fi
        # compile `atom-dislocation` script
        (set -x; gfortran -O3 "Dislocation.f90" -o "atom-dislocation")
        # generates atoms.`REFERENCE_STRUCTURE`.`DISLOCATION_TYPE`.`STRUCTURE_TYPE` file
        echo "Ignore consequent errors..."
        (set -x;
            ./atoms.sh ${ATOM_DISLOCATION_ARGUMENTS[0]} ${ATOM_DISLOCATION_ARGUMENTS[1]} $ELEMENT_NUM ${ATOM_DISLOCATION_ARGUMENTS[3]} ${ATOM_DISLOCATION_ARGUMENTS[4]}
        )

        ### modify `atom-dislocation` script
        # determine if encoded temperature is of type float or int
        if [[ $temp =~ ^[+-]?[0-9]*\.[0-9]+$ ]]
        then
            sed -i "s%^variable temp equal [[:digit:]]*\.*[[:digit:]]*[^[:blank:]#]*%variable temp equal $temp%" \
                "./DisVelocity.in"
        else
            sed -i "s%^variable temp equal [[:digit:]]*\.*[[:digit:]]*[^[:blank:]#]*%variable temp equal $temp.%" \
                "./DisVelocity.in"
        fi # should be float
        # determine if encoded stress is of type float or int
        if [[ $sigma_bar =~ ^[+-]?[0-9]*\.[0-9]+$ ]]
        then
            sed -i "s%^variable sigma equal [[:digit:]]*\.*[[:digit:]]*[^[:blank:]#]*%variable sigma equal $sigma_bar%" \
                "./DisVelocity.in"
        else
            sed -i "s%^variable sigma equal [[:digit:]]*\.*[[:digit:]]*[^[:blank:]#]*%variable sigma equal $sigma_bar.%" \
                "./DisVelocity.in"
        fi # should be float
        # update other variables
        sed -i "s%^variable material string [[:alpha:]]*[^[:blank:]#]*%variable material string $ELEMENT_NAME%" \
            "./DisVelocity.in"
        sed -i "s%^variable atom_file string atoms[\.[[:alpha:]]]*[^[:blank:]#]*%variable atom_file string atoms.$reference_structure.$dislocation_type.$structure_type%" \
            "./DisVelocity.in"
        sed -i "s%^variable equilTime equal [[:digit:]]*\.*[[:digit:]]*[^[:blank:]#]*%variable equilTime equal $equilTime%" \
            "./DisVelocity.in"
        sed -i "s%^variable runTime equal [[:digit:]]*\.*[[:digit:]]*[^[:blank:]#]*%variable runTime equal $runTime%" \
            "./DisVelocity.in"

        ### execute LAMMPS with input parameters
        nohup mpirun -np 50 lmp_rescale-amd-userlib -in "./DisVelocity.in"

        ### post-process LAMMPS
        mkdir "$temp/$sigma"
        mv "atom-dislocation" "./$temp/$sigma/"
        mv "atoms.$reference_structure.$dislocation_type.$structure_type" "./$temp/$sigma/"
        cp "DisVelocity.in" "./$temp/$sigma/"
        mv "dump.equilibration" "./$temp/$sigma/"
        mv "dump.minimize" "./$temp/$sigma/"
        mv "dump.shear" "./$temp/$sigma/"
        mv "dump.shear.unwrap" "./$temp/$sigma/"
        mv "log.lammps" "./$temp/$sigma/"
        # finish `sigma`
    done # finish `temp`
done # finish atomistics





# that's all folks