#!/bin/bash -x



######################### `atom-dislocation` ##########################
### define file variables and command line arguments
# reference structure
# 1) FCC 2) BCC 3) HCP
REFERENCE_STRUCTURE="fcc" # crystal structure

# element name
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
ELEMENT_NAME="Cu" # Periodic Table identifier of element

ELEMENT_NUM="40 20 2" # number of `ELEMENT_NAME` atoms in x, y, z directions

# dislocation type
DISLOCATION_TYPE="Edge" # edge or screw

# structure type
STRUCTURE_TYPE="PAD" # 1) cylinder; 2) PAD; 3) Perfect FCC crystal

# define according to: seq FIRST STEP LAST
declare -a TEMP=($(seq 150 50 500))        # desired temperature
# define according to: seq FIRST STEP LAST
declare -a SIGMA=($(seq 25 25 300))          # applied stress in MPa
equilTime=10000     # number of increment to equilibrate the temperature
runTime=100000      # number of increment to calibrate the velocity




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - END OF HUMAN EDITABLE SECTION - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #




### collect arguments into single tuple
ATOM_DISLOCATION_ARGUMENTS=(
    $REFERENCE_STRUCTURE # 0
    $ELEMENT_NAME # 1
    $ELEMENT_NUM # 2
    $DISLOCATION_TYPE # 3
    $STRUCTURE_TYPE # 4
)
for temp in "${TEMP[@]}"; do
    mkdir "./$temp"
    for sigma in "${SIGMA[@]}"; do
        sigma_bar=$(echo "scale=9; $sigma/1000*10" | bc)


        ### modify arguments tuple to replace strings with integer selections
        # reference structure selection
        case $(echo $REFERENCE_STRUCTURE | tr '[:upper:]' '[:lower:]') in

            "fcc")
                ATOM_DISLOCATION_ARGUMENTS[0]="1"
                ;;

            "bcc")
                ATOM_DISLOCATION_ARGUMENTS[0]="2"
                ;;

            "hcp")
                ATOM_DISLOCATION_ARGUMENTS[0]="3"
                ;;
        esac
        # element name selection
        case $ELEMENT_NAME in

            "Ni" | "Fe")
                ATOM_DISLOCATION_ARGUMENTS[1]="1"
                ;;

            "Cu" | "Nb")
                ATOM_DISLOCATION_ARGUMENTS[1]="2"
                ;;

            "Al" | "Ta")
                ATOM_DISLOCATION_ARGUMENTS[1]="3"
                ;;

            "Mo")
                ATOM_DISLOCATION_ARGUMENTS[1]="4"
                ;;

            "W (Zhou)")
                ATOM_DISLOCATION_ARGUMENTS[1]="5"
                ;;

            "W (Auckland)")
                ATOM_DISLOCATION_ARGUMENTS[1]="6"
                ;;
        esac
        case $(echo $DISLOCATION_TYPE | tr '[:upper:]' '[:lower:]') in

            "edge")
                ATOM_DISLOCATION_ARGUMENTS[3]="1"
                ;;

            "screw")
                ATOM_DISLOCATION_ARGUMENTS[3]="2"
                ;;
        esac
        case $(echo $STRUCTURE_TYPE | tr '[:upper:]' '[:lower:]') in

            "cylinder")
                ATOM_DISLOCATION_ARGUMENTS[4]="1"
                ;;

            "pad")
                ATOM_DISLOCATION_ARGUMENTS[4]="2"
                ;;

            "perfect fcc crystal")
                ATOM_DISLOCATION_ARGUMENTS[4]="3"
                ;;
        esac


        gfortran -O3 "Dislocation.f90" -o "atom-dislocation"
        # generates atom.`REFERENCE_STRUCTURE`.`DISLOCATION_TYPE`.`STRUCTURE_TYPE` file
        ./atoms.sh ${ATOM_DISLOCATION_ARGUMENTS[0]} ${ATOM_DISLOCATION_ARGUMENTS[1]} $ELEMENT_NUM ${ATOM_DISLOCATION_ARGUMENTS[3]} ${ATOM_DISLOCATION_ARGUMENTS[4]}


        ### ignore consequent errors


        if [[ $temp =~ ^[+-]?[0-9]*\.[0-9]+$ ]]
        then
            sed -i "s%^variable temp equal [[:digit:]]*\.*[[:digit:]]*[^[:blank:]#]*%variable temp equal $temp%" "./DisVelocity.in"
        else
            sed -i "s%^variable temp equal [[:digit:]]*\.*[[:digit:]]*[^[:blank:]#]*%variable temp equal $temp.%" "./DisVelocity.in"
        fi

        if [[ $sigma_bar =~ ^[+-]?[0-9]*\.[0-9]+$ ]]
        then
            sed -i "s%^variable sigma equal [[:digit:]]*\.*[[:digit:]]*[^[:blank:]#]*%variable sigma equal $sigma_bar%" "./DisVelocity.in"
        else
            sed -i "s%^variable sigma equal [[:digit:]]*\.*[[:digit:]]*[^[:blank:]#]*%variable sigma equal $sigma_bar.%" "./DisVelocity.in"
        fi

        sed -i "s%^variable material string [[:alpha:]]*[^[:blank:]#]*%variable material string $ELEMENT_NAME%" "./DisVelocity.in"

        sed -i "s%^variable atom_file string atoms[\.[[:alpha:]]]*[^[:blank:]#]*%variable atom_file string atoms.$(echo $REFERENCE_STRUCTURE | tr '[:upper:]' '[:lower:]').$(echo $DISLOCATION_TYPE | tr '[:upper:]' '[:lower:]').$(echo $STRUCTURE_TYPE | tr '[:upper:]' '[:lower:]')%" "./DisVelocity.in"

        sed -i "s%^variable equilTime equal [[:digit:]]*\.*[[:digit:]]*[^[:blank:]#]*%variable equilTime equal $equilTime%" "./DisVelocity.in"

        sed -i "s%^variable runTime equal [[:digit:]]*\.*[[:digit:]]*[^[:blank:]#]*%variable runTime equal $runTime%" "./DisVelocity.in"



        # run LAMMPS
        nohup mpirun -np 50 lmp_rescale-amd-userlib -in "./DisVelocity.in"

        # post-process LAMMPS
        mkdir "./$temp/$sigma"
        mv "./atom-dislocation" "./$temp/$sigma/atom-dislocation"
        mv "./atoms.$(echo $REFERENCE_STRUCTURE | tr '[:upper:]' '[:lower:]').$(echo $DISLOCATION_TYPE | tr '[:upper:]' '[:lower:]').$(echo $STRUCTURE_TYPE | tr '[:upper:]' '[:lower:]')" "./$sigma/atoms.$(echo $REFERENCE_STRUCTURE | tr '[:upper:]' '[:lower:]').$(echo $DISLOCATION_TYPE | tr '[:upper:]' '[:lower:]').$(echo $STRUCTURE_TYPE | tr '[:upper:]' '[:lower:]')"
        cp "./DisVelocity.in" "./$temp/$sigma/DisVelocity.in"
        mv "./dump.equilibration" "./$temp/$sigma/dump.equilibration"
        mv "./dump.minimize" "./$temp/$sigma/dump.minimize"
        mv "./dump.shear" "./$temp/$sigma/dump.shear"
        mv "./dump.shear.unwrap" "./$temp/$sigma/dump.shear.unwrap"
        mv "./log.lammps" "./$temp/$sigma/log.lammps"
    done # finish `sigma`
done # finish `temp`




# end of file
