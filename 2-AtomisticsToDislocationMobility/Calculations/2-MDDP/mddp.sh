#!/bin/bash -x



####################### INPUT VARIABLES #######################
### general
REFERENCE_STRUCTURE="bcc" # crystal structure
ELEMENT_NAME="Fe" # Periodic Table identifier of element


### `DDinput`
FRS_PAIRS="1" # number of Frank-Read Source (FRS) pairs
SOURCE_LENGTH="1 1000" # source length [b]
SEPARATION_DISTANCE="1 1000" # separation distance between FRS [b]
GLIDE_PLANE="20" # 1, 2, 3, 4, 5, 6, 20
RANDOM_SEED="123456" # 6-digit number in [100001-200001]
MAX_SEGMENT_LENGTH="100" # maximum segment length [b]
STRAIN_RATE="1000" # strain rate [1/s]
LOADING_DIRECTION="3" # loading direction, continuum notation
LOADING_CONDITION="0" # 1) creep; 0) constant
TIMESTEP="1.e-12" # [s]





# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - END OF HUMAN EDITABLE SECTION - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #





set +x # turn script tracing off
execution_dir=$(pwd) # where script executes from



#################### DISLOCATION MOBILITY #####################
### collect arguments into single tuple
FCC_ARGUMENTS=( # fcc
    "2" # 0
    "1" # 1
    $FRS_PAIRS # 2
    # $SOURCE_LENGTH # 3
    # $SEPARATION_DISTANCE # 4
    $GLIDE_PLANE # 5
    $RANDOM_SEED # 6
    "$(date +%Y)" # 7
    $MAX_SEGMENT_LENGTH # 8
    $STRAIN_RATE # 9
    $LOADING_DIRECTION # 10
    $LOADING_CONDITION # 11
    $TIMESTEP # 12
)
BCC_ARGUMENTS=( # bcc
    "2" # 0
    "1" # 1
    $FRS_PAIRS # 2
    # $SOURCE_LENGTH # 3
    # $SEPARATION_DISTANCE # 4
    $GLIDE_PLANE # 5
    $RANDOM_SEED # 6
    "$(date +%Y)" # 7
    $MAX_SEGMENT_LENGTH # 8
    $STRAIN_RATE # 9
    $LOADING_DIRECTION # 10
    $LOADING_CONDITION # 11
)


### generate `DDinput` file
reference_structure=$(echo $REFERENCE_STRUCTURE | tr '[:upper:]' '[:lower:]')
if [[ "$reference_structure" == "fcc" ]]; then
    if [ -f "datain" ]; then
        FCC_ARGUMENTS[0]="1" # allows script execution if present
    else
        FCC_ARGUMENTS[0]="2" # prohibits script execution if not
    fi
    clear
    # contains own procession procedure
    (set -x;
        ./fcc.sh ${FCC_ARGUMENTS[0]} ${FCC_ARGUMENTS[1]} ${FCC_ARGUMENTS[2]} \
            $SOURCE_LENGTH $SEPARATION_DISTANCE ${FCC_ARGUMENTS[3]} \
            ${FCC_ARGUMENTS[4]} ${FCC_ARGUMENTS[5]} ${FCC_ARGUMENTS[6]} \
            ${FCC_ARGUMENTS[7]} ${FCC_ARGUMENTS[8]} ${FCC_ARGUMENTS[9]} \
            ${FCC_ARGUMENTS[10]}
    )
elif [[ "$reference_structure" == "bcc" ]]; then
    if [ -f "datain" ]; then
        BCC_ARGUMENTS[0]="1" # allows script execution if present
    else
        BCC_ARGUMENTS[0]="2" # prohibits script execution if present
    fi
    ( # get missing library to execute `./BCCdata` binary
        wget http://archive.ubuntu.com/ubuntu/pool/universe/g/gcc-6/libgfortran3_6.4.0-17ubuntu1_amd64.deb
        sudo dpkg -i libgfortran3_6.4.0-17ubuntu1_amd64.deb # install library
        mv "libgfortran3_6.4.0-17ubuntu1_amd64.deb" "../../Files/"
    ) &> /dev/null
    clear
    # contains own procession procedure
    (set -x;
        ./bcc.sh ${BCC_ARGUMENTS[0]} ${BCC_ARGUMENTS[1]} ${BCC_ARGUMENTS[2]} \
            $SOURCE_LENGTH $SEPARATION_DISTANCE ${BCC_ARGUMENTS[3]} \
            ${BCC_ARGUMENTS[4]} ${BCC_ARGUMENTS[5]} ${BCC_ARGUMENTS[6]} \
            ${BCC_ARGUMENTS[7]} ${BCC_ARGUMENTS[8]} ${BCC_ARGUMENTS[9]}
    )
    sudo apt-get --fix-broken -y install &> /dev/null
else
    echo "Variable REFERENCE_STRUCTURE=$REFERENCE_STRUCTURE \
        not understood. Must be either 'fcc' or 'bcc'."
    exit
fi


### execute `./MDDP08-2008`
clear; (set -x; ./run.sh "$execution_dir")
sleep 10s # let previous process spin up
pid=$(pgrep MDDP08-2008) # get pid of `MDDP08-2008` process
echo "Kill the 'MDDP08-2008' process ('PID=$pid') whenever desired."
echo "Run the 'stress_strain.py' script at any time for live viewing."





# that's all folks