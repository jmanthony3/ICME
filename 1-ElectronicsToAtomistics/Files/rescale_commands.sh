#!/bin/bash -x





# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - END OF HUMAN EDITABLE SECTION - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #





set +x # turn script tracing off



##################### CALCULATE GSFE CURVE ####################
### glob list of `gsfe_*.in` files
touch "inputs.txt"
for input in gsfe_*.in; do # process for each file
    echo "$input" >> "inputs.txt" # append list
done # end list definition


### glob numbers from list of `gsfe_*.in` files
declare -a gsfes=() # start with empty tuple
for line in "inputs.txt"; do
    read -a elem <<< "$line"
    name=${elem[0]}
    gsfes+=($(sed -n "s%.in%%p" "$name")) # append tuple
done # end of tuple definition
rm "inputs.txt" # remove temporary list


### execute QE with input parameters
for gsfe in ${gsfes[@]}; do
    (set -x; mpirun -np 16 pw.x -in $gsfe.in > $gsfe.out)
done





# that's all folks