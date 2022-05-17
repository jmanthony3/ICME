#!/bin/bash -x




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - END OF HUMAN EDITABLE SECTION - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #




set +x
rm "inputs.txt"
touch "inputs.txt"
# process for each file
for input in gsfe_*.in; do
    echo "$input" >> "inputs.txt"
done
declare -a gsfes=()
for line in "inputs.txt"; do
    read -a elem <<< "$line"
    name=${elem[0]}
    gsfes+=($(sed -n "s%.in%%p" "$name"))
done
rm "inputs.txt"
for gsfe in ${gsfes[@]}; do
    mpirun -np 16 pw.x -in $gsfe.in > $gsfe.out
done




# end of file