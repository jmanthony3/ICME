#!/bin/bash -x



########################### `*.sx` ############################
REFERENCE_STRUCTURE="bcc" # crystal structure





# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - END OF HUMAN EDITABLE SECTION - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #





set +x # turn script tracing off
execution_dir=$(pwd) # where script executes from
who=$(whoami) # current user




####################### SETUP ENVIRONMENT #####################
### test directory
# copy/paste necessary files for example
cd "$execution_dir/Files"
rm -r "test/" 2> /dev/null ; mkdir "test" 2> /dev/null
rm -r "test/RescaleUpload" 2> /dev/null ; mkdir "test/RescaleUpload" 2> /dev/null
# copy (in/out)put files to `./test/RescaleUpload/`
reference_structure=$(echo $REFERENCE_STRUCTURE | tr '[:upper:]' '[:lower:]')
if [[ "$reference_structure" == "fcc" ]]; then
    cp "fcc.sx" "./test/RescaleUpload/"
elif [[ "$reference_structure" == "bcc" ]]; then
    cp "bcc.sx" "./test/RescaleUpload/"
else
    echo "Variable REFERENCE_STRUCTURE=$REFERENCE_STRUCTURE \
        not understood. Must be either 'fcc' or 'bcc'."
    exit
fi
cp "numbers.inc" "./test/RescaleUpload/"
cp "params_xtal.inc" "./test/RescaleUpload/"
cp "rve.single.inp" "./test/RescaleUpload/"
cp "texture.txti" "./test/RescaleUpload/"
cp "umat_xtal.f" "./test/RescaleUpload/"
cp "test.xtali" "./test/RescaleUpload/"
cp "rescale_commands.sh" "./test/RescaleUpload/"


### populate "../Calculations/0-Scripts/" folder
mkdir "../Calculations/0-Scripts/" 2> /dev/null
cp "rescale_commands.sh" "../Calculations/0-Scripts/"





# that's all folks