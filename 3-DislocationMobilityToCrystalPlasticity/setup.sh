#!/bin/bash -x




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - END OF HUMAN EDITABLE SECTION - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #




set +x
execution_dir=$(pwd)
who=$(whoami)



# copy/paste necessary files for example
cd "$execution_dir/Files/"
rm -r "test/" || mkdir "test/"
rm -r "test/RescaleUpload/" || mkdir "test/RescaleUpload/"
cp "bcc.sx" "./test/RescaleUpload/"
cp "numbers.inc" "./test/RescaleUpload/"
cp "params_xtal.inc" "./test/RescaleUpload/"
cp "rve.single.inp" "./test/RescaleUpload/"
cp "texture.txti" "./test/RescaleUpload/"
cp "umat_xtal.f" "./test/RescaleUpload/"
cp "test.xtali" "./test/RescaleUpload/"
cp "rescale_commands.sh" "./test/RescaleUpload/"

### populate "../Calculations/0-Scripts/" folder
cp "rescale_commands.sh" "../Calculations/0-Scripts/"





# that's all folks