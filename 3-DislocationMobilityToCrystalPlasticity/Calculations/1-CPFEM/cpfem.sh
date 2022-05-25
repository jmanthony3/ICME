#!/bin/bash -x




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - END OF HUMAN EDITABLE SECTION - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #




set +x
execution_dir=$(pwd)
who=$(whoami)



# copy/paste necessary files for example
rm -r "RescaleUpload/" || mkdir "RescaleUpload/"
cp "../../Files/bcc.sx" "./RescaleUpload/"
cp "../../Files/numbers.inc" "./RescaleUpload/"
cp "../../Files/params_xtal.inc" "./RescaleUpload/"
cp "../../Files/rve.single.inp" "./RescaleUpload/"
cp "../../Files/texture.txti" "./RescaleUpload/"
cp "../../Files/umat_xtal.f" "./RescaleUpload/"
cp "../../Files/test.xtali" "./RescaleUpload/"
cp "../0-Scripts/rescale_commands.sh" "./RescaleUpload/"





# that's all folks