#!/bin/bash -x



####################### MEAM PARAMETERS #######################
ELEMENT_NAME="Fe" # Periodic Table identifier of element





# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - END OF HUMAN EDITABLE SECTION - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #





set +x # turn script tracing off



###################### CREATE ENVIRONMENT #####################
rm -r "RescaleUpload/" || mkdir "RescaleUpload"
mkdir "RescaleDownload"; mkdir "PositionFrameData"
# copy input files to `RescaleUpload/`
cp "../0-Scripts/$ELEMENT_NAME.meam" "./RescaleUpload/"
cp "../0-Scripts/Dislocation.f90" "./RescaleUpload/"
cp "../0-Scripts/DisVelocity.in" "./RescaleUpload/"
cp "../0-Scripts/library.meam" "./RescaleUpload/"
cp "../0-Scripts/rescale_commands.sh" "./RescaleUpload/"
cp "../0-Scripts/atoms.sh" "./RescaleUpload/"





# that's all folks