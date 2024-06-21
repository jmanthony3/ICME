#!/bin/bash -x



####################### MEAM PARAMETERS #######################
ELEMENT_NAME="Fe" # Periodic Table identifier of element





# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - END OF HUMAN EDITABLE SECTION - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #





set +x # turn script tracing off



###################### CREATE ENVIRONMENT #####################
# rm -r "RescaleUpload/" 2> /dev/null
mkdir "RescaleUpload" 2> /dev/null
(mkdir "RescaleDownload"; mkdir "PositionFrameData") 2> /dev/null
# copy input files to `RescaleUpload/`
cp "../../../1-ElectronicsToAtomistics/Calculations/$ELEMENT_NAME.parameter.meam" "./RescaleUpload/$ELEMENT_NAME.meam"
cp "../../../1-ElectronicsToAtomistics/Calculations/$ELEMENT_NAME.library.meam" "./RescaleUpload/library.meam"
cp "../../Files/Dislocation.f90" "./RescaleUpload/"
cp "../../Files/DisVelocity.in" "./RescaleUpload/"
cp "../0-Scripts/rescale_commands.sh" "./RescaleUpload/"
cp "../0-Scripts/atoms.sh" "./RescaleUpload/"





# that's all folks