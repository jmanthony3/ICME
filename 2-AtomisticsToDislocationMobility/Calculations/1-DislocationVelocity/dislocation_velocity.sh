#!/bin/bash -x




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - END OF HUMAN EDITABLE SECTION - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #




set +x
### collect arguments into single tuple
rm -r "RescaleUpload/"
mkdir "RescaleUpload"
cd "../0-Scripts/"
cp "Cu.meam" "../1-DislocationVelocity/RescaleUpload/"
cp "Dislocation.f90" "../1-DislocationVelocity/RescaleUpload/"
cp "DisVelocity.in" "../1-DislocationVelocity/RescaleUpload/"
cp "library.meam" "../1-DislocationVelocity/RescaleUpload/"
cp "rescale_commands.sh" "../1-DislocationVelocity/RescaleUpload/"
cp "atoms.sh" "../1-DislocationVelocity/RescaleUpload/"




# end of file