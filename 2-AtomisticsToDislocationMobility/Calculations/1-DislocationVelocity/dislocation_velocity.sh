#!/bin/bash -x




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - END OF HUMAN EDITABLE SECTION - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #




set +x
rm -r "RescaleUpload/"
mkdir "RescaleUpload"
mkdir "RescaleDownload"
mkdir "PositionFrameData"
cd "../0-Scripts/"
cp "Cu.meam" "../1-DislocationVelocity/RescaleUpload/"
cp "Dislocation.f90" "../1-DislocationVelocity/RescaleUpload/"
cp "DisVelocity.in" "../1-DislocationVelocity/RescaleUpload/"
cp "library.meam" "../1-DislocationVelocity/RescaleUpload/"
cp "rescale_commands.sh" "../1-DislocationVelocity/RescaleUpload/"
cp "atoms.sh" "../1-DislocationVelocity/RescaleUpload/"

(set -x; python3 -m pip install engineering_notation joby_m_anthony_iii pandas sympy)




# end of file