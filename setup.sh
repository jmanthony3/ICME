#!/bin/bash -x




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - END OF HUMAN EDITABLE SECTION - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #




set +x
execution_dir=$(pwd)
who=$(whoami)
cd "1-ElectronicsToAtomistics/"
./setup.sh
cd "../2-AtomisticsToDislocationMobility/"
./setup.sh