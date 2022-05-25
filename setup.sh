#!/bin/bash -x




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - END OF HUMAN EDITABLE SECTION - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #




set +x
execution_dir=$(pwd)
who=$(whoami)



### setup bridges
# bride: 1
cd "1-ElectronicsToAtomistics/"
./setup.sh
# bride: 2
cd "../2-AtomisticsToDislocationMobility/"
./setup.sh
# bride: 3
cd "../3-DislocationMobilityToCrystalPlasticity/"
./setup.sh