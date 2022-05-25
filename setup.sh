#!/bin/bash -x





# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - END OF HUMAN EDITABLE SECTION - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #





set +x
execution_dir=$(pwd)
who=$(whoami)



########################### BRIDGES ###########################
### setup bridges
echo "Setting up '1-ElectronicsToAtomistics' bridge..."
cd "1-ElectronicsToAtomistics"; ./setup.sh # 1
clear; echo "Setting up '2-AtomisticsToDislocationMobility' bridge..."
cd "../2-AtomisticsToDislocationMobility"; ./setup.sh # 2
clear; echo "Setting up '3-DislocationMobilityToCrystalPlasticity' bridge..."
cd "../3-DislocationMobilityToCrystalPlasticity"; ./setup.sh # 3
echo "Done."





# that's all folks