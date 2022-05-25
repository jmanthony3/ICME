#!/bin/bash -x





# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - END OF HUMAN EDITABLE SECTION - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #





set +x
execution_dir=$(pwd)
who=$(whoami)



########################### BRIDGES ###########################
### setup bridges
cd "1-ElectronicsToAtomistics"; ./setup.sh # 1
cd "../2-AtomisticsToDislocationMobility"; ./setup.sh # 2
cd "../3-DislocationMobilityToCrystalPlasticity"; ./setup.sh # 3





# that's all folks