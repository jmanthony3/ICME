#!/bin/bash -x





# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - END OF HUMAN EDITABLE SECTION - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #





set +x
execution_dir=$(pwd)
who=$(whoami)
mkdir "$execution_dir/logs"
mkdir "$execution_dir/logs/1-ElectronicsToAtomistics"
mkdir "$execution_dir/logs/2-AtomisticsToDislocationMobility"
mkdir "$execution_dir/logs/3-DislocationMobilityToCrystalPlasticity"



########################### BRIDGES ###########################
### setup bridges
# bridge: 1
clear; echo "Setting up '1-ElectronicsToAtomistics' bridge..."
(set -x; $execution_dir/1-ElectronicsToAtomistics/setup.sh) \
    1>"$execution_dir/logs/1-ElectronicsToAtomistics/setup.log" \
    2>"$execution_dir/logs/1-ElectronicsToAtomistics/errors.log"

# bridge: 2
clear; echo "Setting up '2-AtomisticsToDislocationMobility' bridge..."
(set -x; $execution_dir/2-AtomisticsToDislocationMobility/setup.sh) \
    1>"$execution_dir/logs/2-AtomisticsToDislocationMobility/setup.log" \
    2>"$execution_dir/logs/2-AtomisticsToDislocationMobility/errors.log"

# bridge: 3
echo "Setting up '3-DislocationMobilityToCrystalPlasticity' bridge..."
cd "../3-DislocationMobilityToCrystalPlasticity"
(set -x; $execution_dir/3-DislocationMobilityToCrystalPlasticity/setup.sh) \
    1>"$execution_dir/logs/3-DislocationMobilityToCrystalPlasticity/setup.log" \
    2>"$execution_dir/logs/3-DislocationMobilityToCrystalPlasticity/errors.log"

echo "Done."





# that's all folks