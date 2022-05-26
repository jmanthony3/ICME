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
echo "Setting up '1-ElectronicsToAtomistics' bridge..."
echo "+ $execution_dir/1-ElectronicsToAtomistics/setup.sh"
cd "$execution_dir/1-ElectronicsToAtomistics"
(set -x; ./setup.sh) \
    1>"$execution_dir/logs/1-ElectronicsToAtomistics/setup.log" \
    2>"$execution_dir/logs/1-ElectronicsToAtomistics/errors.log"

# bridge: 2
clear; echo "Setting up '2-AtomisticsToDislocationMobility' bridge..."
echo "+ $execution_dir/2-AtomisticsToDislocationMobility/setup.sh"
cd "$execution_dir/2-AtomisticsToDislocationMobility"
(set -x; ./setup.sh) \
    1>"$execution_dir/logs/2-AtomisticsToDislocationMobility/setup.log" \
    2>"$execution_dir/logs/2-AtomisticsToDislocationMobility/errors.log"

# bridge: 3
clear; echo "Setting up '3-DislocationMobilityToCrystalPlasticity' bridge..."
echo "+ $execution_dir/3-DislocationMobilityToCrystalPlasticity"
cd "$execution_dir/3-DislocationMobilityToCrystalPlasticity"
(set -x; ./setup.sh) \
    1>"$execution_dir/logs/3-DislocationMobilityToCrystalPlasticity/setup.log" \
    2>"$execution_dir/logs/3-DislocationMobilityToCrystalPlasticity/errors.log"

echo "Done."





# that's all folks