#!/bin/bash -x





# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - END OF HUMAN EDITABLE SECTION - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #





set +x # turn script tracing off
execution_dir=$(pwd)  # where script executes from
who=$(whoami) # current user
mkdir "$execution_dir/logs"
read -sp 'Password for sudo commands: ' pswd



########################### BRIDGES ###########################
### setup bridges
# make sure distro is current
(set -x; (echo $pswd) 2> /dev/null | sudo -S apt-get update)

# bridge: 1
echo "Setting up '1-ElectronicsToAtomistics' bridge..."
echo "+ $execution_dir/1-ElectronicsToAtomistics/setup.sh"
cd "$execution_dir/1-ElectronicsToAtomistics"
./setup.sh $pswd | tee "$execution_dir/logs/1-ElectronicsToAtomistics.log"

# bridge: 2
clear; echo "Setting up '2-AtomisticsToDislocationMobility' bridge..."
echo "+ $execution_dir/2-AtomisticsToDislocationMobility/setup.sh"
cd "$execution_dir/2-AtomisticsToDislocationMobility"
./setup.sh $pswd | tee "$execution_dir/logs/2-AtomisticsToDislocationMobility.log"

# bridge: 3
clear; echo "Setting up '3-DislocationMobilityToCrystalPlasticity' bridge..."
echo "+ $execution_dir/3-DislocationMobilityToCrystalPlasticity"
cd "$execution_dir/3-DislocationMobilityToCrystalPlasticity"
./setup.sh $pswd | tee "$execution_dir/logs/3-DislocationMobilityToCrystalPlasticity.log"

clear; echo "Done."





# that's all folks