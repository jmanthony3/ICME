#!/bin/bash -x

# remove installed packages
sudo apt-get -y remove bc pv gfortran build-essential \
    quantum-espresso python3-pip python2 curl expect p7zip-full

# remove other install programs
rm -rf ~/.juliaup
rm -rf /usr/bin/.pw.x
rm -rf ~/QuantumEspresso
rm -rf ~/Ovito

# clean up registry
sudo apt-get -y autoremove

# clean out generated directories/files
## ./1-ElectronicsToAtomistics
rm -rf ./1-ElectronicsToAtomistics/Calculations/0-Scripts
rm -rf ./1-ElectronicsToAtomistics/Files/temp
rm -rf ./1-ElectronicsToAtomistics/Files/test
rm -rf ./1-ElectronicsToAtomistics/logs
