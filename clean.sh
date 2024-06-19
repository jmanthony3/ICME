#!/bin/bash -x

execution_dir=$(pwd)  # where script executes from

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
cd ./1-ElectronicsToAtomistics/Calculations/1-EnergyOffset
rm -rf !("energy_offset.sh")
cd ../2-KPointStudy
rm -rf !("kpoint_study.sh"|"which_ecutwfc.jl"|"which_ecutwfc.py")
cd ../3-GSFE
rm -rf !("gsfe_create.sh"|"gsfe_process.sh")
cd $execution_dir
rm -rf ./1-ElectronicsToAtomistics/Files/temp
rm -rf ./1-ElectronicsToAtomistics/Files/test
rm -rf ./1-ElectronicsToAtomistics/logs

## ./2-AtomisticsToDislocationMobility
rm -rf ./2-AtomisticsToDislocationMobility/Calculations/0-Scripts
cd ./2-AtomisticsToDislocationMobility/Calculations/1-DislocationVelocity
rm -rf !("dislocation_velocity.sh"|"dislocation_position.jl"|"dislocation_position.py"|"dislocation_velocity.jl"|"dislocation_velocity.py")
cd ../2-MDDP
rm -rf !("bcc.sh"|"fcc.sh"|"mddp.sh"|"run.sh"|"stress_strain.jl"|"stress_strain.py")
cd $execution_dir
rm -rf ./2-AtomisticsToDislocationMobility/Files/MDDP
rm -rf ./2-AtomisticsToDislocationMobility/Files/test
rm -rf ./2-AtomisticsToDislocationMobility/logs

## ./3-DislocationMobilityToCrystalPlasticity
rm -rf ./3-DislocationMobilityToCrystalPlasticity/Calculations/0-Scripts
cd ./3-DislocationMobilityToCrystalPlasticity/Calculations/1-CPFEM
rm -rf !("cpfem.sh"|"voce_hardening_fit.jl"|"voce_hardening_fit.py")
cd $execution_dir
rm -rf ./3-DislocationMobilityToCrystalPlasticity/Files/test
