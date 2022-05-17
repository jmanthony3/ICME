#!/bin/bash -x



# define path into where Quantum Espresso will install
OVITO_INSTALL_LOC=~/Ovito
# encode name and version of of tarball
OVITO_VERSION="ovito-basic-3.7.2-x86_64"
NUM_PROC=$(nproc) # grabs all cores available by default




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - END OF HUMAN EDITABLE SECTION - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #




set +x
execution_dir=$(pwd)
who=$(whoami)



### installing lammps
# https://docs.lammps.org/Install_linux.html
(set -x
    sudo add-apt-repository ppa:gladky-anton/lammps
    sudo add-apt-repository ppa:openkim/latest
    sudo apt-get update
    sudo apt-get install lammps-stable
    sudo apt-get update
    sudo apt-get install lammps-stable-doc
    sudo apt-get install lammps-stable-data
    sudo apt-get install openkim-models
set +x)



### installing ovito
# https://www.ovito.org/linux-downloads/
# https://www.ovito.org/manual/installation.html
mkdir "$OVITO_INSTALL_LOC"
cp "Files/$OVITO_VERSION.tar.xz" "$OVITO_INSTALL_LOC/$OVITO_VERSION.tar.xz"
cd "$OVITO_INSTALL_LOC"
# Unzip with `tar xJfv ovito-X.X.X.tar.xz`.
tar xJfv "$OVITO_VERSION.tar.xz"

# Update environment variables.
echo "Updating environment variables for $who..."
echo "export PATH=\"$OVITO_INSTALL_LOC/$OVITO_VERSION/bin:$PATH\"" >> "~/.bashrc"
(set -x; source ~/.bashrc)



# copy/paste necessary files for example
cd "./Files"
cp "./Cu.meam" "../Calculations/RescaleUpload/Cu.meam"
cp "./Dislocation.f90" "../Calculations/RescaleUpload/Dislocation.f90"
cp "./DisVelocity.in" "../Calculations/RescaleUpload/DisVelocity.in"
cp "./library.meam" "../Calculations/RescaleUpload/library.meam"