#!/bin/bash -x



# define path into where Ovito will install
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



# ### installing lammps
# # https://docs.lammps.org/Install_linux.html
# (set -x;
#     sudo add-apt-repository ppa:gladky-anton/lammps
#     sudo add-apt-repository ppa:openkim/latest
#     sudo apt-get update
#     sudo apt-get install lammps-stable
#     sudo apt-get update
#     sudo apt-get install lammps-stable-doc
#     sudo apt-get install lammps-stable-data
#     sudo apt-get install openkim-models
# )



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
echo "export PATH=\"$OVITO_INSTALL_LOC/$OVITO_VERSION/bin:\$PATH\"" >> ~/.bashrc
(set -x; source ~/.bashrc)

# copy/paste necessary files for example
cd "$execution_dir/Files"
mkdir "test"
mkdir "test/RescaleUpload"
cp "../Calculations/0-Scripts/Cu.meam" "./test/RescaleUpload/Cu.meam"
cp "../Calculations/0-Scripts/Dislocation.f90" "./test/RescaleUpload/Dislocation.f90"
cp "../Calculations/0-Scripts/DisVelocity.in" "./test/RescaleUpload/DisVelocity.in"
cp "../Calculations/0-Scripts/library.meam" "./test/RescaleUpload/library.meam"
cp "../Calculations/0-Scripts/rescale_commands.sh" "./test/RescaleUpload/rescale_commands.sh"
cp "../Calculations/0-Scripts/atoms.sh" "./test/RescaleUpload/atoms.sh"

cp "../Calculations/0-Scripts/datain" "../Calculations/2-MDDP/"