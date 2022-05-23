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
cp "Cu.meam" "./test/RescaleUpload/"
cp "Dislocation.f90" "./test/RescaleUpload/"
cp "DisVelocity.in" "./test/RescaleUpload/"
cp "library.meam" "./test/RescaleUpload/"
cp "rescale_commands.sh" "./test/RescaleUpload/"
cp "atoms.sh" "./test/RescaleUpload/"

cp "datain" "../Calculations/0-Scripts/"
cp "datain" "../Calculations/2-MDDP/"
cp "MDDP_BCC_HW2.7z" "../Calculations/2-MDDP/"
