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
# Unzip with `tar xJfv ovito-X.X.X.tar.xz`.
echo "+ cd $OVITO_INSTALL_LOC"
echo "+ tar xJfv $OVITO_VERSION.tar.xz"
(set -x;
    cd "$OVITO_INSTALL_LOC"
    tar xJfv "$OVITO_VERSION.tar.xz"
)&> "$execution_dir/ovito_untar.log"

# Update environment variables.
echo "Updating environment variables for $who..."
echo "export PATH=\"$OVITO_INSTALL_LOC/$OVITO_VERSION/bin:\$PATH\"" >> ~/.bashrc
(set -x;
    source ~/.bashrc
    sleep 5s
    ovito --version
)

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

cp "datain" "../Calculations/2-MDDP/"
cp "data" "../Calculations/2-MDDP/"


wget http://archive.ubuntu.com/ubuntu/pool/universe/g/gcc-6/libgfortran3_6.4.0-17ubuntu1_amd64.deb
sudo dpkg -i libgfortran3_6.4.0-17ubuntu1_amd64.deb

(set -x;
    sudo apt install p7zip-full
    7z x "MDDP.7z"
)
cd "MDDP"
chmod +x "BCCdata"
chmod +x "FCCdata"
chmod +x "MDDP08-2008"
cp "BCCdata" "../../Calculations/2-MDDP/"
cp "FCCdata" "../../Calculations/2-MDDP/"
cp "MDDP08-2008" "../../Calculations/2-MDDP/"