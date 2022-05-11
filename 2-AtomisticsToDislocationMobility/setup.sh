#!/bin/bash -x



# define path into where Quantum Espresso will install
OVITO_INSTALL_LOC="~/Ovito"
# encode name and version of of tarball
OVITO_VERSION="ovito-basic-3.7.2-x86_64"

execution_dir=$(pwd)
NUM_PROC=2



### ================ ###
### Installing Ovito
### ================ ###
# https://www.ovito.org/linux-downloads/
# https://www.ovito.org/manual/installation.html

# Copy/Paste the `Files/qe-X.X.X.tar.gz`
# archive into a working directory.
mkdir "$OVITO_INSTALL_LOC"
cp "Files/$OVITO_VERSION.tar.gz" "$OVITO_INSTALL_LOC/$OVITO_VERSION.tar.gz"
cd "$OVITO_INSTALL_LOC"
# Unzip with `tar -xzvf ovito-X.X.X.tar.gz`.
tar xJfv "$OVITO_VERSION.tar.gz"

# To set path to `pw.x` as environment variable
# change PATH as needed to QE `/bin/` folder.
# cd ~ && gedit .bashrc && export PATH="$OVITO_INSTALL_LOC/q-e-$OVITO_VERSION/bin:$PATH"
echo "export PATH=\"$OVITO_INSTALL_LOC/$OVITO_VERSION/bin:$PATH\"" >> "~/.bashrc"

# Update environment variables.
source .bashrc



# copy/paste necessary files for example
cd "./Files"
cp "./Cu.meam" "../Calculations/RescaleUpload/Cu.meam"
cp "./Dislocation.f90" "../Calculations/RescaleUpload/Dislocation.f90"
cp "./DisVelocity.in" "../Calculations/RescaleUpload/DisVelocity.in"
cp "./library.meam" "../Calculations/RescaleUpload/library.meam"