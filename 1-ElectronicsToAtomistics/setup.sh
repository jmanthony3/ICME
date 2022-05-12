#!/bin/bash -x



# define install path for Quantum Espresso (QE)
QUANTUM_ESPRESSO_INSTALL_LOC=~/QuantumEspresso
# encode name and version of tarball: qe-X.X.X
QUANTUM_ESPRESSO_VERSION="qe-6.0.0"
# number of processors to use in test case
NUM_PROC=2




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - END OF HUMAN EDITABLE SECTION - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #




execution_dir=$(pwd)
### installing quantum espresso
# add `cmake`, `gcc`, `gfortran`, and `make` capabilities.
sudo apt update
sudo apt install build-essential

# copy/paste the `Files/qe-X.X.X.tar.gz` archive into a working directory
mkdir "$QUANTUM_ESPRESSO_INSTALL_LOC"
cp "$execution_dir/Files/$QUANTUM_ESPRESSO_VERSION*.tar.gz" "$QUANTUM_ESPRESSO_INSTALL_LOC/$QUANTUM_ESPRESSO_VERSION.tar.gz"
cd "$QUANTUM_ESPRESSO_INSTALL_LOC"
# unzip with `tar -xzvf qe-X.X.X.tar.gz`
tar -xzvf "$QUANTUM_ESPRESSO_VERSION.tar.gz"

# `cd` into that extracted folder and execute `./configure && make all`
cd "q-e-$QUANTUM_ESPRESSO_VERSION"
./configure
make all

# set `pw.x` as environment variable change PATH as needed to QE `/bin/` folder
echo "export PATH=\"$QUANTUM_ESPRESSO_INSTALL_LOC/q-e-$QUANTUM_ESPRESSO_VERSION/bin:$PATH\"" >> ~/.bashrc

# Update environment variables.
source ~/.bashrc



### test execution of `EvA_EvV_plot.py`
sudo apt install quantum-espresso
cd "$execution_dir/Files"
mkdir "test"
mpirun -np $NUM_PROC pw.x -in "Cu.in" > "./test/Cu.out" # executes the input parameters with QE
gfortran -O2 "evfit.f" -o "evfit" # compiles `evfit.f` outputs `ev_curve`
cp "Cu.in" "fcc.ev.in" # create appropriate input file to `ev_curve`
chmod +x ./ev_curve # makes file executable
./ev_curve fcc 3.628 # reference structure, lattice parameter
python3 "EvA_EvV_plot.py" # generate plots
mv "fcc.ev.in" "./test/fcc.ev.in"
mv "EvsA" "./test/EvsA"
mv "EvsV" "./test/EvsV"
mv "SUMMARY" "./test/SUMMARY"
mv "pw_ev.out" "./test/pw_ev.out"
mv "Name_of_EvA.pdf" "./test/Name_of_EvA.pdf"
mv "Name_of_EvV.pdf" "./test/Name_of_EvV.pdf"
mv "Name_of_Combined.pdf" "./test/Name_of_Combined.pdf"
rm -rf "temp/"



### test execution of `gsfe_curve.py`
sudo add-apt-repository universe
sudo apt update
sudo apt install python2
sudo apt install curl
curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py
sudo python2 get-pip.py
pip2 --version
pip2 install numpy
# reference structure, lattice parameter, and block motion
python2 "gsfe_curve.py" fcc 3.615 partial &
pid=$!
echo "Killing the 'gsfe_curve.py' process ('$pid') because this will take too long."
kill $pid
mv "gsfe.in" "./test/gsfe.in"
mv "gsfe.out" "./test/gsfe.out"