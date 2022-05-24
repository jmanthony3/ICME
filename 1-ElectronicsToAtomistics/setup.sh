#!/bin/bash -x



# define install path for Quantum Espresso (QE)
QUANTUM_ESPRESSO_INSTALL_LOC=~/QuantumEspresso
# encode name and version of tarball: qe-X.X.X
QUANTUM_ESPRESSO_VERSION="qe-6.0.0"
# number of processors to use in test case
NUM_PROC=$(nproc) # grabs all cores available by default




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - END OF HUMAN EDITABLE SECTION - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #




set +x
execution_dir=$(pwd)
who=$(whoami)

### installing quantum espresso
# add `cmake`, `gcc`, `gfortran`, and `make` capabilities.
echo "Updating distro and including cmake, gcc, gfortran, and make..."
(set -x;
    sudo apt update
    sudo apt install build-essential
)

# copy/paste the `Files/qe-X.X.X.tar.gz` archive into a working directory
mkdir "$QUANTUM_ESPRESSO_INSTALL_LOC"
cp "$execution_dir/Files/$QUANTUM_ESPRESSO_VERSION"*".tar.gz" "$QUANTUM_ESPRESSO_INSTALL_LOC/"
cd "$QUANTUM_ESPRESSO_INSTALL_LOC"
# unzip with `tar -xzvf qe-X.X.X.tar.gz`
(set -x; tar -xzvf "$QUANTUM_ESPRESSO_VERSION.tar.gz")

# `cd` into that extracted folder and execute `./configure && make all`
(set -x;
    cd "q-e-$QUANTUM_ESPRESSO_VERSION"
    ./configure
    make all %> "$execution_dir/quantum_espresso_make.log"
)

# set `pw.x` as environment variable change PATH as needed to QE `/bin/` folder
echo "export PATH=\"$QUANTUM_ESPRESSO_INSTALL_LOC/q-e-$QUANTUM_ESPRESSO_VERSION/bin:\$PATH\"" >> ~/.bashrc

# Update environment variables.
echo "Updating environment variables for $who..."
(set -x;
    source ~/.bashrc
    # pw.x --version
)



### test execution of `EvA_EvV_plot.py`
echo "Installing by 'sudo apt' only for mpi dependencies..."
(set -x; sudo apt install quantum-espresso)
cd "$execution_dir/Files"
mkdir "test"
echo "Executing QE according to Cu.in..."
(set -x;
    # executes the input parameters with QE
    mpirun -np $NUM_PROC pw.x -in "Cu.in" > "./test/Cu.out"
    # compiles `evfit.f` outputs `ev_curve`
    gfortran -O2 "evfit.f" -o "evfit"
)
cp "Cu.in" "fcc.ev.in" # create appropriate input file to `ev_curve`
chmod +x ./ev_curve # makes file executable
(set -x; ./ev_curve fcc 3.628) # reference structure, lattice parameter
echo "Ensuring pip3 capabilities for matplotlib and numpy..."
(set -x;
    sudo apt install python3-pip
    python3 -m pip install matplotlib numpy
    python3 "EvA_EvV_plot.py" # generate plots
)
mv "evfit" "./test/"
mv "fcc.ev.in" "./test/"
mv "EvsA" "./test/"
mv "EvsV" "./test/"
mv "SUMMARY" "./test/"
mv "pw_ev.out" "./test/"
mv "Name_of_EvA.pdf" "./test/"
mv "Name_of_EvV.pdf" "./test/"
mv "Name_of_Combined.pdf" "./test/"
mv "evfit.4" "./test/"
mv "pw_ev.out" "./test/"
rm -r "temp/"



### test execution of `gsfe_curve.py`
echo "Installing Python 2..."
(set -x;
    sudo add-apt-repository universe
    sudo apt update
    sudo apt install python2
    sudo apt install curl
    curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py
    sudo python2 get-pip.py
    pip2 --version
    pip2 install numpy
)
# reference structure, lattice parameter, and block motion
(set -x; python2 "gsfe_curve.py" fcc 3.615 partial &)
sleep 5s
pid=$(pgrep pw.x)
echo "Killing the 'gsfe_curve.py' process ('PID=$pid') because this will take too long..."
kill $pid
mv "gsfe.in" "./test/gsfe.in"
mv "gsfe.out" "./test/gsfe.out"
mv "GSFE_SUMMARY" "./test/GSFE_SUMMARY"
rm -r "temp/"



### populate "../Calculations/0-Scripts/" folder
cp "ev_curve" "../Calculations/0-Scripts/"
cp "EvA_EvV_plot.py" "../Calculations/0-Scripts/"
cp "evfit.f" "../Calculations/0-Scripts/"
cp "gsfe_curve.py" "../Calculations/0-Scripts/"
cp "OutputFileCreator.py" "../Calculations/0-Scripts/"
cp "OutputSummarizer.py" "../Calculations/0-Scripts/"
cp "rescale_commands.sh" "../Calculations/0-Scripts/"





# that's all folks