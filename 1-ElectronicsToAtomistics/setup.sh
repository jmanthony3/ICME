#!/bin/bash -x



#################### INSTALLATION VARIABLES ###################
# define install path for Quantum Espresso (QE)
QUANTUM_ESPRESSO_INSTALL_LOC=~/QuantumEspresso
# encode name and version of tarball: qe-X.X.X
QUANTUM_ESPRESSO_VERSION="qe-6.0.0"
# number of processors to use in test case
NUM_PROC=$(nproc) # grabs all cores available by default





# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - END OF HUMAN EDITABLE SECTION - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #





set +x # turn script tracing off
execution_dir=$(pwd) # where script executes from
who=$(whoami) # current user
mkdir "$execution_dir/logs"



####################### INSTALL SOFTWARE ######################
### installing quantum espresso
# add `cmake`, `gcc`, `gfortran`, and `make` capabilities.
echo "Updating distro and including cmake, gcc, gfortran, and make..."
(set -x;
    sudo apt-get update
    sudo apt-get install build-essential
)
sudo apt-get install pv

# copy tarball into installation directory
mkdir "$QUANTUM_ESPRESSO_INSTALL_LOC"
cp "$execution_dir/Files/$QUANTUM_ESPRESSO_VERSION"*".tar.gz" \
    "$QUANTUM_ESPRESSO_INSTALL_LOC/"

# unzip with `tar -xzvf qe-X.X.X.tar.gz`
(set -x;
    cd "$QUANTUM_ESPRESSO_INSTALL_LOC"
    tar -xzvf "$QUANTUM_ESPRESSO_VERSION.tar.gz"
# show progress of untar and write log
) 2> "$execution_dir/logs/quantum_espresso_untar.log" | pv -pterb --size 257353 > "$execution_dir/logs/quantum_espresso_untar.log"

# `cd` into that extracted folder and execute `./configure`
(set -x;
    cd "$QUANTUM_ESPRESSO_INSTALL_LOC/q-e-$QUANTUM_ESPRESSO_VERSION"
    ./configure
# show progress of configure and write log
) 2> "$execution_dir/logs/quantum_espresso_configure.log" | pv -pterb --size 5813 > "$execution_dir/logs/quantum_espresso_configure.log"

# `cd` into that extracted folder and execute `make all`
(set -x;
    cd "$QUANTUM_ESPRESSO_INSTALL_LOC/q-e-$QUANTUM_ESPRESSO_VERSION"
    make all
# show progress of make and write log
) 2> "$execution_dir/logs/quantum_espresso_make.log" | pv -pterb --size 752541 > "$execution_dir/logs/quantum_espresso_make.log"

# set `pw.x` as environment variable; change PATH as needed to QE `/bin/` folder
echo "Setting 'pw.x' as environment variable..."
echo "export PATH=\"$QUANTUM_ESPRESSO_INSTALL_LOC/q-e-$QUANTUM_ESPRESSO_VERSION/bin:\$PATH\"" >> ~/.bashrc

# update environment variables for user
echo "Updating environment variables for $who..."
(set -x; source ~/.bashrc)


### test execution of `EvA_EvV_plot.py`
# for mpi dependency
echo "Installing by 'sudo apt-get' only for mpi dependencies..."
(set -x; sudo apt-get install quantum-espresso)

# navigate back to appropriate directory
cd "$execution_dir/Files"
mkdir "test" # make test folder
echo "Executing QE according to Cu.in..."
(set -x;
    # execute QE with input parameters
    mpirun -np $NUM_PROC pw.x -in "Cu.in" > "./test/Cu.out" 2> /dev/null
    # compiles `evfit.f` outputs `ev_curve`
    gfortran -O2 "evfit.f" -o "evfit" 2> /dev/null
)

# create EvsA and EvsV curves
cp "Cu.in" "fcc.ev.in" # create appropriate input file to `ev_curve`
(set -x; ./ev_curve fcc 3.628 2> /dev/null) # reference structure, lattice parameter
echo "Ensuring pip3 capabilities for matplotlib and numpy..."
(set -x;
    sudo apt-get install python3-pip # install pip3
    python3 -m pip install matplotlib numpy # install modules
    python3 "EvA_EvV_plot.py" # generate plots
)
# move (in/out)put files to `./test/`
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
rm -r "temp/" # remove calculations temporary folder


### test execution of `gsfe_curve.py`
# install python2
echo "Installing Python 2..."
(set -x;
    sudo add-apt-repository universe
    sudo apt-get update
    sudo apt-get install python2
    sudo apt-get install curl
    curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py
    sudo python2 get-pip.py
    pip2 --version
    pip2 install numpy
)
# python2 gsfe_curve.py: reference structure, lattice parameter, block motion
(set -x; python2 "gsfe_curve.py" fcc 3.615 partial &)
sleep 10s # let previous process spin up
pid=$(pgrep pw.x) # get pid of `pw.x` process
echo "Killing the 'gsfe_curve.py' process ('PID=$pid') \
    because this will take too long..."
kill $pid
# move (in/out)put files to `./test/`
mv "gsfe.in" "./test/"
mv "gsfe.out" "./test/"
mv "GSFE_SUMMARY" "./test/"
rm -r "temp/" # remove calculations temporary folder


### populate `../Calculations/0-Scripts/`
cp "ev_curve" "../Calculations/0-Scripts/"
cp "EvA_EvV_plot.py" "../Calculations/0-Scripts/"
cp "evfit.f" "../Calculations/0-Scripts/"
cp "gsfe_curve.py" "../Calculations/0-Scripts/"
cp "OutputFileCreator.py" "../Calculations/0-Scripts/"
cp "OutputSummarizer.py" "../Calculations/0-Scripts/"
cp "rescale_commands.sh" "../Calculations/0-Scripts/"





# that's all folks