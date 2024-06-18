#!/bin/bash -x

sudo apt-get -y remove bc pv gfortran build-essential \
    quantum-espresso python3-pip python2 curl expect p7zip-full

rm -rf ~/.juliaup
rm -rf /usr/bin/.pw.x
rm -rf ~/QuantumEspresso
rm -rf ~/Ovito

sudo apt-get autoremove
