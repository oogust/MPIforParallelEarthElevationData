#!/bin/bash

export BUILD_OPENBLAS=1
export BUILD_SCALAPACK=1
export BLAS_SIZE=8
export SCALAPACK_SIZE=8

for hostname in $(cat $OAR_NODEFILE | uniq)
do 
    oarsh $hostname "sudo-g5k apt-get -y install libscalapack-openmpi-dev" &
done
