#!/bin/bash

#PBS -A dp033
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=1

module load gsl/intel/2.4
module load intel/compilers/18.0.3
module load intel/mpi/18.0.3

cd $PBS_O_WORKDIR

gunzip -v phixsdata.txt.gz transitiondata.txt.gz ratecoeff.dat.gz packets*.out.gz

./exspec

gzip -v emission*.out absorption.out

gzip -v packets*.out
mkdir packets
mv packets*.out.gz packets/
