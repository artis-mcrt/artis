#!/bin/bash

#PBS -P fm5
#PBS -q express
#PBS -l walltime=12:00:00
#PBS -l mem=8GB
#PBS -l ncpus=1
#PBS -l wd
#PBS -m abe

# ncpus must be a factor of the cores per node (16)
# mem is total memory (all cores combined)
# on raijin, max 128GB/node = 8GB per core and MPI task (2% of nodes)
# max 64GB/node = 4GB per core (31% of nodes)
# max 32GB/node = 2GB per core (66% of nodes)
# normal queue walltime limits
# 48 hrs for 1-255 cores
# 24 hrs for 256-511 cores
# 10 hrs for 512-1023 cores
# 5 hours for 1024-56960 cores

ulimit -l 2097152

module load gsl/2.5
module load intel-mpi/2018.3.222
module load intel-cc/2018.3.222

gunzip -v phixsdata.txt.gz transitiondata.txt.gz ratecoeff.dat.gz
gunzip packets*.out.gz

./exspec

if [ -f spec.out ]; then
  gzip -v phixsdata.txt.gz transitiondata.txt.gz ratecoeff.dat.gz
  gzip -v packets*.out.gz
  mkdir packets
  mv packets*.out.gz packets/
fi