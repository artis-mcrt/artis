#!/bin/bash
# each AMD node has 192GB RAM and 2x64 core CPUs
#--ntasks-per-node=48 will give 192/48=4GB per MPI rank
#--ntasks-per-core=1 gives 1.5GB per rank (probably fine for 1D models)
sbatch -J $(basename $(exec pwd)) --ntasks=960 --ntasks-per-node=48 --constraint=amd --partition=long --time=24:00:00 --mail-type=ALL --mail-user=${USER}@gsi.de --no-requeue -- artis/scripts/artis-virgo-slurmjob.sh
