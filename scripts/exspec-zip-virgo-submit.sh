#!/bin/bash -x
sbatch -J exspec_$(basename $(pwd)) --ntasks=1 --cpus-per-task 48 --mem-per-cpu=4000MB --partition=long --time=48:00:00 --constraint="[9654,9555]" --mail-type=ALL --mail-user=${USER}@gsi.de -- artis/scripts/exspec-zip-virgo-slurmjob.sh
