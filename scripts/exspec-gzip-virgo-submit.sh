#!/bin/bash
sbatch -J exspec_$(basename $(exec pwd)) --ntasks=1 --cpus-per-task 48 --mem-per-cpu=4000 --partition=long --time=48:00:00 --constraint=amd,epyc,9654 --mail-type=ALL --mail-user=${USER}@gsi.de -- artis/scripts/exspec-gzip-virgo-slurmjob.sh
