#!/bin/bash
sbatch -J exspec_$(basename $(exec pwd)) --ntasks=1 --cpus-per-task 16 --mem-per-cpu=4096MB --partition=long --time=48:00:00 --mail-type=ALL --mail-user=${USER}@gsi.de -- artis/scripts/exspec-gzip-virgo-slurmjob.sh
