#!/bin/bash
sbatch -J ex_$(basename $(exec pwd)) --ntasks=1 --cpus-per-task 8 --mem-per-cpu=2048MB --partition=long --time=48:00:00 --mail-type=ALL --mail-user=${USER}@gsi.de -- artis/scripts/exspec-gzip-virgo-slurmjob.sh
