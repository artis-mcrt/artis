#!/bin/bash
sbatch -J $(basename $(exec pwd)) --ntasks=960 --ntasks-per-node=32 --partition=long --time=24:00:00 --mail-type=ALL --mail-user=${USER}@gsi.de -- artis/scripts/artis-virgo-slurmjob.sh
