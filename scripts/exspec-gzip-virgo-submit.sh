sbatch -J artis --ntasks=1 --partition=long --time=02:00:00 --mail-type=ALL --mail-user=${USER}@gsi.de -- artis/scripts/exspec-gzip-virgo-slurmjob.sh
