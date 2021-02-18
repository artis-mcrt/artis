sbatch -J artis --ntasks=960 --partition=main --time=08:00:00 --mail-type=ALL --mail-user=${USER}@gsi.de -- artis/scripts/artis-virgo-slurmjob.sh
