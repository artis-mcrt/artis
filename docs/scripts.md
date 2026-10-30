# Bundled scripts

The `scripts` folder holds the job scripts of the clusters that we use and some helper tools.

## The job scripts

- `artis-<cluster>.sh`: the Slurm job scripts of `sn3d`. On Leonardo and Virgo, `artis-<cluster>-submit.sh` submits
  `artis-<cluster>-slurmjob.sh`.
- `submit.sh`: it finds the cluster from the host name and submits the `sn3d` job script of that cluster.
- `exspec-zip-<cluster>.sh`: the post-processing job. It runs `run-exspec-if-needed.sh` and then `exspec-after.sh`. On
  Leonardo and Virgo, `exspec-zip-<cluster>-submit.sh` submits `exspec-zip-<cluster>-slurmjob.sh`.
- `run-exspec-if-needed.sh`: it runs `exspec` if `emission.out` is absent.
- `exspec-after.sh`: after the last timestep of the model, the script does these steps:
    - It removes the restart files.
    - It merges the direction-resolved files with `mergeangleres.py`.
    - It compresses the large text files with zstd.
    - It archives the logs with `tar_rm_logs.sh`.
    - It makes the parquet caches of artistools with `artistools cache make`, e.g. of the estimators and the packets.
- `corehours-before.sh` and `corehours-after.sh`: the job scripts source these files to log the CPU allocation, the
  wall time, and the core hours.

## The helper tools

- `slurmjobhoursleft.py`: print the hours that remain in a Slurm job. The job scripts give this value to `sn3d -w`.
- `mergeangleres.py`: merge the files of the direction bins into `light_curve_res.out`, `spec_res.out`, and
  `specpol_res.out`.
- `combine_allrank_out_files.py`: make the combined files of all ranks (see [Output files](outputfiles.md)).
- `tar_rm_logs.sh`: the script archives the log files of the ranks into `logfiles.tar.zst`. It then removes them,
  except `output_0-0.txt`.
- `rsyncnooutput.sh`: copy a model folder with rsync, without the output files.
- `aliases.sh`: shell aliases for Slurm and the `MAKEFLAGS` setting. Replace the placeholder email address in the file
  with your address. Then source the file from your `.bashrc`.
- `clean.sh`: remove all output files, keep the input files, and reset the simulation to the beginning. The script also
  removes the `job_from_ts*` folders of `sn3d` and the `*.slurm` folders of older versions.
- `sumcorehoursslurm.py`: calculate the summed core hours of all jobs from the Slurm job output files.
- `make_options_docs.py`: write the [Compile-time options](options.md) page from `artisoptions_doc.md` and the presets.
