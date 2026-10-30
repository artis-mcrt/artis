# Output files

A simulation uses two types of folder:

- The model folder is the top-level folder. It holds the input files, the restart files, and the files of the whole
  simulation.
- A job folder is a subfolder of the model folder. `sn3d` writes one job folder in each Slurm job. `sn3d` names the
  folder from the start timestep of the job, e.g. `job_from_ts0000` for a new simulation and `job_from_ts0008` for a
  job that resumes at timestep 8. Rank 0 writes the name of the folder to the standard output.

A build with libzstd writes each `.out` file zstd compressed, e.g. `spec.out.zst`, so that no separate compression step
is necessary. The logs, `input.txt`, and the restart files stay plain. Every file gets zstd level 9.
[artistools](https://github.com/artis-mcrt/artistools) reads both forms.

## The files of each job

Each job writes these files into its job folder:

- `output_n-0.txt`: a log file for each MPI rank n.
- the estimators of each cell at each timestep: the temperatures, the ionisation, and the heating and cooling rates.
  The option `WRITE_COMBINED_ALLRANK_OUT_FILES` selects one of two forms, and every preset sets the option:
    - with the option, `sn3d` writes the text of all ranks into one combined file, `estimators_allranks.out`;
    - without the option, each rank writes its own file, `estimators_nnnn.out`.

    The option gives the same choice for the nlte and the radfield files, e.g. `nlte_allranks.out`. With the constant
    `LOG_MACROATOM` in `macroatom.cc`, it also gives `macroatom_allranks.out`.
    `scripts/combine_allrank_out_files.py` makes the combined files from the files of the ranks, and it keeps the
    files of the ranks unless you give `--rm`. Run the script with uv, e.g. `uv run artis/scripts/combine_allrank_out_files.py`.
  The script needs Python 3.14 or a later version, and uv gets that version from the metadata of the script.

## The files of the whole simulation

A run writes these files into the model folder:

- `packets/packets00_nnnn.out`: the Monte Carlo packets from each rank, which `exspec` can turn into spectra and light
  curves again.
- `light_curve.out`, `spec.out`, and the other spectrum files that [Post-processing with exspec](exspec.md) lists:
  `sn3d` writes the light curves and spectra at each timestep, and the emission, absorption, and direction-resolved
  files at the last requested timestep.
- `deposition.out`: the radioactive energy deposition rate as a function of time. The columns of the deposited
  luminosities come from the Monte Carlo estimators, also when `PARTICLE_THERMALISATION_SCHEME` is
  `INSTANTFULLDEPOSITION` and the heating of the cells uses the analytic emission rates. For a beta-plus decay, the
  `Qdot_ana_erg/s/g` column uses the energy Q_EC, which includes the annihilation of the positron.
- `linelist.out`: the list of the atomic lines. The first row is a header that starts with "#" and names the columns:
  the vacuum wavelength in Angstroms, the atomic number, the ion stage, the upper level index, and the lower level
  index. The level indices start at 0 in each ion. Each later row is one line, in the order of the line index that the
  emission and absorption codes of the packets use. Only the first job of a simulation writes this file.
- `gridsave_ts*.tmp` and `packets_*_ts*.tmp`: restart files that hold the state at the start of a timestep, after its
  grid update. A later job continues from that timestep and propagates its packets.

`sn3d` also keeps an `output_0-0.txt` symlink in the model folder that points to the rank-0 log of the current job.

:::{warning}
A restart file must only be consistent with the binary that wrote it. A resumed run needs the same
`artisoptions.h` and the same source version as the run that wrote the files.
:::

## A new simulation removes the files of the previous one

A new simulation first removes the output files, the restart files, and the `job_from_ts*` and `*.slurm` folders of the
previous simulation. It removes the same files as `scripts/clean.sh`, except `slurm-*.out`, `machine.file.*`, and
`core.*`, which can belong to the current job.

## Plotting and analysis

[artistools](https://github.com/artis-mcrt/artistools) reads these files. Some commands for a first look:

```sh
artistools describe <model folder>
artistools check <model folder>
artistools plotlightcurves <model folder>
artistools plotspectra -t <days> <model folder>
artistools plotestimators <model folder>
```

`artistools filterpackets` removes each packet except the escaped r-packets from the packet files, which makes the
files smaller. `artistools corehours` sums the core hours of the jobs from the `output_0-0.txt` log of each job.
