# Runs on a cluster

Each simulation has its own {term}`model folder`. Keep the source code and the git metadata in a subfolder `artis` of the model
folder, not only the executables. The code then finds its data files in `artis/data`.

## Which version to use

The repository has two branches and a tag for each release:

- `main` holds the newest code. The version [latest](https://artis-mcrt.readthedocs.io/en/latest/) of this site
  describes it.
- `release` is the branch for production runs. It holds the code of the newest release tag.
- Each release has a tag with its date, e.g. `v2026.11.9`. The version
  [stable](https://artis-mcrt.readthedocs.io/en/stable/) of this site describes the newest tag.

Use the branch `release` for a production run, and read the version `stable` of this site. For the results of a
paper, use one tag, so that a later change of `release` does not change the code of the model folder. Run these
commands in the new model folder, with the tag of your choice. Do this before you change `input.txt` and before you
run `make`.

```sh
git -C artis checkout v2026.11.9
cp artis/scripts/input-template.txt input.txt
```

The second command copies the template of the tag, because the lines of `input.txt` can change between versions. The
tags before `v2026.11.9` have no template.

The log of rank 0 of each job starts with the git version of the build. Record that version with the results:

```sh
grep -E -m 2 'git (version|status):' output_0-0.txt
```

A version that ends with `-dirty`, or a `git status` line that names files, says that the clone had changes outside
the commits of git.

## Make the model folder

`scripts/new-model-folder.sh` of a clone makes the model folder. It clones the branch `release` into the subfolder
`artis` of the new folder, and it copies the template of `input.txt`. Make a first clone if you have none:

```sh
git clone https://github.com/artis-mcrt/artis.git ~/artis
~/artis/scripts/new-model-folder.sh mymodel
cd mymodel/artis
```

An HPC cluster gives the compiler and the MPI library through a module or a spack command. For the systems that we use,
the top of the Slurm script in `scripts/artis-*.sh` names the module specifications that work. For Open MPI, select the
compiler with `export OMPI_CXX=g++`.

Next, symlink an options {term}`preset` to `artisoptions.h`:

```sh
ln -s artisoptions_classic.h artisoptions.h
```

[Option presets](presets.md) says what each {term}`preset` is for, and [Compile-time options](options.md) describes every
option. Set the number of packets of all {term}`ranks <rank>` together (`NUM_PACKETS`) for your rank count, e.g. with
`vim artisoptions.h`. The value of a preset is not always correct for a production run (see [Option presets](presets.md#which-preset)).

Next, compile and go up a level to the model folder:

```sh
make
cd ..
```

Then make sure that the model folder holds a full set of input files: `model.txt`, `abundances.txt`, an atomic
database (see [Atomic and model data](data.md)), and `input.txt`. The `input.txt` of the template holds the values of
a test model. Set its lines 1 to 3, 18, and 19 for your model: the {term}`timesteps <timestep>`, the times, and the {term}`LTE` and grey
timesteps (see [Input files](inputfiles.md)).

## Submit the job

Export these variables in your `.bashrc`. Set only the variables that your cluster needs, because they apply to each
`sbatch` command of your shell.

```sh
export EMAIL=your_email_address
export SBATCH_ACCOUNT=your_project_account
export SBATCH_PARTITION=your_partition
```

`EMAIL` gives the address for the messages of Slurm. `sbatch` reads `SBATCH_ACCOUNT` and `SBATCH_PARTITION` at each
submission. A job gets the environment of its submission by default (the option `--export=ALL` of `sbatch`), so the
job scripts submit the next job with the same values. These variables replace the `#SBATCH` lines of the job script,
and an option on the command line of `sbatch` replaces the variables (see the section "INPUT ENVIRONMENT VARIABLES" of
the [sbatch manual](https://slurm.schedmd.com/sbatch.html)).

Then run `artis/scripts/submit.sh` in the model folder:

```sh
artis/scripts/submit.sh
```

The script finds the cluster from the host name. It knows Cosma8, Kelvin2, Leonardo, and Virgo, and it gives `EMAIL`
to `sbatch`. Each job script gives `EMAIL` to `sbatch` again when it submits the next job. On Leonardo and Virgo, the
script calls the submit script of the cluster, `artis/scripts/artis-<cluster>-submit.sh`, because `sbatch` must get
some options on the command line there:

- The Virgo submit script gives the rank count, the memory, the partition, and the wall time. `SBATCH_PARTITION`
  thus does not change the partition there.
- The Leonardo submit script gives only the job name and the email address. The other values come from the `#SBATCH`
  lines of `scripts/artis-leonardo-slurmjob.sh`.

The job script holds the other values, e.g. the rank count and the wall time. Each model folder has its own clone, so
an edit of the job script changes only the values of that model folder. The `git status` line of the log then names
the changed script.

Keep the default name `slurm-<job ID>.out` for the output file of the job. As an alternative, put that file outside
the model folder. A new simulation removes all other files of the model folder that have the ending `.out`.

## What the job scripts do

The job scripts resubmit themselves until the simulation is finished. They then queue a post-processing job
(`scripts/exspec-zip-*.sh`) that merges and compresses the output files. A job that stops at a `timestep_finish` before
the last timestep of the model queues no job. `sn3d` writes the spectra and light curves
itself, so this job runs `exspec` only if `emission.out` is absent. The post-processing job also converts the packets
and the {term}`estimator` files to parquet format with [artistools](https://github.com/artis-mcrt/artistools), which it
installs automatically with uv.

## Port the job scripts to a new cluster

The scripts of `scripts/` name the modules and the partitions of four clusters. For a different cluster, copy
`scripts/artis-cosma8.sh` to a new file, e.g. `scripts/artis-mycluster.sh`. Change these lines of the copy:

- the `#SBATCH` lines of the partition, the account, the rank count, and the wall time, unless the variables of
  [Submit the job](#submit-the-job) set them;
- the `module` lines;
- the two `sbatch` commands at the end. They name `artis/scripts/artis-cosma8.sh` and
  `artis/scripts/exspec-zip-cosma8.sh`, so a copy that keeps them submits the jobs of Cosma8.

A job script must do these steps:

1. Set the Slurm values, e.g. the rank count, the wall time, the partition, and the account.
2. Load a compiler and an MPI library, or set `OMPI_CXX`.
3. Go to the model folder. Slurm gives it in `SLURM_SUBMIT_DIR`.
4. Run `make sn3d exspec` in `artis`, so that `exspec` has the same source version as `sn3d`.
5. Get the hours that remain with `scripts/slurmjobhoursleft.py`, and give them to `sn3d -w`.
6. Stop if the status of `mpirun` is not zero, because that is a crash. `sn3d` gives zero when it writes
   `RESTART_NEEDED`.
7. Search `output_0-0.txt` for `RESTART_NEEDED`. If the text is there, submit the same job script again.
8. If that text is absent, search `output_0-0.txt` for `No need for restart`. `sn3d` writes it after the last timestep
   of the model. If the text is there, submit the post-processing job.

Copy `scripts/exspec-zip-cosma8.sh` for the post-processing job, and change its modules. Add the host name of the
cluster to `scripts/submit.sh` if you use that script.

[Bundled scripts](scripts.md) lists every script, and [Post-processing with exspec](exspec.md) says what `exspec`
does.

Next: [Analyse the results](analysis.md).
