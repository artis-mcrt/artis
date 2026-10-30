# Runs on a cluster

Each simulation has its own model folder. Keep the source code and the git metadata in a subfolder `artis` of the model
folder, not only the executables. The code then finds its data files in `artis/data`.

Make the model folder. Then clone the release branch into it:

```sh
mkdir mymodel
cd mymodel
git clone --branch release https://github.com/artis-mcrt/artis.git
cd artis
```

An HPC cluster gives the compiler and the MPI library through a module or a spack command. For the systems that we use,
the top of the Slurm script in `scripts/artis-*.sh` names the module specifications that work. For Open MPI, select the
compiler with `export OMPI_CXX=g++`.

Next, symlink an options preset to `artisoptions.h`:

```sh
ln -s artisoptions_classic.h artisoptions.h
```

You will likely want to change the number of packets of all ranks together (`NUM_PACKETS`). Use a text editor, e.g.
`vim artisoptions.h`. [Option presets](presets.md) says what each preset is for, and
[Compile-time options](options.md) describes every option.

Most presets have `NUM_PACKETS` values for production runs with approximately 1000 ranks. The `nltephotospheric`
preset has a much smaller value. Increase `NUM_PACKETS` for a production run with that preset. Each rank keeps its
share of the packets in memory. Decrease `NUM_PACKETS` for a run with fewer ranks.

Next, compile and go up a level to the model folder:

```sh
make
cd ..
```

Then make sure that the model folder holds a full set of input files: `model.txt`, `abundances.txt`, and an atomic
database (see [Atomic and model data](data.md)). Configure the timesteps in `input.txt` (see
[Input files](inputfiles.md)).

## Submit the job

Edit the job script of your cluster, e.g. the rank count, the account, and your email address. Then queue it with
`sbatch`. For example, on the Cosma8 cluster:

```sh
vim artis/scripts/artis-cosma8.sh
sbatch artis/scripts/artis-cosma8.sh
```

On Leonardo and Virgo, run `artis/scripts/artis-<cluster>-submit.sh` in place of `sbatch`. The Virgo submit script
holds the rank count. These two scripts take the email address for the Slurm messages from the `EMAIL` environment
variable (see `scripts/aliases.sh`).

On each of these clusters, `artis/scripts/submit.sh` finds the cluster from the host name and submits the applicable
job script. Run it in the model folder in place of `sbatch`.

## What the job scripts do

The job scripts resubmit themselves until the simulation is finished. They then queue a post-processing job
(`scripts/exspec-zip-*.sh`) that merges and compresses the output files. `sn3d` writes the spectra and light curves
itself, so this job runs `exspec` only if `emission.out` is absent. The post-processing job also converts the packets
and the estimator files to parquet format with [artistools](https://github.com/artis-mcrt/artistools), which it
installs automatically with uv.

[Bundled scripts](scripts.md) lists every script, and [Post-processing with exspec](exspec.md) says what `exspec`
does.
