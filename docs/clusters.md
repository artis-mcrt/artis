# Runs on a cluster

Each simulation has its own {term}`model folder`. Keep the source code and the git metadata in a subfolder `artis` of the model
folder, not only the executables. The code then finds its data files in `artis/data`.

Make the model folder. Then clone the default branch `main` into it, and copy the template of `input.txt`:

```sh
mkdir mymodel
cd mymodel
git clone https://github.com/artis-mcrt/artis.git
cp artis/scripts/input-template.txt input.txt
cd artis
```

`scripts/new-model-folder.sh` of a clone that you have does the same steps, e.g. from a model folder next to the new
one:

```sh
artis/scripts/new-model-folder.sh ../mymodel
cd ../mymodel/artis
```

An HPC cluster gives the compiler and the MPI library through a module or a spack command. For the systems that we use,
the top of the Slurm script in `scripts/artis-*.sh` names the module specifications that work. For Open MPI, select the
compiler with `export OMPI_CXX=g++`.

Next, symlink an options {term}`preset` to `artisoptions.h`:

```sh
ln -s artisoptions_classic.h artisoptions.h
```

Set the number of packets of all {term}`ranks <rank>` together (`NUM_PACKETS`) with a text editor, e.g. `vim artisoptions.h`. [Option presets](presets.md) says what each preset is for, and
[Compile-time options](options.md) describes every option.

Most presets have `NUM_PACKETS` values for production runs with approximately 1000 ranks. The `nltephotospheric`
preset has a much smaller value. Increase `NUM_PACKETS` for a production run with that preset. Each rank keeps its
share of the packets in memory. Decrease `NUM_PACKETS` for a run with fewer ranks.

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

Edit the job script of your cluster, e.g. the rank count, the account, and your email address. Then queue it with
`sbatch`. For example, on the Cosma8 cluster:

```sh
vim artis/scripts/artis-cosma8.sh
sbatch artis/scripts/artis-cosma8.sh
```

On Leonardo and Virgo, run `artis/scripts/artis-<cluster>-submit.sh` in place of `sbatch`. The Virgo submit script
holds the rank count. These two scripts take the email address for the Slurm messages from the `EMAIL` environment
variable. Export it in your `.bashrc`.

On Cosma8, Kelvin2, Leonardo, and Virgo, `artis/scripts/submit.sh` finds the cluster from the host name and submits
the applicable job script. Run it in the model folder in place of `sbatch`.

## What the job scripts do

The job scripts resubmit themselves until the simulation is finished. They then queue a post-processing job
(`scripts/exspec-zip-*.sh`) that merges and compresses the output files. `sn3d` writes the spectra and light curves
itself, so this job runs `exspec` only if `emission.out` is absent. The post-processing job also converts the packets
and the {term}`estimator` files to parquet format with [artistools](https://github.com/artis-mcrt/artistools), which it
installs automatically with uv.

## Port the job scripts to a new cluster

The scripts of `scripts/` name the modules and the partitions of four clusters. For a different cluster, copy
`scripts/artis-cosma8.sh` and change the names. A job script must do these steps:

1. Set the Slurm values, e.g. the rank count, the wall time, the partition, and the account.
2. Load a compiler and an MPI library, or set `OMPI_CXX`.
3. Go to the model folder. Slurm gives it in `SLURM_SUBMIT_DIR`.
4. Run `make sn3d exspec` in `artis`, so that `exspec` has the same source version as `sn3d`.
5. Get the hours that remain with `scripts/slurmjobhoursleft.py`, and give them to `sn3d -w`.
6. Stop if the status of `mpirun` is not zero, because that is a crash. `sn3d` gives zero when it writes
   `RESTART_NEEDED`.
7. Search `output_0-0.txt` for `RESTART_NEEDED`. If the text is there, submit the same job script again.
8. If the text is absent and `packets/packets00_0000.out*` exists, the model is complete. Submit the post-processing
   job.

Copy `scripts/exspec-zip-cosma8.sh` for the post-processing job, and change its modules. Add the host name of the
cluster to `scripts/submit.sh` if you use that script.

[Bundled scripts](scripts.md) lists every script, and [Post-processing with exspec](exspec.md) says what `exspec`
does.
