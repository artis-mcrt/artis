# Input files

These files go in the model folder, which should always contain the ARTIS source folder (or a symlink to it) named
`artis`. The physics data that comes with the code is then found automatically in `artis/data` and needs no copy (see
[Atomic and model data](data.md)).

A build with libzstd also reads a text input file in zstd compressed form, e.g. `model.txt.zst`,
`transitiondata.txt.zst`, or `packets/packets00_0000.out.zst` for `exspec`. The plain file has priority when both
exist.

`sn3d` writes the final packet files of a job into the folder `packets/`, and the virtual packet files into
`vpackets/`, `vspecpol/`, and `vpkt_grid/`. `exspec` reads the packet files there.

## input.txt

Run-time configuration with:

- the random number seed, which must be a fixed value for reproducible runs;
- the number of timesteps;
- the first and last timestep of this job. Normally the last timestep should be set to the number of timesteps. Long
  simulations are split over resubmitted jobs by the wall-time mechanism (`sn3d -w`), which advances the start timestep
  on each restart but never the last one. An earlier last timestep makes the simulation stop there, e.g. to inspect or
  post-process partial results, until `input.txt` is edited to continue. The last restart files of such a run hold the
  start of its last timestep. The continuation job therefore propagates that timestep again and writes its spectra and
  light curves again;
- the start and end time in days;
- whether the run continues from the restart files of a previous job;
- the number of pure LTE timesteps;
- the optically-thick condition that switches cells to a grey opacity treatment;
- `nprocs_exspec`: the number of packet files that `exspec` will read, one for each `sn3d` rank. `sn3d` sets this line
  itself. At the start of a new simulation, `sn3d` writes `input.txt` again with a comment on each line and with this
  value, and it keeps the same content in `input-newrun.txt`.

The format is positional, so every line must be present and in order, including the lines that are no longer used. See
[tests/classicmode_3d_inputfiles/input-newrun.txt](https://github.com/artis-mcrt/artis/blob/main/tests/classicmode_3d_inputfiles/input-newrun.txt)
for a complete example with a comment on every line.

If `input.txt` is absent at the start of a run, rank 0 restores it from `input-newrun.txt`.

## vpkt.txt

Required when virtual packets are enabled (`VPKT_ON` in `artisoptions.h`). It sets the observer directions, the number
of spectra per observer, and the time and wavelength ranges to record. See
[tests/classicmode_3d_inputfiles/vpkt.txt](https://github.com/artis-mcrt/artis/blob/main/tests/classicmode_3d_inputfiles/vpkt.txt)
for an example. `artistools makevpktinput` writes the file for a model.

## model.txt

Grid parameters, cell densities, and nuclear composition. The optional column `q` holds the trapped radiation energy
per mass [erg/g] at the snapshot time of the model. The value must already include the adiabatic losses before the
snapshot time. `INITIAL_PACKETS_ON` uses it for the initial temperature and for the initial packets.

The cell ids start at 0 or at 1. `sn3d` finds the first id from the first cell of `model.txt`. `abundances.txt` and
`Ye.txt` must use the same ids. The log messages always give the 0-based cell index.

## abundances.txt

Required file with the per-cell elemental mass fractions, which include elements whose isotopic abundances are not
given in `model.txt`.

For a 1D or a 2D model, `sn3d` divides the values of each cell by their sum, so the values can also be proportional to
the mass fractions, e.g. element densities. For a 3D model, `sn3d` uses the values without normalisation. A 3D cell
whose values do not sum to 1 within 2 percent gets a warning in the log.

## adata.txt

One block per ion consisting of:

Header:<br/>
`[atomic number] [ionisation stage] [level count] [ionisation energy in eV]`

Level count rows of:<br/>
`[level index] [energy in eV] [statistical weight] [level transition count] [level name or description]`

The level indices start at 0 or at 1. `sn3d` finds the first index from the first level of `adata.txt`.
`transitiondata.txt` and `phixsdata_v2.txt` must use the same numbering. `phixsdata.txt` always numbers the levels
from 1. `sn3d` finds the first index of `autoion.txt` separately: the numbering starts at 0 if a level number in the
file is 0, and at 1 otherwise. The highest autoionising level of each ion must be the highest level of the ion. The log
messages always give the 0-based level index.

## transitiondata.txt

One block per ion consisting of:

Header:<br/>
`[atomic number] [ionisation stage] [transition count]`

Transition count rows of:<br/>
`[lower level index] [upper level index] [A value] [collision strength] [forbidden flag]`

In case the collision strength is not available, it will be -1.0 for permitted and -2.0 for forbidden transitions.

A legacy four-column format (`[transition index] [lower level index] [upper level index] [A value]`, with no collision
strength or forbidden columns) is also accepted. The format is detected from the column count of the first row of each
block.

## compositiondata.txt

It sets a filter on the elements, ion stages, and energy levels that the code reads from the atomic data files. The
code ignores the text to the right of a `#` character.

## phixsdata.txt or phixsdata_v2.txt

Photoionisation cross sections in v1 format (arbitrary energy tables) or v2 format (regularly spaced energy grid).

## gamma_*.txt

Gamma decay spectra containing energies [MeV] and probabilities.
