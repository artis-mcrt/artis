# Post-processing with exspec

:::{note}
`exspec` is optional. `sn3d` writes the light curves and the spectra itself, so a standard simulation needs no
post-processing.
:::

`make` builds `exspec` together with `sn3d`. It combines the packet files of all {term}`ranks <rank>` into spectra and light curves,
with the same code that `sn3d` uses, so it gives the same files for a simulation that ran to the last timestep of the
model. `exspec` writes all timesteps of the model. It reads the packet files of `packets/`, which `sn3d`
writes at the last requested {term}`timestep`. A simulation that stopped earlier has no packet files, and `exspec` then stops
with an error. Two situations need `exspec`:

- You need the files again with other settings, e.g. after a change of the frequency range. You then build the code
  again and run `exspec`, and you need no new simulation.
- The {term}`model folder` holds the packet files but no spectra, e.g. from a run of an older version of ARTIS.
  `scripts/run-exspec-if-needed.sh` runs `exspec` when `emission.out` is absent.

The requested end of the timestep range limits the latest escape time of the packets. The packet files can also
hold escapes from earlier jobs. The spectra and the observer light curve use the arrival time, which includes the
light travel time. A far-side packet can escape before the requested end and arrive after it. For example, an
escape at 9.5 days can contribute to an arrival bin at 10.5 days.

A near-side packet that escapes after the requested end can arrive before it. The packet files cover escape
times through the requested end only. Arrival bins that are incomplete therefore start before the requested end. The complete range of arrival
times ends earlier than the requested end.

After a partial run, later arrival bins can also hold nonzero flux. An incomplete bin can hold nonzero flux.
Extend the simulation to obtain the contributions of the later escapes.

Run it in the model folder with any number of ranks from one up to the number of packet files:

```sh
mpirun -np 8 artis/exspec
```

`exspec` reads the same `input.txt`, model, and atomic data files as `sn3d`, so it must run in the same folder. The
`nprocs_exspec` line of `input.txt` gives the number of packet files, which is the number of `sn3d` ranks and not the
number of `exspec` ranks. Rank 0 writes the log to `exspec.txt`. `exspec` changes no input file.

It writes `light_curve.out`, `spec.out`, `emission.out`, `emissiontrue.out`, and `absorption.out`. A run with
`KEEP_ESCAPED_GAMMAS` also gets `gamma_light_curve.out` and `gamma_spec.out`, and a run with `POL_ON` also gets
`specpol.out`, `emissionpol.out`, and `absorptionpol.out`. A 2D or a 3D model also gives the direction-resolved
files. [Output files](outputfiles.md#the-format-of-the-files) gives the format of each file.
