# Post-processing with exspec

:::{note}
`exspec` is optional. `sn3d` writes the light curves and the spectra itself, so a standard simulation needs no
post-processing.
:::

`make` builds `exspec` together with `sn3d`. Both programs use the same calculation for the spectra and the light
curves. `sn3d` uses the packet values in memory and rounds them when it writes the packet files. `exspec` reads
the rounded values from all {term}`ranks <rank>`, so its spectra and light curves can differ from those of `sn3d`.

`exspec` writes all timesteps of the model. It reads the files in `packets/`, which `sn3d` writes at the last
requested {term}`timestep`. If the simulation stops earlier, `exspec` stops because it cannot find the packet files.

To change the frequency bins for the saved packets, rebuild `exspec` with the new frequency limits.
Run `exspec` to write the spectra again. [NU_MIN_R](options.md#nu_min_r) and [NU_MAX_R](options.md#nu_max_r)
also set the frequency limits for blackbody emission in `sn3d`. To apply new limits to thermal emission,
rebuild `sn3d` with those limits. Start a new simulation.

An older {term}`model folder` can hold packet files but no spectra. `scripts/run-exspec-if-needed.sh` runs
`exspec` if the folder has neither `emission.out` nor `emission.out.zst`.

The requested end of the timestep range limits the latest escape time of the packets. The packet files can also
hold escapes from earlier jobs. The spectra and the observer light curve use the arrival time, which includes the
light travel time. A far-side packet can escape before the requested end and arrive after it. For example, an
escape at 9.5 days can contribute to an arrival bin at 10.5 days.

A near-side packet that escapes after the requested end can arrive before it. The packet files cover escape
times through the requested end only. Arrival bins that are incomplete therefore start before the requested end. The complete range of arrival
times ends earlier than the requested end.

After a partial run, later arrival bins can also hold nonzero flux. An incomplete bin can hold nonzero flux.
Extend the simulation to obtain the contributions of the later escapes.

Run `exspec` in the model folder. Use one rank or more, up to the number of packet files.

```sh
mpirun -np 1 artis/exspec
```

`exspec` reads the same `input.txt`, model, and atomic data files as `sn3d`.
The `nprocs_exspec` line of `input.txt` gives the number of `sn3d` ranks that wrote the packet files.

:::{warning}
If `input.txt` is absent, `exspec` copies `input-newrun.txt` to `input.txt`. This copy holds the parameters for a new
simulation. In this case, restore the `input.txt` for your restart files before you resume `sn3d`.
Check that `timestep_start` matches those files and that `simulation_continued_from_saved` is 1.
A new simulation [removes previous output and restart files](outputfiles.md#a-new-simulation-removes-the-files-of-the-previous-one).
:::

Rank 0 writes the log to `exspec.txt`. It writes one line to the standard output at the start and one at the end.
`artis/exspec -h` gives the help.
`artis/exspec --version` gives the git version, the preset, and the compiler of the build. An unknown option or
argument stops `exspec` before it reads a file.

It writes `light_curve.out`, `spec.out`, `emission.out`, `emissiontrue.out`, and `absorption.out`. A run with
`KEEP_ESCAPED_GAMMAS` also gets `gamma_light_curve.out` and `gamma_spec.out`, and a run with `POL_ON` also gets
`specpol.out`, `emissionpol.out`, and `absorptionpol.out`. A 2D or a 3D model also gives the direction-resolved
files. [Output files](outputfiles.md#the-format-of-the-files) gives the format of each file.
