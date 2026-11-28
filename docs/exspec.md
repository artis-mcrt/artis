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
times through the requested end only. An arrival bin can therefore have nonzero flux and still lack contributions
from later escapes. Extend the simulation to include those contributions.

Run `exspec` in the model folder. Use one rank or more, up to the number of packet files.

```sh
mpirun -np 1 artis/exspec
```

`exspec` replaces the spectra and light curves in the model folder. Keep a copy of the original files if you need
to compare the two calculations.

Use the model and atomic data files from the simulation. The line indices in the packet files refer to the
original line list. Build `exspec` from the same source version as `sn3d`.
Keep the values of `POL_ON` and `SAMPLE_RPKT_EMISSION`. These options select columns in the packet files.
The `nprocs_exspec` line of `input.txt` gives the number of `sn3d` ranks that wrote the packet files.

:::{warning}
If `input.txt` and `input.txt.zst` are both absent, `exspec` copies a plain `input-newrun.txt` to `input.txt` when it
finds one.
Before you resume `sn3d`, [check the mode and timestep](troubleshooting.md#the-run-starts-again-from-the-first-timestep)
against your restart files. A new simulation [removes the previous output and restart files](outputfiles.md#a-new-simulation-removes-the-files-of-the-previous-one).
:::

Rank 0 writes the log to `exspec.txt`. It writes one line to the standard output at the start and one at the end.
`artis/exspec -h` gives the help.
`artis/exspec --version` gives the git version, the preset, and the compiler of the build. An unknown option or
argument stops `exspec` before it reads a file.

`exspec` writes its output in the model folder.

| Quantity | Files | Builds that write the files |
| --- | --- | --- |
| Light curve | `light_curve.out` | All builds |
| Spectrum | `spec.out` | All builds |
| Spectra by atomic process | `emission.out`, `emissiontrue.out`, `absorption.out` | All builds |
| Gamma rays | `gamma_light_curve.out`, `gamma_spec.out` | `KEEP_ESCAPED_GAMMAS` |
| Polarisation | `specpol.out`, `emissionpol.out`, `absorptionpol.out` | `POL_ON` |

[Output files](outputfiles.md#the-format-of-the-files) gives the columns and units.

## Read the spectra of each direction

For a 2D or a 3D input model, `sn3d` and `exspec` write the files of each direction bin in `speclc_angle_res/`.
You can read a single bin directly. Run the merge script in the model folder to prepare the files for artistools.

```sh
python3 artis/scripts/mergeangleres.py
```

Run the command again if `sn3d` or `exspec` writes new files for the direction bins.
The [output file reference](outputfiles.md#the-direction-resolved-files) gives the bin numbers, file layout, and flux normalisation.

For a spectrum along a direction from `vpkt.txt`, use the [virtual packet files](outputfiles.md#the-virtual-packet-files).
`sn3d` calculates those spectra during the simulation.
