# Post-processing with exspec

:::{note}
`exspec` is optional. `sn3d` writes the light curves and the spectra itself, so a standard simulation needs no
post-processing.
:::

`make` builds `exspec` together with `sn3d`. It combines the packet files of all {term}`ranks <rank>` into spectra and light curves,
with the same code that `sn3d` uses, so it gives the same files. It reads the packet files of `packets/`, which `sn3d`
writes at the last requested {term}`timestep`. A simulation that stopped earlier has no packet files, and `exspec` then gives
nothing. Two situations need `exspec`:

- You need the files again with other settings, e.g. after a change of the frequency range or of the `MNUBINS` bin
  count of `spectrum_lightcurve.cc`. You then build the code again and run `exspec`, and you need no new simulation.
- The {term}`model folder` holds the packet files but no spectra, e.g. from a run of an older version of ARTIS.
  `scripts/run-exspec-if-needed.sh` runs `exspec` when `emission.out` is absent.

Run it in the model folder with any number of ranks from one up to the number of packet files:

```sh
mpirun -np 8 artis/exspec
```

`exspec` reads the same `input.txt`, model, and atomic data files as `sn3d`, so it must run in the same folder. The
`nprocs_exspec` line of `input.txt` gives the number of packet files, which is the number of `sn3d` ranks and not the
number of `exspec` ranks. Each `exspec` rank reads a block of the packet files, and the ranks of one node share the
spectra arrays. Only rank 0 writes a log, to `exspec.txt`. `exspec` changes no input file.

It writes `light_curve.out`, `spec.out`, `emission.out`, `emissiontrue.out`, and `absorption.out`. A run with
`KEEP_ESCAPED_GAMMAS` also gets `gamma_light_curve.out` and `gamma_spec.out`, and a run with `POL_ON` also gets
`specpol.out`, `emissionpol.out`, and `absorptionpol.out`. A 2D or a 3D model gives a direction-resolved version of
each of these files in the `speclc_angle_res` folder, but the gamma files have no such form.
`scripts/mergeangleres.py` merges the files of the direction bins into `light_curve_res.out`, `spec_res.out`, and
`specpol_res.out`.
