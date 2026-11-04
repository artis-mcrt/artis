# Post-processing with exspec

:::{note}
`exspec` is optional. `sn3d` writes the light curves and the spectra itself, so a standard simulation needs no
post-processing.
:::

`make` builds `exspec` together with `sn3d`. It combines the packet files of all ranks into spectra and light curves,
with the same code that `sn3d` uses, so it gives the same files. Two situations need it:

- You want the files again with other settings, e.g. after a change of the frequency range or of the `MNUBINS` bin
  count of `spectrum_lightcurve.cc`. You then build the code again and run `exspec`, and you need no new simulation.
- The simulation stopped before its last requested timestep. `sn3d` writes the light curves and the spectra at each
  timestep, but it writes `emission.out`, `emissiontrue.out`, `absorption.out`, and the direction-resolved files
  only at the last requested timestep. `scripts/run-exspec-if-needed.sh` tests for `emission.out`, and it runs
  `exspec` in that case.

Run it in the model folder with any number of ranks from one up to the number of packet files:

```sh
mpirun -np 8 artis/exspec
```

`exspec` reads the same `input.txt`, model, and atomic data files as `sn3d`, so it must run in the same folder. The
`nprocs_exspec` line of `input.txt` gives the number of packet files, which is the number of `sn3d` ranks and not the
number of `exspec` ranks. Each `exspec` rank reads a block of the packet files, and the ranks of one node share the
spectra arrays. Only rank 0 writes a log, to `exspec.txt`. `exspec` changes no input file.

It writes `light_curve.out`, `spec.out`, `emission.out`, `emissiontrue.out`, and `absorption.out`. It also writes
`gamma_light_curve.out` and `gamma_spec.out` when `KEEP_ESCAPED_GAMMAS` is set, and `specpol.out`, `emissionpol.out`,
and `absorptionpol.out` when `POL_ON` is set. Direction-resolved versions of these go in the `speclc_angle_res`
folder. `scripts/mergeangleres.py` merges the files of the direction bins into `light_curve_res.out`, `spec_res.out`,
and `specpol_res.out`.
