# Post-processing with exspec

As well as `sn3d`, `make` builds `exspec`, which combines the packet files from all ranks into spectra and light
curves. `sn3d` writes the same files itself, so `exspec` is necessary only to make them again from the packet files,
e.g. after a change of the frequency range or of the `MNUBINS` bin count of `spectrum_lightcurve.cc`.

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
