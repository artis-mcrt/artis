# Make your own model

This page makes a small 1D model from a density profile, checks it, and runs it. [Input files](inputfiles.md) gives
every rule of the formats. [Atomic and model data](data.md) names the converters of artistools for the output of a
hydrodynamic simulation or of a nuclear network.

The model of this page has 1 solar mass of iron-group ejecta in 20 shells. The density is exponential in the velocity,
and the inner shells hold Ni56. It is a test of the method and not a model of a real transient.

## 1. Make the model folder

Run these commands in a clone of `https://github.com/artis-mcrt/artis.git`, e.g. the clone of the Quickstart. The
commit of that clone must be on GitHub, because the script clones from the URL of `origin`.

The script makes the model folder and clones the branch `release` into its subfolder `artis`, so that the options of
section 4 belong to this model only. The `git` command then checks out the commit of your clone in place of `release`. The `cp` command copies the template of `input.txt` of
that commit, because the lines of `input.txt` can change between versions. The last two commands download the atomic
data of the tests and unpack it:

```sh
scripts/new-model-folder.sh ~/ownmodel
git -C ~/ownmodel/artis checkout --detach "$(git rev-parse HEAD)"
cd ~/ownmodel
cp artis/scripts/input-template.txt input.txt
curl -fLO "https://github.com/artis-mcrt/artis/releases/download/$(cat artis/tests/atomic_data_release.txt)/atomicdata_feconi.tar.xz"
tar -xf atomicdata_feconi.tar.xz
```

`artis/tests/atomic_data_release.txt` names the release that holds the archive.

The archive supplies atomic data for iron (Fe), cobalt (Co), and nickel (Ni).
This model uses nickel for the radioactive material and iron for the stable material.
The decay of Ni56 also produces cobalt.

## 2. Write the model files

Put this script into the {term}`model folder` as `make_model.py`:

```python
"""Write model.txt and abundances.txt of a 1D model with an exponential density profile."""

import math
from pathlib import Path

MSUN = 1.989e33  # g
DAY = 86400.0  # s
KM = 1.0e5  # cm

t_model_days = 1.0
ncells = 20
v_max_kms = 15000.0
v_scale_kms = 3000.0  # the e-folding velocity of the density
mass_total_msun = 1.0
v_ni56_max_kms = 6000.0  # the Ni56 is inside this velocity
x_ni56_inner = 0.7

t_model_s = t_model_days * DAY
v_outer_kms = [v_max_kms * (i + 1) / ncells for i in range(ncells)]
v_inner_kms = [0.0, *v_outer_kms[:-1]]

# volume of each shell at t_model, and an unscaled density at the middle velocity of the shell
volume_cm3 = [
    4.0 / 3.0 * math.pi * ((v_o * KM * t_model_s) ** 3 - (v_i * KM * t_model_s) ** 3)
    for v_i, v_o in zip(v_inner_kms, v_outer_kms, strict=True)
]
rho_unscaled = [
    math.exp(-0.5 * (v_i + v_o) / v_scale_kms) for v_i, v_o in zip(v_inner_kms, v_outer_kms, strict=True)
]

# scale the density to the total mass
mass_unscaled = sum(r * v for r, v in zip(rho_unscaled, volume_cm3, strict=True))
rho_g_cm3 = [r * mass_total_msun * MSUN / mass_unscaled for r in rho_unscaled]

x_ni56 = [x_ni56_inner if v_o <= v_ni56_max_kms else 0.0 for v_o in v_outer_kms]

with Path("model.txt").open("w", encoding="utf-8") as model_file:
    model_file.write(f"{ncells}\n{t_model_days}\n")
    model_file.write("#inputcellid vel_r_max_kmps logrho X_Fegroup X_Ni56 X_Co56 X_Fe52 X_Cr48\n")
    for cell_id in range(ncells):
        logrho = math.log10(rho_g_cm3[cell_id])
        model_file.write(f"{cell_id} {v_outer_kms[cell_id]} {logrho:.6f} 1.0 {x_ni56[cell_id]} 0.0 0.0 0.0\n")

with Path("abundances.txt").open("w", encoding="utf-8") as abundances_file:
    for cell_id in range(ncells):
        x_by_atomic_number = [0.0] * 28  # Z = 1 to 28
        x_by_atomic_number[28 - 1] = x_ni56[cell_id]
        x_by_atomic_number[26 - 1] = 1.0 - x_ni56[cell_id]
        abundances_file.write(f"{cell_id} " + " ".join(str(x) for x in x_by_atomic_number) + "\n")

mass_ni56_msun = sum(
    r * v * x for r, v, x in zip(rho_g_cm3, volume_cm3, x_ni56, strict=True)
) / MSUN
print(f"total mass {mass_total_msun:.3f} Msun, Ni56 mass {mass_ni56_msun:.3f} Msun")
```

Then run it:

```sh
python3 make_model.py
```

The script gives this line:

```text
total mass 1.000 Msun, Ni56 mass 0.260 Msun
```

The script uses format version 1. It writes two files.

- `model.txt` has the number of cells, the snapshot time in days, a header line, and one row for each shell.
- `abundances.txt` has one row for each cell. The values start at hydrogen, so Fe is value 26 and Ni is value 28.

In the first shell, `X_Ni56` is 0.7. The nickel entry in `abundances.txt` also contains 0.7.
These values describe the same mass. The iron entry contains the remaining 0.3.
`X_Fegroup` is 1.0 because all the mass in this shell belongs to iron-group elements.

## 3. Check the model

```sh
artistools inputmodel describe .
```

Compare the output with the values for this model.

| Quantity in the output | Expected value |
| --- | --- |
| Snapshot time | 1 day |
| Number of shells | 20 |
| `vmax` | 15,000 km/s |
| `M_tot_rho` (sum of density times volume) | 1.000 solar mass |
| `M_tot_elem` (sum of the element masses) | 1.000 solar mass |
| Fe and Ni at the snapshot time | 0.740 and 0.260 solar masses |

For a model made from a hydrodynamic simulation, compare these masses with the masses in the original snapshot.
Use the density and the shell volumes at the same time. In homologous expansion, multiplying the time by two
multiplies each volume by eight and divides each density by eight.

`artistools plotdensity .` plots the enclosed mass and the mass per unit velocity against the velocity.

## 4. Select the options and the times

Select the {term}`preset` that is nearest to the model (see [Option presets](presets.md)), and set a small packet count:

```sh
cd artis
ln -s artisoptions_classic.h artisoptions.h
```

Change `NUM_PACKETS` in `artisoptions.h` to `1e5`, and then build:

```sh
make
cd ..
```

`artisoptions.h` is a symlink, so `git -C artis restore artisoptions_classic.h` in the model folder removes the change.

Change these lines of `input.txt`, and keep the comment of each line:

```text
20                       #  1: ntimesteps ...
000 004                  #  2: timestep_start timestep_finish ...
2 30                     #  3: tmin_days tmax_days ...
3                        # 18: num_lte_timesteps
3.0 3                    # 19: optical_depth_is_thick num_grey_timesteps
```

`tmin_days` must be at or after the snapshot time of `model.txt`, or the run stops. This job runs the first four
{term}`timesteps <timestep>`. [Plan a simulation](planning.md) says what each value does.

## 5. Check the input and run

```sh
mpirun -np 4 artis/sn3d -c
```

The standard output has the line `input check: complete` when the check finds no error. Read the warnings of the
check:

```sh
grep warning inputcheck/output_0-0.txt
```

[Troubleshooting](troubleshooting.md#warnings-that-need-a-check-of-the-input) gives the warnings that need a change of
the input. A run with the Fe, Co, and Ni data of this page gives a warning that compares the levels of the ground term
of Co III with its photoionisation targets. The tests give the same warning.

Then start the run:

```sh
mpirun -np 4 artis/sn3d
```

Look for the lines `timestep 3 is complete` in `output_0-0.txt`. To run the remaining timesteps, set
`timestep_finish` in `input.txt` to 20 and start `sn3d` again (see [Plan a simulation](planning.md)).

## What a real model changes

- **The composition.** Add each element with a significant mass to `abundances.txt`, and use an atomic database that
  holds these elements (see [Atomic and model data](data.md)).
- **The nuclides.** Add a column `X_<nuclide>` to `model.txt` for each other radioactive nuclide, e.g. `X_Ni57` and
  `X_Co57`.
- **The packets and the shells.** Increase `NUM_PACKETS`, and use enough shells to resolve the density and the
  composition.
- **The geometry.** A multidimensional model can specify its velocity edges and its composition in one file.
  See the [3D spherical example](inputfiles.md#format-version-2).

Next: [Plan a simulation](planning.md).
