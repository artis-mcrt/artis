# Analyse the results

[artistools](https://github.com/artis-mcrt/artistools) reads the output files of a {term}`model folder`. The
[Quickstart](quickstart.md) says how to install it. Each command takes the model folder as its last argument, and
`.` is the current folder. Each plot command writes a PDF file into the current folder. `--help` after a command lists
its options.

## Start with the checks

```sh
artistools describe <model folder>
artistools check <model folder>
artistools logsummary <model folder>
```

`describe` gives the grid, the {term}`timesteps <timestep>`, and the options of the run. It also names the command that reads each output
file. `check` gives one row for each check of the results. `logsummary` counts each kind of warning in the logs of
{term}`rank` 0.

## Commands for each task

| Task | Command |
| --- | --- |
| the bolometric light curve | `artistools plotlightcurves <model folder>` |
| the light curves of some filters | `artistools plotlightcurves --filter B V <model folder>` |
| the deposition of the decay products | `artistools plotlightcurves --deposition gamma betaplus <model folder>` |
| the spectrum at one time | `artistools plotspectra -t <days> <model folder>` |
| the spectrum by emission and absorption process | `artistools plotspectra -t <days> --emissionabsorption <model folder>` |
| the emission of each element | `artistools plotspectra -t <days> --showemission --groupby element <model folder>` |
| the spectrum from the packet files | `artistools plotspectra -t <days> --frompackets <model folder>` |
| the spectra of each timestep as text | `artistools writespectra <model folder>` |
| an estimator of one cell against the timestep | `artistools plotestimators Te -x timestep --cell <cell> <model folder>` |
| estimators of all cells at one time | `artistools plotestimators -t <days> -p Te nne <model folder>` |
| the populations of one ion | `artistools plotnltepops Fe --ion-stages 2 --cell <cell> -t <days> <model folder>` |
| the radiation field of a cell | `artistools plotradfield --cell <cell> -t <days> <model folder>` |
| the deposition rates of the cells | `artistools deposition <model folder>` |
| the time that each rank took | `artistools plotlogfiles <model folder>` |
| the core hours of all jobs | `artistools corehours <model folder>` |

`plotspectra` accepts only the times at which the whole model can have sent light to the observer. The row
`arrival times` of `artistools check` gives that range. `plotnltepops` needs a run with {term}`NLTE` levels.

## Direction-resolved results

A 2D or a 3D model gives the results of each observer direction bin. Each `--plotviewingangle` number is one bin:

```sh
artistools plotspectra -t <days> --plotviewingangle 0 <model folder>
artistools plotlightcurves --plotviewingangle 0 10 <model folder>
artistools plotspherical <model folder>
```

`plotspherical` draws a direction map from the escaped packets.

## Is the packet count sufficient?

```sh
artistools check --frompackets <model folder>
```

The row `packet noise` gives the largest relative noise of a timestep and the noise at the peak. For each timestep,
that noise is the square root of the sum of the squared packet energies, divided by the sum of the packet energies.
The check reads the packet files, so it is slow. With packets of equal energy, 100 times more packets give 10 times
less noise (see [Plan a simulation](planning.md#the-packets)).

The {term}`estimators <estimator>` have their own noise, and the check above does not measure it. Run the model again with a different
packet count, and compare T_e and the ionisation of the cells. The results are sufficient when the difference is small
for your purpose.

## Is the timestep count sufficient?

Run the model again with more timesteps. Compare the light curve, the spectra, and T_e of the cells at the same times.
Also examine these points:

- Plot T_e of a cell against the timestep, and look at the first timesteps after the {term}`LTE` and the grey timesteps. The
  gas state needs some timesteps to become stable (see [Plan a simulation](planning.md#the-first-timesteps)).
- A warning `timestep ... has a width of ... times its midpoint time` in the log says that a timestep is too wide for
  the time-dependent solutions.
- The row of the temperature limits of `artistools check` gives the cells at `MINTEMP` or `MAXTEMP`.
