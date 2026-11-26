# How a simulation works

`sn3d` follows energy packets through ejecta in {term}`homologous expansion`. This page describes the method.
[Running sn3d](running.md) describes how you start a run, and the [Glossary](glossary.md) gives each term.

## The packets

ARTIS uses indivisible energy packets (Lucy (2002), A&A, 384, 725-735, doi:10.1051/0004-6361:20011756). The code
never splits a packet into smaller packets, and it follows each packet until the packet escapes or the simulation
ends. A packet changes its type as its energy changes its form:

- each {term}`rank` samples the decay energy of the whole model. With [INITIAL_PACKETS_ON](options.md#initial_packets_on),
  the energy also includes decays before the first timestep and the initial energy from the model. By default, each
  {term}`pellet` receives an equal share and moves with the ejecta until its release time.
- the decay makes a {term}`gamma packet`, or a packet of positrons, electrons, alpha particles, or fission fragments;
- a gamma packet gives its energy to the gas by Compton scattering, by photoelectric absorption, or by pair
  production;
- with the Spencer-Fano solver, a fast lepton can heat, ionise, or excite the gas. The solver gives the
  probabilities for its packet. In a cell with {term}`grey opacity`, or with another [NT_SCHEME](options.md#nt_scheme),
  all deposited lepton energy heats the gas. Deposited energy from alpha particles and fission fragments always
  heats the gas.
- a {term}`k-packet` holds thermal energy. A cooling channel of the cell converts it into radiation or into a {term}`macro atom`. In
  an optically thick cell, the k-packet emits an {term}`r-packet` from a Planck function instead;
- a macro atom makes internal transitions until it emits radiation or becomes a k-packet again;
- an r-packet is ultraviolet, optical, or infrared radiation. It moves until it reaches a line, a continuum
  absorption, or an electron scattering.

The diagram gives the changes of the packet type. A gamma packet and an r-packet can also escape:

```{mermaid}
flowchart TD
    pellet["pellet"] -->|"decay"| gamma["gamma packet"]
    pellet -->|"decay, with transport"| particle["packet of positrons, electrons, or alpha particles"]
    pellet -->|"decay, instant deposition"| deposited["deposited energy"]
    gamma -->|"Compton, photoelectric, pair production"| deposited
    gamma -->|"Compton scattering, pair annihilation"| gamma
    particle -->|"energy loss"| deposited
    deposited -->|"heating"| kpkt["k-packet"]
    deposited -->|"ionisation, excitation"| ma["macro atom"]
    kpkt -->|"collisional cooling"| ma
    ma -->|"collisional"| kpkt
    kpkt -->|"free-free, free-bound"| rpkt["r-packet"]
    ma -->|"radiative"| rpkt
    rpkt -->|"line, bound-free"| ma
    rpkt -->|"free-free, bound-free"| kpkt
    rpkt -->|"leaves the grid"| escape["escaped packet"]
```

An r-packet that leaves the grid escapes. The escaped packets make the light curves and the spectra. With
[VPKT_ON](options.md#vpkt_on), the code also sends {term}`virtual packets <virtual packet>` towards fixed observer directions, which gives the
spectrum of each direction with less noise.

## The two grids

The input model gives the {term}`model grid cells <model grid cell>`. Each one has a density, a composition, and its own solution of the gas
state. The packets move through the {term}`propagation cells <propagation cell>`, and more than one propagation cell can share one model grid
cell. A 1D model on a 3D Cartesian grid is the usual example (see
[GRID_TYPE_OVERRIDE](options.md#grid_type_override)). The arrays of the gas state hold only the model grid cells that
contain matter, because an empty cell needs no solution.

## The timestep loop

`sn3d` repeats these steps for each {term}`timestep`:

```{mermaid}
flowchart TD
    update["1. grid update of the cells of each rank"] --> exchange["2. the ranks exchange the cell state"]
    exchange --> restart["3. restart files, from the second timestep of a job"]
    restart --> wall{"-w: enough wall time?"}
    wall -->|"no"| stop["job stops with RESTART_NEEDED"]
    wall -->|"yes"| propagate["4. propagation of the packets of each rank"]
    propagate --> estimators["5. the ranks sum the estimators"]
    estimators --> output["6. output of the timestep"]
    output --> next{"last requested timestep?"}
    next -->|"no"| update
    next -->|"yes"| packets["packet files, then the job ends"]
```

1. **The grid update.** Each rank solves the state of its own model grid cells. The first update of a new run uses
   the initial temperatures and the Saha ionisation balance. Later updates use the {term}`estimators <estimator>`
   from the previous timestep to calculate the temperature, the populations, and the opacities.
2. **The exchange of the cell state.** The ranks send their cells to the other ranks. Every rank then holds the whole
   grid, because a packet of any rank can enter any cell.
3. **The {term}`restart files <restart file>`.** From the second timestep of a job, `sn3d` writes the grid and the packets and updates
   `input.txt`. A run with the `-w` option then tests the remaining wall time, and it stops cleanly when too little
   time remains.
4. **The propagation.** Each rank moves its own packets to the end of the timestep.
5. **The estimators.** `sn3d` adds the estimators of all ranks, so that each rank holds the complete sums.
6. **The output.** `sn3d` writes the deposition, the light curves, and the spectra. At the last requested timestep it
   also writes the packet files.

A run that starts at a later timestep reads the restart files of the earlier run first. The restart files must come
from the same binary, because they hold the arrays of that build (see [Output files](outputfiles.md)).

## The estimators

A packet adds to the estimators of each cell that it enters, e.g. to the mean intensity of the radiation field and to
the photoionisation rates. The packets of all ranks enter the same cells, so the sums of one rank are incomplete
until step 5 above. The grid update of the next timestep then reads them.

Outside an {term}`LTE` update, a cell that used detailed transport in the previous timestep gets a dilute blackbody
fit to its radiation field. The code can also fit each frequency bin. The radiation field gives photoionisation rates
and bound-free heating. For a continuum with an active detailed estimator, the packet paths give the photoionisation
rate directly. The heating of that continuum still uses the radiation field. See
[the radiation field model](methods.md#the-radiation-field-model).

[Plan a simulation](planning.md#the-ranks) describes the work and the memory of each rank.
