# How a simulation works

`sn3d` follows energy packets through ejecta in {term}`homologous expansion`. This page describes the method.
[Running sn3d](running.md) describes how you start a run, and the [Glossary](glossary.md) gives each term.

## The packets

ARTIS uses indivisible energy packets (Lucy (2002), A&A, 384, 725-735, doi:10.1051/0004-6361:20011756). The code
never splits a packet into smaller packets, and it follows each packet until the packet escapes or the simulation
ends. A packet changes its type as its energy changes its form:

- each {term}`pellet` holds the same energy, which is the decay energy of the whole model between the first and the last
  time of the simulation, divided by the number of packets of the {term}`rank`.
  A pellet moves with the ejecta until the time of its decay.
  [INITIAL_PACKETS_ON](options.md#initial_packets_on) adds the initial energy of the cells to that sum. With
  [UNIFORM_PELLET_ENERGIES](options.md#uniform_pellet_energies) false, the decay times are uniform instead, and the
  energy of each pellet follows the decay power at its own time;
- the decay makes a {term}`gamma packet`, or a packet of positrons, electrons, alpha particles, or fission fragments;
- a gamma packet gives its energy to the gas by Compton scattering, by photoelectric absorption, or by pair
  production;
- the energy of a fast lepton goes to the heating, the ionisation, and the excitation of the gas. The Spencer-Fano
  solver gives the three fractions (see [NT_SCHEME](options.md#nt_scheme)). The alpha particles and the fission
  fragments give all of their energy to the heating, so they need no such solution;
- a {term}`k-packet` holds thermal energy. A cooling channel of the cell converts it into radiation or into a {term}`macro atom`. In
  an optically thick cell, the k-packet emits an {term}`r-packet` from a Planck function instead. A
  k-packet also loses energy: the diffusion step takes the adiabatic loss of the radiation, and the thermal balance
  of the cell gives the share that becomes radiation;
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

1. **The grid update.** Each rank solves the state of its own model grid cells. The solution gives the electron
   temperature, the ion and level populations, and the opacities, and it uses the {term}`estimators <estimator>` of the previous
   timestep.
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

The radiation field model fits a diluted blackbody to the full spectrum of each cell, and it can also fit one to each
frequency bin (Shingles et al. (2020), MNRAS, 492, 2029-2043, section 2.2, doi:10.1093/mnras/stz3412). The rates of
the photoionisation and of the heating use these fits.

## The work of each rank

Each rank updates a range of the model grid cells, and `modelgridrankassignments.out` gives the range of each rank.
Each rank also keeps an equal share of [NUM_PACKETS](options.md#num_packets). The ranks on one node share the large
read-mostly arrays, e.g. the atomic data, so a node holds one copy of them. [Running sn3d](running.md) describes the
ranks, the threads, and the memory of a run.
