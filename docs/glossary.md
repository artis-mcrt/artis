# Glossary

The pages of this site use these terms. [How a simulation works](simulation.md) puts them together.

:::{glossary}
estimator
  A sum that the packets build while they move through a cell, e.g. the mean intensity of the radiation field or a
  photoionisation rate. Each rank sums the contributions of its own packets, and `sn3d` then adds the sums of all
  ranks. The grid update of the next timestep reads them.

gamma packet
  A packet that holds the energy of its pellet as gamma rays. It gives that energy to the gas by Compton scattering,
  by photoelectric absorption, or by pair production. It can also escape.

job folder
  The subfolder `job_from_ts<n>` that `sn3d` writes in one Slurm job, where n is the first timestep of that job. It
  holds the log file of each rank and thread. See [Output files](outputfiles.md).

k-packet
  A packet of thermal kinetic energy. The code adds the cooling rate of each process of the cell, which are the
  free-free emission, the free-bound emission, the collisional excitation, and the collisional ionisation. It then
  samples one of these channels, which converts the k-packet into radiation or into a macro atom. In an optically thick
  cell, the k-packet emits an r-packet from a Planck function instead, and it samples no channel. A k-packet also loses
  energy, because the diffusion step takes the adiabatic loss of the radiation, and the thermal balance of the cell
  gives the share that becomes radiation.

macro atom
  The state of a packet between an absorption and the next emission. The macro atom makes internal radiative and
  collisional transitions between the levels and the ions, until it emits radiation or becomes a k-packet. Lucy
  (2002), A&A, 384, 725-735, doi:10.1051/0004-6361:20011756 and Lucy (2003), A&A, 403, 261-275,
  doi:10.1051/0004-6361:20030357 describe the scheme.

model folder
  The top-level folder of one simulation. It holds the input files, the restart files, and the job folders.

model grid cell
  A cell of the input model. `model.txt` gives its density and its composition. More than one propagation cell can
  share one model grid cell. The code calls its index `modelgridindex`, or `mgi`.

packet
  An indivisible energy packet. The code never splits a packet, and `sn3d` follows it until it escapes or the
  simulation ends. Some events change the comoving energy of a packet, e.g. the adiabatic loss of a k-packet. Lucy
  (2002), A&A, 384, 725-735, doi:10.1051/0004-6361:20011756 describes the method.

pellet
  The first state of a packet. Each pellet holds the same energy by default, which is the decay energy of the whole
  model between the first and the last time of the simulation, divided by the number of packets, and the decay time of
  each pellet is a sample of the decay curve ([UNIFORM_PELLET_ENERGIES](options.md#uniform_pellet_energies)). The pellet
  moves with the ejecta until the time of its decay. The decay then makes a gamma packet, or a packet of positrons,
  electrons, alpha particles, or fission fragments.

preset
  One of the `artisoptions_*.h` files. Each preset holds the compile-time options for one type of model. See
  [Option presets](presets.md).

propagation cell
  A cell of the grid that the packets move through. The code calls its index `cellindex`. The grid can hold more
  propagation cells than model grid cells, e.g. a 3D Cartesian grid over a 1D model.

r-packet
  A packet of ultraviolet, optical, or infrared radiation. It moves until it reaches a line, a continuum absorption,
  or an electron scattering. The r-packets that escape make the light curves and the spectra.

rank
  One MPI process. Each rank keeps its own share of the packets, and it updates its own range of the model grid cells.
  See [Running sn3d](running.md).

Sobolev approximation
  The approximation that a packet meets a line at one point. The expansion of the ejecta shifts the comoving frequency
  of the packet through the line there. The default mode of ARTIS gives each line its own Sobolev optical depth.

timestep
  One step of the time grid. `sn3d` updates the grid and then moves every packet to the end of the timestep.
  `timesteps.out` gives the start, the midpoint, and the width of each timestep.

virtual packet
  A packet that the code sends towards an observer direction at an emission or a scattering event of a real packet. It
  gives the spectrum of that direction with less noise. Bulla, Sim & Kromer (2015), MNRAS, 450, 967-981,
  doi:10.1093/mnras/stv657 describe the scheme.
:::
