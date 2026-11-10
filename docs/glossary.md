# Glossary

The pages of this site use these terms. [How a simulation works](simulation.md) puts them together.

:::{glossary}
estimator
  A sum that the packets build while they move through a cell, e.g. the mean intensity of the radiation field or a
  photoionisation rate. Each rank sums the contributions of its own packets, and `sn3d` then adds the sums of all
  ranks. The grid update of the next timestep reads them.

expansion opacity
  A binned opacity that sums the lines of a wavelength bin, with a weight from the Sobolev optical depth of each line.
  See [BINNED_OPACITY_LINE_WEIGHT](options.md#binned_opacity_line_weight).

gamma packet
  A packet that holds the energy of its pellet as gamma rays. It gives that energy to the gas by Compton scattering,
  by photoelectric absorption, or by pair production. It can also escape.

grey opacity
  An opacity with no dependence on the frequency. In the first `num_grey_timesteps` timesteps, a thick cell gets a grey
  opacity (see [Plan a simulation](planning.md#the-first-timesteps)).

homologous expansion
  The motion of ejecta where each mass element keeps a constant velocity. The position is then the velocity multiplied
  by the time, and the density decreases with the time to the power of -3.

job folder
  The subfolder `job_from_ts<n>` that `sn3d` writes in one Slurm job, where n is the first timestep of that job. It
  holds the log file of each rank and thread. See [Output files](outputfiles.md).

k-packet
  A packet of thermal kinetic energy. A cooling process of the cell converts it into radiation or into a macro atom
  (see [Methods and references](methods.md#the-macro-atom-and-the-k-packets)).

LTE
  Local thermodynamic equilibrium. The Saha and the Boltzmann relations give the populations.

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

NLTE
  A state that is not in LTE. A rate matrix gives the populations of the NLTE levels.

packet
  An indivisible energy packet. The code never splits a packet, and `sn3d` follows it until it escapes or the
  simulation ends. Some events change the comoving energy of a packet, e.g. the adiabatic loss of a k-packet. Lucy
  (2002), A&A, 384, 725-735, doi:10.1051/0004-6361:20011756 describes the method.

pellet
  The first state of a packet. The pellet moves with the ejecta until the time of its decay. The decay then makes a
  gamma packet, or a packet of positrons, electrons, alpha particles, or fission fragments.

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

restart file
  A file `gridsave_ts*.tmp` or `packets_*_ts*.tmp` that holds the state of a run at the start of a timestep. A resumed
  run must use the binary that wrote the file.

Sobolev approximation
  The approximation that a packet meets a line at one point. The expansion of the ejecta shifts the comoving frequency
  of the packet through the line there. The default mode of ARTIS gives each line its own Sobolev optical depth.

Spencer-Fano solution
  The energy spectrum of the fast electrons of a cell. It gives the fractions of the deposited energy that go to the
  heating, the ionisation, and the excitation (see [Methods and references](methods.md#the-non-thermal-electrons)).

superlevel
  One row of the NLTE rate matrix that holds the levels of an ion above its NLTE levels.

thermalisation
  The conversion of energy into thermal energy of the gas, e.g. of the energy of the decay particles.

timestep
  One step of the time grid. `sn3d` updates the grid and then moves every packet to the end of the timestep.
  `timesteps.out` gives the start, the midpoint, and the width of each timestep.

virtual packet
  A packet that the code sends towards an observer direction at an emission or a scattering event of a real packet. It
  gives the spectrum of that direction with less noise. Bulla, Sim & Kromer (2015), MNRAS, 450, 967-981,
  doi:10.1093/mnras/stv657 describe the scheme.

Ye
  The electron fraction, which is the number of electrons per baryon. `Ye.txt` or a column of `model.txt` gives it.
:::
