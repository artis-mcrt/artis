# Output files

A simulation uses two types of folder:

- The {term}`model folder` is the top-level folder. It holds the input files, the {term}`restart files <restart file>`, and the files of the whole
  simulation.
- A {term}`job folder` is a subfolder of the model folder. `sn3d` writes one job folder in each Slurm job. `sn3d` names the
  folder from the start {term}`timestep` of the job, e.g. `job_from_ts0000` for a new simulation and `job_from_ts0008` for a
  job that resumes at timestep 8. Rank 0 writes the name of the folder to the standard output, so the log of a Slurm
  job names its folder.

A model folder after two jobs has this layout:

```text
mymodel/                              the model folder
├── artis/                            the source code and the programs
├── input.txt, input-newrun.txt       the run parameters
├── model.txt, abundances.txt         the model of the ejecta
├── adata.txt, transitiondata.txt, phixsdata_v2.txt, compositiondata.txt
├── spec.out.zst, light_curve.out.zst, deposition.out.zst, ...
├── output_0-0.txt                    a symlink to the log of rank 0
├── packets/                          the packet files of the last requested timestep
├── job_from_ts0000/                  the job folder of the first job
│   ├── output_0-0.txt, output_1-0.txt, ...
│   └── estimators_allranks.out.zst
└── job_from_ts0008/                  the job folder of the second job
```

During a run, the model folder also holds the restart files `gridsave_ts*.tmp` and `packets_*_ts*.tmp`.

A build with libzstd writes each `.out` file zstd compressed, e.g. `spec.out.zst`, so that no separate compression
step is necessary. The logs, `input.txt`, `artis.pid`, `syn_dir.txt`, and the restart files stay plain. Every file
gets zstd level 9.
[artistools](https://github.com/artis-mcrt/artistools) reads both forms.

The first two sections of this page list the files. [The format of the files](#the-format-of-the-files) gives the
columns and the units of each file.

## The files of each job

Each job writes these files into its job folder:

- `output_n-t.txt`: a log file for each MPI {term}`rank` n and each thread t of that rank. A run without threads writes
  only the files of thread 0.
- the {term}`estimators <estimator>` of each cell at each timestep: the temperatures, the ionisation, and the heating and cooling rates.
  The option `WRITE_COMBINED_ALLRANK_OUT_FILES` selects one of two forms, and every {term}`preset` sets the option:
    - with the option, `sn3d` writes the text of all ranks into one combined file, `estimators_allranks.out`;
    - without the option, each rank writes its own file, `estimators_nnnn.out`. A rank that updates no cell writes
      no file.

    The option gives the same choice for the nlte and the radfield files, e.g. `nlte_allranks.out`. With the constant
    `LOG_MACROATOM` in `macroatom.cc`, it also gives `macroatom_allranks.out`.
    `scripts/combine_allrank_out_files.py` makes the combined files from the files of the ranks, and it keeps the
    files of the ranks unless you give `--rm`. Run the script with uv, e.g.
    `uv run artis/scripts/combine_allrank_out_files.py`.
  The script needs Python 3.14 or a later version, and uv gets that version from the metadata of the script.

## The files of the whole simulation

A run writes these files into the model folder:

- `packets/packets00_nnnn.out`: the Monte Carlo packets from each rank, which `exspec` can turn into spectra and light
  curves again.
- the files of the {term}`virtual packets <virtual packet>`, which a run with `VPKT_ON` writes at the last requested timestep:
    - `vspecpol/vspecpol_nnnn.out`: the spectrum of each observer direction, with the Stokes parameters I, Q, and U;
    - `vpkt_grid/vpkt_grid_nnnn.out`: the emission over the velocity grid, when `vpkt.txt` switches that map on;
    - `vpackets/vpackets_nnnn.out`: one row for each contribution of a virtual packet, with `VPKT_WRITE_CONTRIBS`.
- `light_curve.out`, `spec.out`, and the other spectrum files that [Post-processing with exspec](exspec.md) lists.
  `sn3d` writes the light curves and spectra at each timestep. It writes the emission and absorption files at the
  last requested timestep. A 2D or a 3D model also gives the direction-resolved files there.
- `deposition.out`: the radioactive energy deposition rate as a function of time. The columns of the deposited
  luminosities come from the Monte Carlo estimators, also when `PARTICLE_THERMALISATION_SCHEME` is
  `INSTANTFULLDEPOSITION` and the heating of the cells uses the analytic emission rates. For a beta-plus decay, the
  `Qdot_ana_erg/s/g` column uses the energy Q_EC, which includes the annihilation of the positron.
- `linelist.out`: the list of the atomic lines. The first row is a header that starts with "#" and names the
  columns:
    - the vacuum wavelength in Angstroms;
    - the atomic number;
    - the ion stage;
    - the upper level index;
    - the lower level index.

  The level indices start at 0 in each ion. Each later row is one line, in the order of the line index that the
  emission and absorption codes of the packets use. Only the first job of a simulation writes this file.
- the reference tables of the setup, which rank 0 writes:
    - `timesteps.out`: each timestep, with its start, its midpoint, and its width in days;
    - `grid.out`: each {term}`propagation cell` that holds matter, with its model cell;
    - `modelgridrankassignments.out`: the first model cell and the cell count of each rank;
    - `nuclides.out`: each nuclide of the decay network, with its atomic number and its mass number;
    - `gammalinelist.out`: each gamma line of each nuclide, with its energy in MeV and its number of photons for
      each decay;
    - `bflist.out`: each bound-free continuum, with its ion, its lower level, and its upper level.
- `syn_dir.txt`: the unit vector of the theta = 0 direction of the direction bins. The code keeps the value `0 0 1`.
- `artis.pid`: the process id of rank 0 and the path of the model folder. `sn3d` and `exspec` stop when this file
  names a process that still runs in the same folder, and each program removes the file at the end.
- `gridsave_ts*.tmp` and `packets_*_ts*.tmp`: restart files that hold the state at the start of a timestep, after its
  grid update. A later job continues from that timestep and propagates its packets. A run with `VPKT_ON` also
  writes `vspecpol_nnnn_ts*.tmp`. The velocity grid map adds `vpkt_grid_nnnn_ts*.tmp`, and `VPKT_WRITE_CONTRIBS`
  adds `vpackets_nnnn_ts*.tmp`. A resumed run reads all of them.

`sn3d` also keeps an `output_0-0.txt` symlink in the model folder that points to the rank-0 log of the current job.

:::{warning}
A restart file must only be consistent with the binary that wrote it. A resumed run needs the same
`artisoptions.h` and the same source version as the run that wrote the files.
:::

## The format of the files

The files are text. The values of a row have a space between them. In this section, the times are in seconds, the
lengths in cm, the energies in erg, and the frequencies in Hz, unless a different unit is given.

### The light curves

`light_curve.out` has one row for each timestep that is complete, and no header:

| Column | Meaning |
| --- | --- |
| 1 | The midpoint of the timestep [days]. |
| 2 | The luminosity of the escaped r-packets in the rest frame [L_sun]. The time of a packet is its arrival time at the observer. |
| 3 | The luminosity of the escaped r-packets in the comoving frame [L_sun]. The time of a packet is its escape time multiplied by sqrt(1 - (`vmax`/c)^2), without the light travel time. |

`gamma_light_curve.out` has the same columns for the escaped {term}`gamma packets <gamma packet>`. Only a run with `KEEP_ESCAPED_GAMMAS`
writes it.

### The spectra

`spec.out` has one header row, and then one row for each of the 1000 frequency bins. The bins have a logarithmic
spacing from `NU_MIN_R` to `NU_MAX_R`.

- The header row holds a 0, and then the midpoint of each complete timestep [days].
- Each other row holds the middle frequency of a bin, and then the flux of each timestep.
- The flux is F_nu [erg/s/cm^2/Hz] for an observer at a distance of 1 Mpc. A packet goes into the timestep of its
  arrival time at the observer.

`gamma_spec.out` has the same layout for the escaped gamma packets, with the frequency range of the gamma rays. Only
a run with `KEEP_ESCAPED_GAMMAS` writes it. `specpol.out` has the timestep columns three times in sequence, for the
Stokes parameters I, Q, and U. Only a run with `POL_ON` writes `specpol.out`, `emissionpol.out`, and
`absorptionpol.out`.

### The emission and absorption files

`emission.out` and `emissiontrue.out` divide the flux of `spec.out` between the atomic processes of the last
emission. A packet whose last emission has no process is in `spec.out` and in neither of these two files, e.g.
after a {term}`thermalisation` with `RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY`. `absorption.out` gives the flux of the escaped packets for each process of their last absorption, in the
frequency bin of that absorption. Its columns therefore do not divide the flux of the same bin of `spec.out`. The
files have no header, no frequency column, and no time column.

- The files have one row for each frequency bin and each timestep of the model. The rows of all timesteps of the
  first frequency bin come first, then the rows of the second bin. The timesteps after the last complete timestep
  hold zeros. `timesteps.out` gives the number of timesteps.
- The unit of the values is the unit of `spec.out`.
- A column belongs to an ion. Its number, from 0, is the element index times the maximum ion count plus the ion
  index. The element index is the row of the element in `compositiondata.txt`, from 0. The ion index counts from the
  lowest ion stage of the element, from 0. The maximum ion count is the largest number of ion stages of one element.
- `emission.out` has a first group of ion columns for the bound-bound emission, a second group for the bound-free
  emission, and then one column for the free-free emission. The process is the last emission of the packet.
- `emissiontrue.out` has the same columns. The process is the last emission of the packet from the thermal energy
  of the gas, so a later scattering or {term}`macro atom` event does not change it.
- `absorption.out` has one group of ion columns. A packet goes into the frequency bin of its last absorption, and
  into the column of the ion of that line. Only a bound-bound absorption gives a value.
- `emissionpol.out` and `absorptionpol.out` hold, for each frequency bin, the rows of Stokes I for all timesteps,
  then the rows of Stokes Q, and then the rows of Stokes U.

### The direction-resolved files

The `speclc_angle_res` folder holds one file of each type for each of the 100 direction bins, e.g.
`spec_res_00.out` to `spec_res_99.out`. The types are:

- `light_curve_res`;
- `spec_res`;
- `emission_res`;
- `emissiontrue_res`;
- `absorption_res`.

A run with `POL_ON` adds `specpol_res`, `emissionpol_res`, and `absorptionpol_res`. Each file has the layout of the
file for all directions, e.g. `spec_res_07.out` has the layout of `spec.out`. The values are for an observer in that
direction, as if the emission were the same in all directions.

The number of a direction bin is 10 times its cos(theta) bin plus its phi bin. The polar angle theta is relative to
the z axis, and the cos(theta) bin 0 starts at cos(theta) = -1. The phi bins do not follow the azimuthal angle phi in
sequence. The bins 0 to 4 go from 2 pi down to pi, and the bins 5 to 9 go from 0 up to pi.

`scripts/mergeangleres.py` makes `light_curve_res.out`, `spec_res.out`, and `specpol_res.out`. `sn3d` and `exspec`
write only the files of the direction bins. A merged file holds the content of the files of all direction bins, one
after the other in the sequence of the bin numbers. Only a run with `POL_ON` has `specpol_res.out`.

### deposition.out

`deposition.out` has a header row that starts with `#` and names the columns. Then it has one row for each complete
timestep. A column with `Lsun` in its name is a luminosity [L_sun], and a column with `erg/s/g` is a power for each
gram of the model.

| Column | Meaning |
| --- | --- |
| `ts` | The timestep number. |
| `tmid_days`, `tmid_s` | The midpoint of the timestep in days and in seconds. |
| `total_dep_Lsun` | The deposition of the gamma rays, the positrons, the electrons, the alpha particles, and the fission fragments together. |
| `gammadep_Lsun` | The deposition of the gamma rays, from the paths of the packets. |
| `gammadep_discrete_Lsun` | The deposition of the gamma rays, from the absorption events of the packets. |
| `positrondep_Lsun`, `elecdep_Lsun`, `alphadep_Lsun` | The deposition of the positrons, the electrons, and the alpha particles, from the paths of the packets. |
| `positrondep_discrete_Lsun`, `elecdep_discrete_Lsun`, `alphadep_discrete_Lsun` | The deposition of the same particles, from the absorption events of the packets. |
| `spfission_dep_discrete_Lsun` | The deposition of the fission fragments. Only a model with spontaneous fission has this column. |
| `eps_gamma_Lsun`, `eps_elec_Lsun`, `eps_alpha_Lsun` | The emission of gamma-ray energy, and of the kinetic energy of the electrons and of the alpha particles, from the decays of the pellets. |
| `eps_positron_ana_Lsun`, `eps_elec_ana_Lsun`, `eps_alpha_ana_Lsun` | The emission of the kinetic energy of the positrons, the electrons, and the alpha particles, from the decay data and the abundances at the midpoint. |
| `eps_spfission_ana_Lsun` | The emission of the energy of the fission fragments, from the decay data and the abundances. Only a model with spontaneous fission has this column. |
| `eps_erg/s/g` | The emission of the pellets in gamma rays and in all particles. |
| `Qdot_betaminus_ana_erg/s/g`, `Qdotalpha_ana_erg/s/g`, `Qdotspfission_ana_erg/s/g` | The energy generation rate of one decay type, with the energy of the neutrinos, from the decay data and the abundances. Only a model with spontaneous fission has the fission column. |
| `Qdot_ana_erg/s/g` | The energy generation rate of all decay types, with the energy of the neutrinos. |

### The estimator files

An estimator file is free text. It has one block for each timestep and each model cell. The file of one rank,
`estimators_nnnn.out`, holds only the cells that the rank updates, and the combined file holds the cells of all
ranks. The first line of a block has pairs of a name and a value:

```text
timestep 0 modelgridindex 1 titeration 0 TR 17366.9 Te 17366.9 W 1 TJ 17366.9 grey_depth 667.454 thick 1 nne 2.05621e+08 Ye 0.434379 tdays    0.43
```

| Name | Meaning |
| --- | --- |
| `timestep`, `modelgridindex` | The timestep and the model cell, both from 0. |
| `titeration` | Always 0. |
| `TR`, `W` | The radiation temperature [K] and the dilution factor of the radiation field model. |
| `Te` | The electron temperature [K]. |
| `TJ` | The temperature [K] of a blackbody with the mean intensity of the cell. |
| `grey_depth` | The grey optical depth from the cell to the outer boundary. |
| `thick` | 0 for a cell with the detailed opacities, 1 for a cell with the grey opacity treatment, and 2 for a cell that is thick only for the virtual packets. |
| `nne` | The number density of the free electrons [cm^-3]. |
| `Ye` | The electron fraction. |
| `tdays` | The midpoint of the timestep [days]. |

An empty cell has only the line `timestep <n> modelgridindex <n> EMPTYCELL`.

The other lines of a block start with a label. A line for an element has `Z=<atomic number>` and then pairs of an
ion stage and a value, e.g. `2: 4.685e+05`. A block has these lines only for the elements with a mass fraction
above zero in the cell:

| Label | Meaning |
| --- | --- |
| `populations` | The number density of each ion stage [cm^-3]. `SUM` is the number density of the element. Then each nuclide of the element has its number density, and `<element>_otherstable` is the density of the stable nuclides that `model.txt` does not give. |
| `gamma_R` | The photoionisation rate coefficient of each ion stage [1/s]. This line, `gamma_R_integral`, and `gamma_NT` have no value for the highest ion stage of the element. |
| `gamma_R_integral` | The same coefficient from the radiation field model, without the bound-free estimators. Only a run with `DETAILED_BF_ESTIMATORS_FIRST_TIMESTEP` has this line. |
| `gamma_NT` | The non-thermal ionisation rate coefficient [1/s]. Only a run with a non-thermal solver has this line. |
| `corrphotoionrenorm` | The renormalisation factor of the tabulated photoionisation coefficient of each ion stage. Only a run with `USE_LUT_PHOTOION` has this line. |
| `gammaestimator` | The ionisation rate of each ion stage, divided by the population of its ground level [1/s]. Only a run with `USE_LUT_PHOTOION` has this line. This line and `corrphotoionrenorm` hold only the ion stages with a photoionisation from the ground level, and they have no value for the highest ion stage. |
| `heating_bf`, `heating_coll`, `heating_ff`, `cooling_ff`, `cooling_fb`, `cooling_coll` | The heating and the cooling rate of each ion [erg/s/cm^3]. Only a run with `WRITE_ION_HEATING_COOLING_RATES` has these lines. A cell with no thermal balance in the grid update has none of them. This is a cell in an LTE timestep, a cell with the grey opacity, and each cell in the first timestep of a job. `heating_bf` and `cooling_fb` have no value for the highest ion stage. `heating_coll` also needs `COLLISIONAL_HEATING_FROM_LEVELPOPS`. |
| `emission_ana:` | The emission rate of the gamma rays and of each particle type [erg/s/cm^3], from the decay data and the abundances. |
| `deposition:` | The deposition rate of the gamma rays and of each particle type [erg/s/cm^3]. |
| `heating:` | The free-free, the bound-free, and the collisional heating rate, and the heating rate of the deposition [erg/s/cm^3]. `heating_dep/total_dep` is the fraction of the deposition that heats the gas. |
| `cooling:` | The free-free, the free-bound, the collisional, and the adiabatic cooling rate [erg/s/cm^3]. |

### The nlte files

A run with no {term}`NLTE` levels writes no nlte file. An nlte file has a header row and then one row for each level of the
NLTE solution:

| Column | Meaning |
| --- | --- |
| `timestep`, `modelgridindex` | The timestep and the model cell, both from 0. |
| `Z`, `ionstage` | The atomic number and the ion stage. |
| `level` | The level index, from 0. The value -1 is the superlevel of the ion. |
| `n_LTE` | The number density of the level from the Boltzmann equation [cm^-3]. |
| `n_NLTE` | The number density of the level from the NLTE solution [cm^-3]. |
| `ion_popfrac` | `n_NLTE` divided by the number density of the ion. |

The {term}`LTE` timesteps write no rows. In the other timesteps, an element without an NLTE solution in a cell has no rows
for that cell:

- an element with no NLTE levels;
- an element with no mass in the cell;
- an element that fell back to LTE;
- each element of a cell with the {term}`grey opacity`.

### The radfield files

Only a run with `MULTIBIN_RADFIELD_FIRST_TIMESTEP` writes a radfield file. It has a header row and then one row for
each bin of the radiation field model. A cell in an LTE timestep, and a cell with the grey opacity, have no rows,
because they get no fit. The first timestep of each job also has no rows, because the grid update of that timestep
makes no fit.

| Column | Meaning |
| --- | --- |
| `timestep`, `modelgridindex` | The timestep and the model cell, both from 0. |
| `bin_num` | The bin index, from 0. The value -1 is the fit to the full spectrum. A value below -1 is a line with a detailed estimator. |
| `nu_lower`, `nu_upper` | The frequency limits of the bin. |
| `nuJ` | The integral of the frequency times the mean intensity over the bin [erg/s^2/cm^2/sr]. |
| `J` | The integral of the mean intensity over the bin [erg/s/cm^2/sr]. |
| `J_nu_avg` | `J` divided by the width of the bin [erg/s/cm^2/sr/Hz]. |
| `T_R`, `W` | The temperature [K] and the dilution factor of the diluted blackbody of the bin. |

The row of the full spectrum (`bin_num` -1) holds 0 in `nu_lower`, in `nu_upper`, and in `J_nu_avg`. The row of a
line with a detailed estimator holds the frequency of the line in `nu_lower` and in `nu_upper`, and the estimator of
the line in `J_nu_avg`. Its other values are -1.

### The macroatom file

`macroatom_allranks.out` is a file for a developer, which only a build with the constant `LOG_MACROATOM` of
`macroatom.cc` writes. It has a header row and then one row for each radiative de-excitation that ends a macro
atom. The table gives the columns in the sequence of the file:

| Column | Meaning |
| --- | --- |
| `timestep`, `modelgridindex` | The timestep and the model cell, both from 0. |
| `Z` | The atomic number. |
| `ionstage_in`, `ionstage_out` | The ion stage of the macro atom at its activation, and at the de-excitation. |
| `level_in`, `level_out` | The level index of the macro atom at its activation, and at the de-excitation. |
| `activline` | The line index of the absorption that activated the macro atom, or a negative value for another activation. |
| `nu_cmf_in`, `nu_cmf_out` | The frequency of the packet in the comoving frame [Hz], at the activation and after the de-excitation. |
| `nu_rf_in`, `nu_rf_out` | The frequency of the packet in the rest frame [Hz], at the activation and after the de-excitation. |

### The packet files

A packet file has a header row that starts with `#` and names the columns. Then it has one row for each packet of
the rank. A run without `KEEP_ESCAPED_GAMMAS` writes no rows for the escaped gamma packets.

| Column | Meaning |
| --- | --- |
| `number` | The number of the packet, which is unique only inside one rank. |
| `where` | The index of the propagation cell of the packet. `grid.out` gives the model cell. An empty propagation cell has no row in `grid.out`. |
| `type_id` | The type of the packet (see the table below). |
| `posx`, `posy`, `posz` | The position. |
| `dirx`, `diry`, `dirz` | The direction, as a unit vector. |
| `tdecay` | The decay time of the pellet. A packet of the initial energy of the cells holds the start time of the model, and no nuclide decays then. |
| `e_cmf`, `e_rf` | The energy in the comoving frame and in the rest frame. |
| `nu_cmf`, `nu_rf` | The frequency in the comoving frame and in the rest frame. For a pellet and for a packet of particles, `nu_cmf` holds the energy divided by the Planck constant, and `nu_rf` has no meaning. |
| `escape_type_id` | The type of the packet when it escaped, or 0. |
| `escape_time` | The time when the packet left the grid, or -1. |
| `emissiontype` | The process of the last emission (see below). |
| `trueemissiontype` | The process of the last emission from the thermal energy of the gas. |
| `em_posx`, `em_posy`, `em_posz` | The position of the last emission. A scattering of an r-packet and the Compton scattering of a gamma packet also set the position and `em_time`. The position is `nan` for a pellet, and for a packet of particles from a pellet decay, which sets only `em_time`. |
| `absorption_type`, `absorption_freq` | The process of the last absorption (see below), and the rest-frame frequency of the last absorption of an r-packet. `absorption_freq` is 0 until an r-packet absorption sets it. |
| `nscatterings` | The number of electron scatterings since the last emission. An event in a cell with the grey opacity also adds one. |
| `em_time` | The time of the last emission, or of the last scattering. |
| `stokes_q`, `stokes_u` | The Stokes parameters Q/I and U/I. Only a run with `POL_ON` has these columns. |
| `originated_from_particlenotgamma` | 1 if the decay of the pellet made a particle, and 0 if it made a gamma ray. The initial energy of the cells also holds 0, and its `pellet_decaytype` is -1. |
| `trueem_posx`, `trueem_posy`, `trueem_posz`, `trueem_time` | The position and the time of the last emission from the thermal energy of the gas. A packet with no such emission holds `nan` in the position and -1 in the time. |
| `pellet_nucindex` | The nuclide of the pellet decay, as the index of `nuclides.out`, or -1 for the initial energy of the cells. |
| `pellet_decaytype` | The decay type of the pellet, as a value of the list below. The value -1 is the initial energy of the cells. |

The values of `pellet_decaytype`:

- 0: an alpha decay;
- 1: an electron capture;
- 2: a beta-plus decay;
- 3: a beta-minus decay;
- 5: a spontaneous fission.

A run with `SAMPLE_RPKT_EMISSION` adds the columns of the sampled emission:

- `nemissiontype_updates`;
- `sampled_emissiontype`;
- `sampled_absorption_type` and `sampled_absorption_freq`;
- `sampled_em_posx`, `sampled_em_posy`, and `sampled_em_posz`;
- `sampled_em_time`.

The values of `type_id` and of `escape_type_id`:

| Value | Type |
| --- | --- |
| 0 | no type |
| 10 | gamma packet |
| 11 | r-packet |
| 12 | k-packet |
| 20 | energy of non-thermal leptons, before the code divides it between heating, ionisation, and excitation |
| 21, 22, 23 | electrons of a beta-minus decay, positrons of a beta-plus decay, and alpha particles, before the deposition |
| 24 | energy of alpha particles or fission fragments after the deposition |
| 32 | escaped packet |
| 100 | pellet |
| 120 | energy of a pellet that decayed before the start time of the model, or the initial energy of the cells with `INITIAL_PACKETS_ON` |

The values of `emissiontype` and of `trueemissiontype`:

- A value of 0 or more is the index of a line. The index counts the rows of `linelist.out` after its header, from 0.
- A negative value above -9999000 is a bound-free continuum. The index of the continuum in `bflist.out` is
  -1 - the value.
- -9999999 is free-free emission.
- -9999000 means that no emission has set the value.

The values of `absorption_type`:

| Value | Process |
| --- | --- |
| 0 or more | a line, with the same index as an emission type. A packet that had no absorption also holds 0. |
| -1 | free-free absorption |
| -2 | bound-free absorption |
| -3, -4, -5 | Compton scattering, photoelectric absorption, and pair production of a gamma ray |
| -6 | a reserved value of older packet files, for a pellet with no gamma-ray spectrum. The current code does not set it. |
| -7 | a pellet that decayed before the start time of the model, or the initial energy of the cells |
| -10 | a pellet decay to a particle |
| -11 | a bound-bound absorption in a binned opacity, which has no line index |

### The virtual packet files

Each rank writes its own files, and each file holds only the virtual packets of that rank. The types of file have
different normalisations:

- A `vspecpol` file holds the share of its rank, so the sum over the files of all ranks gives the result.
- A `vpkt_grid` file has no division by the rank count, and the packets of each rank hold the energy of the whole
  model. The mean over the files of all ranks therefore gives the result.
- The energies of a `vpackets` file have no division by the rank count either. Divide their sum over all ranks by
  the rank count.

`vspecpol_nnnn.out` has one block for each observer direction and each spectrum of that direction. The blocks of the
spectra of the first direction come first.

- The first row of a block holds a 0, and then the middle of each time bin [days], three times in sequence.
- Each other row holds the middle frequency of a bin. Then it holds the flux of each time bin for Stokes I, then for
  Stokes Q, and then for Stokes U.
- The flux has the unit of `spec.out`.
- The time bins and the frequency bins have a logarithmic spacing. `VSPEC_TIMEBINS`, `VSPEC_TIMEMIN`,
  `VSPEC_TIMEMAX`, `VSPEC_NUBINS`, `VSPEC_NUMIN`, and `VSPEC_NUMAX` of `vpkt.h` set them.

`vpkt_grid_nnnn.out` has one block for each observer direction and each wavelength range of the map. A block has
one row for each cell of the velocity grid, which has `VGRID_NY` times `VGRID_NZ` cells. The block of the first
wavelength range of the first observer direction comes first. Each row holds the two velocity
coordinates of the cell [cm/s] in the plane at a right angle to the observer direction. Then it holds three sums
over the virtual packets, for Stokes I, Q, and U.

Each virtual packet adds its rest-frame energy [erg], multiplied by the probability per unit solid angle of the
direction to the observer and by exp(-tau). The values have no factor of 4 pi.

`vpackets_nnnn.out` has a header row that starts with `#`. Then it has one row for each event of a real packet
whose virtual packet escaped in one direction or more:

| Column | Meaning |
| --- | --- |
| `emissiontype`, `trueemissiontype`, `absorption_type`, `absorption_freq` | The values of the real packet, as in the packet files. |
| `dir<n>_t_arrive_d` | The arrival time [days] of the virtual packet of observer direction n, or -1 if it did not escape. |
| `dir<n>_nu_rf` | The rest-frame frequency of that virtual packet, or -1. |
| `dir<n>_e_rf_<m>` | The rest-frame energy of that virtual packet for spectrum m, multiplied by the probability per unit solid angle of the direction to the observer and by exp(-tau), or 0. |

### The reference tables

| File | Layout |
| --- | --- |
| `timesteps.out` | A header row, and then the timestep number, the start, the midpoint, and the width of each timestep [days]. |
| `grid.out` | The propagation cell index and the model cell index of each propagation cell that holds matter. |
| `modelgridrankassignments.out` | A header row, and then one row for each rank: the rank, its first model cell, its cell count, and its count of cells with matter. |
| `nuclides.out` | A header row, and then the nuclide index, the atomic number, and the mass number of each nuclide. |
| `gammalinelist.out` | A header row, and then one row for each gamma-ray line in the sequence of the energy. The columns are the line index, the nuclide index, the atomic number, and the mass number. Then come the index of the line in its nuclide, the energy [MeV], and the number of photons for each decay. |
| `bflist.out` | The number of bound-free continua, and then one row for each continuum. The first three columns are its index, the element index, and the ion index. Then come the lower level index and the level index of the upper ion. All of them start at 0. |

## A new simulation removes the files of the previous one

A new simulation first removes the output files, the restart files, the `job_from_ts*` folders, and the `*.slurm`
job scripts of the previous simulation. It removes the same files as `scripts/clean.sh`, except `slurm-*.out`,
`machine.file.*`, and `core.*`, which can belong to the current job.

## Plotting and analysis

[artistools](https://github.com/artis-mcrt/artistools) reads these files. Some commands for a first look:

```sh
artistools describe <model folder>
artistools check <model folder>
artistools plotlightcurves <model folder>
artistools plotspectra -t <days> <model folder>
artistools plotestimators <model folder>
```

[Analyse the results](analysis.md) gives the commands for each task.

`artistools filterpackets` removes each packet except the escaped {term}`r-packets <r-packet>` from the packet files, which makes the
files smaller. `artistools corehours` sums the core hours of the jobs from the `output_0-0.txt` log of each job.
