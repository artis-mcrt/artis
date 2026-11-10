# Code map

This page connects the physics to the source files. [Methods and references](methods.md) describes each method, and
[AGENTS.md](https://github.com/artis-mcrt/artis/blob/main/AGENTS.md) gives the conventions of the code. All the C++
source is in the repository root. A comment at the top of each `.cc` file describes that file.

## The programs

| File | Purpose |
| --- | --- |
| `sn3d.cc` | The main simulation and its timestep loop. |
| `exspec.cc` | Makes the spectra and the light curves again from the packet files. |
| `unittests.cc` | The unit tests. |

## The modules

| File | Physics or task |
| --- | --- |
| `input.cc` | Reads `input.txt` and the atomic data. |
| `grid.cc` | Reads the model, and holds the model grid and the propagation grid. |
| `update_grid.cc` | The grid update of each timestep. |
| `update_packets.cc` | Sorts the packets, and moves each packet to the end of the timestep. |
| `packet.cc` | Makes the pellets, and reads and writes the packet files. |
| `decay.cc` | The nuclides, the decay chains, and the decay energies. |
| `gammapkt.cc` | The transport and the absorption of the gamma packets. |
| `nonthermal.cc` | The Spencer-Fano solver and the deposition of the fast particles. |
| `rpkt.cc` | The opacities and the movement of the r-packets. |
| `kpkt.cc` | The cooling rates and the k-packets. |
| `macroatom.cc` | The transitions of a macro atom. |
| `vpkt.cc` | The virtual packets. |
| `radfield.cc` | The estimators and the model of the radiation field. |
| `ratecoeff.cc` | The tables of the photoionisation and of the recombination rate coefficients. |
| `ltepop.cc` | The partition functions and the ionisation balance without NLTE levels. |
| `nltepop.cc` | The rate matrix and the solution of the NLTE populations. |
| `thermalbalance.cc` | The heating rates and the electron temperature. |
| `chargetransfer.cc` | The charge transfer reactions. |
| `spectrum_lightcurve.cc` | The spectra and the light curves of the escaped packets. |
| `stats.cc` | The counters of the packet events. |
| `mpi_logging.cc` | The log files. |

## The sequence of one timestep

`do_timestep()` in `sn3d.cc` calls these functions:

1. `update_grid()` (`update_grid.cc`) solves the state of the cells of the {term}`rank`.
2. `mpi_communicate_grid_properties()` sends the cell state to the other ranks.
3. `save_grid_and_packets()` writes the {term}`restart files <restart file>`, from the second {term}`timestep` of a job.
4. `walltime_sufficient_for_timestep()` tests the remaining wall time.
5. `zero_estimators()` sets the {term}`estimators <estimator>` to zero.
6. `update_packets()` (`update_packets.cc`) moves the packets of the rank.
7. `mpi_reduce_estimators()` sums the estimators of all ranks.
8. `write_deposition_file()` and `write_light_curves_and_spectra()` write the results of the timestep.

For a cell, `update_grid_cell()` normalises the estimators and calls `radfield::fit_parameters()`. In a timestep with
{term}`NLTE` populations, it then repeats these solutions in `update_grid.cc`:

- `nonthermal::solve_spencerfano()` for the non-thermal electrons;
- `call_T_e_finder()` for the electron temperature;
- `solve_nlte_pops_element()` for the populations of each element.

## The sequence of one packet

`do_packet()` in `update_packets.cc` calls one function for the type of the packet. The packet repeats this until it
reaches the end of the timestep or escapes.

| Packet type | Function | Next type |
| --- | --- | --- |
| `TYPE_RADIOACTIVE_PELLET` | `update_pellet()` | gamma, particle, deposited energy, or pre-k-packet |
| `TYPE_GAMMA` | `gammapkt::do_gamma()` | gamma, particle, deposited energy, or escape |
| `TYPE_NONTHERMAL_PREDEPOSIT_*` | `do_nonthermal_predeposit()` | deposited energy, the same type in the next timestep, or escape |
| `TYPE_NTLEPTON_DEPOSITED` | `nonthermal::do_ntlepton_deposit()` | k-packet or macro atom |
| `TYPE_NTALPHA_FISPROD_DEPOSITED` | `nonthermal::do_ntalpha_fisprod_deposit()` | k-packet |
| `TYPE_KPKT` | `kpkt::do_kpkt()` or `kpkt::do_kpkt_blackbody()` | r-packet or macro atom |
| `TYPE_PRE_KPKT` | `kpkt::do_kpkt_blackbody()` | r-packet |
| `TYPE_RPKT` | `do_rpkt()` | macro atom, k-packet, or escape |

A {term}`macro atom` has no step of its own. `do_macroatom()` in `macroatom.cc` runs directly after the event that activates
it, and it ends as an {term}`r-packet` or as a {term}`k-packet`.

## The three cell indices

`cellindex` selects a {term}`propagation cell`, `modelgridindex` selects a {term}`model grid cell`, and `nonemptymgi` selects a model
grid cell that contains matter. The comment above `get_propcell_modelgridindex()` in `grid.h` defines them.
