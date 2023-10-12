#include "packet.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "decay.h"
#include "grid.h"
#include "input.h"
#include "sn3d.h"
#include "vectors.h"

static void place_pellet(const double e0, const int cellindex, const int pktnumber, struct packet *pkt_ptr)
/// This subroutine places pellet n with energy e0 in cell m
{
  /// First choose a position for the pellet. In the cell.
  /// n is the index of the packet. m is the index for the grid cell.
  pkt_ptr->where = cellindex;
  pkt_ptr->number = pktnumber;  /// record the packets number for debugging
  pkt_ptr->prop_time = globals::tmin;
  // pkt_ptr->last_cross = BOUNDARY_NONE;
  pkt_ptr->originated_from_particlenotgamma = false;

  if constexpr (GRID_TYPE == GRID_SPHERICAL1D) {
    const double zrand = rng_uniform();
    const double r_inner = grid::get_cellcoordmin(cellindex, 0);
    const double r_outer = grid::get_cellcoordmax(cellindex, 0);
    // use equal volume probability distribution to select radius
    const double radius = pow(zrand * pow(r_inner, 3) + (1. - zrand) * pow(r_outer, 3), 1 / 3.);
    // assert_always(radius >= r_inner);
    // assert_always(radius <= r_outer);

    get_rand_isotropic_unitvec(pkt_ptr->pos);
    vec_scale(pkt_ptr->pos, radius);

  } else if constexpr (GRID_TYPE == GRID_CYLINDRICAL2D) {
    const double zrand1 = rng_uniform();
    const double rcyl_inner = grid::get_cellcoordmin(cellindex, 0);
    const double rcyl_outer = grid::get_cellcoordmax(cellindex, 0);
    // use equal area probability distribution to select radius
    const double rcyl_rand = sqrt(zrand1 * pow(rcyl_inner, 2) + (1. - zrand1) * pow(rcyl_outer, 2));
    const double theta_rand = rng_uniform() * 2 * PI;
    pkt_ptr->pos[0] = std::cos(theta_rand) * rcyl_rand;
    pkt_ptr->pos[1] = std::sin(theta_rand) * rcyl_rand;

    const double zrand2 = rng_uniform_pos();
    pkt_ptr->pos[2] = grid::get_cellcoordmin(cellindex, 1) + (zrand2 * grid::wid_init(cellindex, 1));

  } else if constexpr (GRID_TYPE == GRID_CARTESIAN3D) {
    for (int axis = 0; axis < 3; axis++) {
      const double zrand = rng_uniform_pos();
      pkt_ptr->pos[axis] = grid::get_cellcoordmin(cellindex, axis) + (zrand * grid::wid_init(cellindex, axis));
    }
  } else {
    assert_always(false);
  }

  // ensure that the random position was inside the cell we selected
  assert_always(grid::get_cellindex_from_pos(pkt_ptr->pos, pkt_ptr->prop_time) == cellindex);

  const int mgi = grid::get_cell_modelgridindex(cellindex);

  decay::setup_radioactive_pellet(e0, mgi, pkt_ptr);

  // initial e_rf is probably never needed (e_rf is set at pellet decay time), but we
  // might as well give it a correct value since this code is fast and runs only once

  // pellet packet is moving with the homologous flow, so dir is proportional to pos
  vec_norm(pkt_ptr->pos, pkt_ptr->dir);  // assign dir = pos / vec_len(pos)
  const double dopplerfactor = doppler_packet_nucmf_on_nurf(pkt_ptr->pos, pkt_ptr->dir, pkt_ptr->prop_time);
  pkt_ptr->e_rf = pkt_ptr->e_cmf / dopplerfactor;

  pkt_ptr->trueemissiontype = EMTYPE_NOTSET;
}

void packet_init(struct packet *pkt)
/// Subroutine that initialises the packets if we start a new simulation.
{
#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  printout("UNIFORM_PELLET_ENERGIES is %s\n", (UNIFORM_PELLET_ENERGIES ? "true" : "false"));

  printout("INITIAL_PACKETS_ON is %s\n", (INITIAL_PACKETS_ON ? "on" : "off"));

  /// The total number of pellets that we want to start with is just
  /// npkts. The total energy of the pellets is given by etot.
  const double etot_tinf = decay::get_global_etot_t0_tinf();

  printout("etot %g (t_0 to t_inf)\n", etot_tinf);

  const double e0_tinf = etot_tinf / globals::npkts;
  printout("packet e0 (t_0 to t_inf) %g erg\n", e0_tinf);

  decay::setup_decaypath_energy_per_mass();

  // Need to get a normalisation factor.
  auto en_cumulative = std::vector<double>(grid::ngrid);

  double norm = 0.0;
  for (int m = 0; m < grid::ngrid; m++) {
    const int mgi = grid::get_cell_modelgridindex(m);
    if (mgi < grid::get_npts_model())  // some grid cells are empty
    {
      double q = decay::get_modelcell_simtime_endecay_per_mass(mgi);
      if constexpr (INITIAL_PACKETS_ON && USE_MODEL_INITIAL_ENERGY) {
        q += grid::get_initenergyq(mgi);
      }

      norm += grid::get_gridcell_volume_tmin(m) * grid::get_rho_tmin(mgi) * q;
    }
    en_cumulative[m] = norm;
  }
  assert_always(norm > 0);

  const double etot = norm;
  /// So energy per pellet is
  const double e0 = etot / globals::npkts;
  printout("packet e0 (in time range) %g erg\n", e0);

  printout("etot %g erg (in time range) erg\n", etot);

  // Now place the pellets in the ejecta and decide at what time they will decay.

  if (globals::npkts > MPKTS) {
    printout("Too many packets. Abort.\n");
    abort();
  }

  printout("Placing pellets...\n");
  for (int n = 0; n < globals::npkts; n++) {
    const double zrand = rng_uniform();
    const double targetval = zrand * norm;

    // first en_cumulative[i] such that en_cumulative[i] > targetval
    auto upperval = std::upper_bound(en_cumulative.cbegin(), en_cumulative.cend(), targetval);
    assert_always(upperval != en_cumulative.end());
    const ptrdiff_t cellindex = std::distance(en_cumulative.cbegin(), upperval);

    place_pellet(e0, cellindex, n, &pkt[n]);
  }

  decay::free_decaypath_energy_per_mass();  // will no longer be needed after packets are set up

  double e_cmf_total = 0.;
  for (int n = 0; n < globals::npkts; n++) {
    pkt[n].interactions = 0;
    e_cmf_total += pkt[n].e_cmf;
  }
  const double e_ratio = etot / e_cmf_total;
  printout("packet energy sum %g should be %g normalisation factor: %g\n", e_cmf_total, etot, e_ratio);
  assert_always(std::isfinite(e_cmf_total));
  e_cmf_total *= e_ratio;
  for (int n = 0; n < globals::npkts; n++) {
    pkt[n].e_cmf *= e_ratio;
    pkt[n].e_rf *= e_ratio;
  }
  printout("total energy that will be freed during simulation time: %g erg\n", e_cmf_total);
}

void write_packets(char filename[], const struct packet *const pkt) {
  // write packets text file
  FILE *packets_file = fopen_required(filename, "w");
  fprintf(packets_file,
          "#number where type_id posx posy posz dirx diry dirz last_cross tdecay e_cmf e_rf nu_cmf nu_rf "
          "escape_type_id escape_time next_trans interactions last_event emissiontype trueemissiontype "
          "em_posx em_posy em_posz absorption_type absorption_freq nscatterings em_time absorptiondirx absorptiondiry "
          "absorptiondirz stokes1 stokes2 stokes3 pol_dirx pol_diry pol_dirz originated_from_positron "
          "true_emission_velocity trueem_time pellet_nucindex\n");
  for (int i = 0; i < globals::npkts; i++) {
    fprintf(packets_file, "%d ", pkt[i].number);
    fprintf(packets_file, "%d ", pkt[i].where);
    fprintf(packets_file, "%d ", pkt[i].type);
    fprintf(packets_file, "%lg %lg %lg ", pkt[i].pos[0], pkt[i].pos[1], pkt[i].pos[2]);
    fprintf(packets_file, "%lg %lg %lg ", pkt[i].dir[0], pkt[i].dir[1], pkt[i].dir[2]);
    fprintf(packets_file, "%d ", pkt[i].last_cross);
    fprintf(packets_file, "%g ", pkt[i].tdecay);
    fprintf(packets_file, "%g ", pkt[i].e_cmf);
    fprintf(packets_file, "%g ", pkt[i].e_rf);
    fprintf(packets_file, "%g ", pkt[i].nu_cmf);
    fprintf(packets_file, "%g ", pkt[i].nu_rf);
    fprintf(packets_file, "%d ", pkt[i].escape_type);
    fprintf(packets_file, "%g ", pkt[i].escape_time);
    fprintf(packets_file, "%d ", pkt[i].next_trans);
    fprintf(packets_file, "%d ", pkt[i].interactions);
    fprintf(packets_file, "%d ", pkt[i].last_event);
    fprintf(packets_file, "%d ", pkt[i].emissiontype);
    fprintf(packets_file, "%d ", pkt[i].trueemissiontype);
    fprintf(packets_file, "%lg %lg %lg ", pkt[i].em_pos[0], pkt[i].em_pos[1], pkt[i].em_pos[2]);
    fprintf(packets_file, "%d ", pkt[i].absorptiontype);
    fprintf(packets_file, "%lg ", pkt[i].absorptionfreq);
    fprintf(packets_file, "%d ", pkt[i].nscatterings);
    fprintf(packets_file, "%g ", pkt[i].em_time);
    fprintf(packets_file, "%lg %lg %lg ", pkt[i].absorptiondir[0], pkt[i].absorptiondir[1], pkt[i].absorptiondir[2]);
    fprintf(packets_file, "%lg %lg %lg ", pkt[i].stokes[0], pkt[i].stokes[1], pkt[i].stokes[2]);
    fprintf(packets_file, "%lg %lg %lg ", pkt[i].pol_dir[0], pkt[i].pol_dir[1], pkt[i].pol_dir[2]);
    fprintf(packets_file, "%d ", static_cast<int>(pkt[i].originated_from_particlenotgamma));
    fprintf(packets_file, "%g ", pkt[i].trueemissionvelocity);
    fprintf(packets_file, "%g ", pkt[i].trueem_time);
    fprintf(packets_file, "%d ", pkt[i].pellet_nucindex);
    fprintf(packets_file, "\n");
  }
  fclose(packets_file);
}

void read_temp_packetsfile(const int timestep, const int my_rank, struct packet *const pkt) {
  // read packets binary file
  char filename[MAXFILENAMELENGTH];
  snprintf(filename, MAXFILENAMELENGTH, "packets_%.4d_ts%d.tmp", my_rank, timestep);

  printout("Reading %s...", filename);
  FILE *packets_file = fopen_required(filename, "rb");
  assert_always(std::fread(pkt, sizeof(struct packet), globals::npkts, packets_file) == (size_t)globals::npkts);
  // read_packets(packets_file);
  fclose(packets_file);
  printout("done\n");
}

auto verify_temp_packetsfile(const int timestep, const int my_rank, const struct packet *const pkt) -> bool {
  // return true if verification is good, otherwise return false

  // read packets binary file
  char filename[MAXFILENAMELENGTH];
  snprintf(filename, MAXFILENAMELENGTH, "packets_%.4d_ts%d.tmp", my_rank, timestep);

  printout("Verifying file %s...", filename);
  FILE *packets_file = fopen_required(filename, "rb");
  struct packet pkt_in;
  bool readback_passed = true;
  for (int n = 0; n < globals::npkts; n++) {
    assert_always(std::fread(&pkt_in, sizeof(struct packet), 1, packets_file) == 1);
    if (pkt_in != pkt[n]) {
      printout("failed on packet %d\n", n);
      printout(" compare number %ld %ld\n", pkt_in.number, pkt[n].number);
      printout(" compare nu_cmf %lg %lg\n", pkt_in.nu_cmf, pkt[n].nu_cmf);
      printout(" compare e_rf %lg %lg\n", pkt_in.e_rf, pkt[n].e_rf);
      readback_passed = false;
    }
  }
  fclose(packets_file);
  if (readback_passed) {
    printout("  verification passed\n");
  } else {
    printout("  verification FAILED\n");
  }
  return readback_passed;
}

void read_packets(const char filename[], struct packet *pkt) {
  // read packets*.out text format file
  std::ifstream packets_file(filename);
  assert_always(packets_file.is_open());

  std::string line;

  int packets_read = 0;
  while (std::getline(packets_file, line)) {
    if (lineiscommentonly(line)) {
      continue;
    }

    packets_read++;
    const int i = packets_read - 1;

    if (i > globals::npkts - 1) {
      printout(
          "ERROR: More data found beyond packet %d (expecting %d packets). Recompile exspec with the correct number "
          "of packets. Run (wc -l < packets00_0000.out) to count them.\n",
          packets_read, globals::npkts);
      abort();
    }

    std::istringstream ssline(line);

    int pkt_type_in = 0;
    ssline >> pkt[i].number >> pkt[i].where >> pkt_type_in;
    pkt[i].type = static_cast<enum packet_type>(pkt_type_in);

    ssline >> pkt[i].pos[0] >> pkt[i].pos[1] >> pkt[i].pos[2];

    ssline >> pkt[i].dir[0] >> pkt[i].dir[1] >> pkt[i].dir[2];

    int last_cross_in = 0;
    ssline >> last_cross_in;
    pkt[i].last_cross = static_cast<enum cell_boundary>(last_cross_in);

    ssline >> pkt[i].tdecay;

    ssline >> pkt[i].e_cmf >> pkt[i].e_rf >> pkt[i].nu_cmf >> pkt[i].nu_rf;

    int escape_type = 0;
    ssline >> escape_type >> pkt[i].escape_time;
    pkt[i].escape_type = static_cast<enum packet_type>(escape_type);

    ssline >> pkt[i].next_trans >> pkt[i].interactions >> pkt[i].last_event;
    assert_always(pkt[i].interactions >= 0);

    ssline >> pkt[i].emissiontype >> pkt[i].trueemissiontype;

    ssline >> pkt[i].em_pos[0] >> pkt[i].em_pos[1] >> pkt[i].em_pos[2];

    ssline >> pkt[i].absorptiontype >> pkt[i].absorptionfreq >> pkt[i].nscatterings;

    ssline >> pkt[i].em_time;

    ssline >> pkt[i].absorptiondir[0] >> pkt[i].absorptiondir[1] >> pkt[i].absorptiondir[2];

    ssline >> pkt[i].stokes[0] >> pkt[i].stokes[1] >> pkt[i].stokes[2];

    ssline >> pkt[i].pol_dir[0] >> pkt[i].pol_dir[1] >> pkt[i].pol_dir[2];

    int int_originated_from_particlenotgamma = 0;
    ssline >> int_originated_from_particlenotgamma;
    pkt[i].originated_from_particlenotgamma = (int_originated_from_particlenotgamma != 0);

    ssline >> pkt[i].trueemissionvelocity;

    ssline >> pkt[i].trueem_time;

    ssline >> pkt[i].pellet_nucindex;
  }

  if (packets_read < globals::npkts) {
    printout(
        "ERROR: Read failed after packet %d (expecting %d packets). Recompile exspec with the correct number of "
        "packets. Run (wc -l < packets00_0000.out) to count them.\n",
        packets_read, globals::npkts);
    abort();
  }

  packets_file.close();
}
