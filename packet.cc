#include "packet.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <format>
#include <fstream>
#include <iostream>
#include <ranges>
#include <span>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#pragma clang unsafe_buffer_usage begin
#include <mpi.h>
#pragma clang unsafe_buffer_usage end

#include "artisoptions.h"
#include "constants.h"
#include "decay.h"
#include "globals.h"
#include "grid.h"
#include "input.h"
#include "random.h"
#include "sn3d.h"
#include "vectors.h"

namespace {

// Place pellet n with energy e_cmf_per_packet in cell m
void place_pellet(const double e_cmf_per_packet, const int cellindex, const int pktnumber, Packet& pkt) {
  // First choose a position for the pellet. In the cell.
  pkt.where = cellindex;
  pkt.number = pktnumber;  // record the packets number for debugging
  pkt.prop_time = globals::tmin;
  pkt.originated_from_particlenotgamma = false;

  if constexpr (GRID_TYPE == GridType::SPHERICAL1D) {
    const double zrand = rng_uniform();
    const double r_inner = grid::get_cellcoordmin(cellindex, 0);
    const double r_outer = grid::get_cellcoordmax(cellindex, 0);
    // use equal volume probability distribution to select radius
    const double radius = pow((zrand * pow(r_inner, 3)) + ((1. - zrand) * pow(r_outer, 3)), 1 / 3.);
    // assert_always(radius >= r_inner);
    // assert_always(radius <= r_outer);

    pkt.pos = vec_scale(get_rand_isotropic_unitvec(), radius);

  } else if constexpr (GRID_TYPE == GridType::CYLINDRICAL2D) {
    const double zrand = rng_uniform_pos();
    const double rcyl_inner = grid::get_cellcoordmin(cellindex, 0);
    const double rcyl_outer = grid::get_cellcoordmax(cellindex, 0);
    // use equal area probability distribution to select radius
    const double rcyl_rand = std::sqrt((zrand * std::pow(rcyl_inner, 2)) + ((1. - zrand) * std::pow(rcyl_outer, 2)));
    const double theta_rand = rng_uniform() * 2 * PI;
    pkt.pos = {std::cos(theta_rand) * rcyl_rand, std::sin(theta_rand) * rcyl_rand,
               grid::get_cellcoordmin(cellindex, 1) + (rng_uniform_pos() * grid::wid_init(cellindex, 1))};

  } else if constexpr (GRID_TYPE == GridType::CARTESIAN3D) {
    for (int axis = 0; axis < 3; axis++) {
      pkt.pos[axis] = grid::get_cellcoordmin(cellindex, axis) + (rng_uniform_pos() * grid::wid_init(cellindex, axis));
    }
  } else {
    assert_always(false);
  }

  // ensure that the random position was inside the cell we selected
  assert_always(grid::get_cellindex_from_pos(pkt.pos, pkt.prop_time) == cellindex);

  const auto nonemptymgi = grid::get_propcell_nonemptymgi(cellindex);

  decay::setup_radioactive_pellet(e_cmf_per_packet, nonemptymgi, pkt);

  // initial e_rf is probably never needed (e_rf is set at pellet decay time), but we
  // might as well give it a correct value since this code is fast and runs only once

  // pellet packet is moving with the homologous flow, so dir is proportional to pos
  pkt.dir = vec_norm(pkt.pos);  // assign dir = pos / vec_len(pos)
  const double dopplerfactor = calculate_doppler_nucmf_on_nurf(pkt.pos, pkt.dir, pkt.prop_time);
  pkt.e_rf = pkt.e_cmf / dopplerfactor;

  pkt.trueemissiontype = EMTYPE_NOTSET;
}

}  // anonymous namespace

void packet_init(std::span<Packet> pkt)
// Subroutine that initialises the packets if we start a new simulation.
{
  MPI_Barrier(MPI_COMM_WORLD);
  printout("UNIFORM_PELLET_ENERGIES is %s\n", (UNIFORM_PELLET_ENERGIES ? "true" : "false"));

  printout("INITIAL_PACKETS_ON is %s\n", (INITIAL_PACKETS_ON ? "on" : "off"));

  // The total number of pellets that we want to start with is just
  // npkts. The total energy of the pellets is given by etot.
  const double etot_tmodel_tinf = decay::get_global_etot_tmodel_tinf();

  printout("etot %g (t_model to t_inf)\n", etot_tmodel_tinf);

  printout("e_cmf per packet (t_model to t_inf) %g erg\n", etot_tmodel_tinf / globals::npkts);

  decay::setup_decaypath_energy_per_mass();

  // Need to get a normalisation factor
  auto en_cumulative = std::vector<double>(grid::ngrid);

  double etot_simtime = 0.;
  for (int m = 0; m < grid::ngrid; m++) {
    const int mgi = grid::get_propcell_modelgridindex(m);
    if (mgi < grid::get_npts_model() && grid::get_numpropcells(mgi) > 0)  // some grid cells are empty
    {
      const auto nonemptymgi = grid::get_nonemptymgi_of_mgi(mgi);
      double q = decay::get_modelcell_simtime_endecay_per_mass(nonemptymgi);
      if constexpr (INITIAL_PACKETS_ON && USE_MODEL_INITIAL_ENERGY) {
        q += grid::get_initenergyq(mgi);
      }

      etot_simtime += grid::get_propcell_volume_tmin(m) * grid::get_rho_tmin(mgi) * q;
    }
    en_cumulative[m] = etot_simtime;
  }
  assert_always(etot_simtime > 0);

  constexpr auto strtimelow{INITIAL_PACKETS_ON ? "tmodel" : "tmin"};
  logprintlnfmt("etot ({} to tmax) {} erg", strtimelow, etot_simtime);

  // So energy per pellet is:
  const double e_cmf_per_packet = etot_simtime / globals::npkts;
  logprintlnfmt("e_cmf per packet ({} to tmax) {} erg", strtimelow, e_cmf_per_packet);

  // Now place the pellets in the ejecta and decide at what time they will decay.

  logprintlnfmt("Placing pellets...");
  const auto allpkts = std::ranges::iota_view{0, globals::npkts};
  std::for_each(allpkts.begin(), allpkts.end(), [&, etot_simtime, e_cmf_per_packet](const int n) {
    pkt[n] = Packet{};
    const double targetval = rng_uniform() * etot_simtime;

    // first i such that en_cumulative[i] > targetval
    const int cellindex = static_cast<int>(std::ranges::upper_bound(en_cumulative, targetval) - en_cumulative.cbegin());
    assert_always(cellindex < grid::ngrid);

    place_pellet(e_cmf_per_packet, cellindex, n, pkt[n]);
  });

  decay::free_decaypath_energy_per_mass();  // will no longer be needed after packets are set up

  double e_cmf_total = std::ranges::fold_left(pkt, 0., [](const double sum, const Packet& p) { return sum + p.e_cmf; });
  const double e_ratio = etot_simtime / e_cmf_total;
  printout("packet energy sum %g should be %g normalisation factor: %g\n", e_cmf_total, etot_simtime, e_ratio);
  assert_always(std::isfinite(e_cmf_total));
  e_cmf_total *= e_ratio;
  for (int n = 0; n < globals::npkts; n++) {
    pkt[n].e_cmf *= e_ratio;
    pkt[n].e_rf *= e_ratio;
  }
  printout("total energy that will be freed during simulation time: %g erg\n", e_cmf_total);
}

// write packets text file
void write_packets(const std::string& filename, std::span<const Packet> pkt) {
  auto packets_file = std::fstream(filename, std::ios::out | std::ios::trunc);
  assert_always(packets_file.is_open());
  packets_file << "#number where type_id posx posy posz dirx diry dirz tdecay e_cmf e_rf nu_cmf nu_rf "
                  "escape_type_id escape_time emissiontype trueemissiontype "
                  "em_posx em_posy em_posz absorption_type absorption_freq nscatterings em_time stokes1 stokes2 "
                  "stokes3 originated_from_particlenotgamma "
                  "trueem_posx trueem_posy trueem_posz trueem_time pellet_nucindex pellet_decaytype\n";

  for (int i = 0; i < globals::npkts; i++) {
    packets_file << pkt[i].number << ' ' << pkt[i].where << ' ' << std::to_underlying(pkt[i].type) << ' ';
    packets_file << pkt[i].pos[0] << ' ' << pkt[i].pos[1] << ' ' << pkt[i].pos[2] << ' ';
    packets_file << pkt[i].dir[0] << ' ' << pkt[i].dir[1] << ' ' << pkt[i].dir[2] << ' ';
    packets_file << pkt[i].tdecay << ' ';
    packets_file << pkt[i].e_cmf << ' ' << pkt[i].e_rf << ' ' << pkt[i].nu_cmf << ' ' << pkt[i].nu_rf << ' ';
    packets_file << std::to_underlying(pkt[i].escape_type) << ' ' << pkt[i].escape_time << ' ';
    packets_file << pkt[i].emissiontype << ' ' << pkt[i].trueemissiontype << ' ';
    packets_file << pkt[i].em_pos[0] << ' ' << pkt[i].em_pos[1] << ' ' << pkt[i].em_pos[2] << ' '
                 << pkt[i].absorptiontype << ' ' << pkt[i].absorptionfreq << ' ' << pkt[i].nscatterings << ' '
                 << pkt[i].em_time << ' ';
    packets_file << pkt[i].stokes[0] << ' ' << pkt[i].stokes[1] << ' ' << pkt[i].stokes[2] << ' ';
    packets_file << static_cast<int>(pkt[i].originated_from_particlenotgamma) << ' ' << pkt[i].trueem_pos[0] << ' '
                 << pkt[i].trueem_pos[1] << ' ' << pkt[i].trueem_pos[2] << ' ' << pkt[i].trueem_time << ' '
                 << pkt[i].pellet_nucindex << ' ' << pkt[i].pellet_decaytype;
    packets_file << '\n';
  }
}

void read_temp_packetsfile(const int timestep, const int my_rank, std::span<Packet> pkt) {
  // read binary packets file
  const auto filename = std::format("packets_{:04d}_ts{:d}.tmp", my_rank, timestep);

  logprintlnfmt("Reading {}", filename);
  auto packets_file = fopen_required_uniqueptr(filename, "rb");
  assert_always(std::fread(pkt.data(), sizeof(Packet), globals::npkts, packets_file.get()) ==
                static_cast<size_t>(globals::npkts));

  printout("done\n");
}

auto verify_temp_packetsfile(const int timestep, const int my_rank, std::span<const Packet> pkt) -> bool {
  // return true if verification is good, otherwise return false

  // read binary packets file
  const auto filename = std::format("packets_{:04d}_ts{:d}.tmp", my_rank, timestep);

  logprintlnfmt("Verifying file {}", filename);
  auto packets_file = fopen_required_uniqueptr(filename, "rb");
  Packet pkt_in;
  bool readback_passed = true;
  for (int n = 0; n < globals::npkts; n++) {
    assert_always(std::fread(&pkt_in, sizeof(Packet), 1, packets_file.get()) == 1);
    if (pkt_in != pkt[n]) {
      printout("failed on packet %d\n", n);
      printout(" compare number %d %d\n", pkt_in.number, pkt[n].number);
      printout(" compare nu_cmf %lg %lg\n", pkt_in.nu_cmf, pkt[n].nu_cmf);
      printout(" compare e_rf %lg %lg\n", pkt_in.e_rf, pkt[n].e_rf);
      readback_passed = false;
    }
  }
  if (readback_passed) {
    printout("  verification passed\n");
  } else {
    printout("  verification FAILED\n");
  }
  return readback_passed;
}

auto read_packets(const std::string& filename, std::span<Packet> packets) -> std::span<Packet> {
  // read packets*.out text format file
  auto packets_file = fstream_required(filename, std::ios::in);

  std::string line;

  int packets_read = 0;
  while (get_noncommentline(packets_file, line)) {
    packets_read++;
    const int i = packets_read - 1;

    if (i > globals::npkts - 1) {
      printout(
          "ERROR: More data found beyond packet %d (expecting %d packets). Recompile exspec with the correct number "
          "of packets. Run (wc -l < packets00_0000.out) to count them.\n",
          packets_read, globals::npkts);
      std::abort();
    }

    std::istringstream ssline(line);

    int pkt_type_in = 0;
    ssline >> packets[i].number >> packets[i].where >> pkt_type_in;
    packets[i].type = static_cast<enum packet_type>(pkt_type_in);

    ssline >> packets[i].pos[0] >> packets[i].pos[1] >> packets[i].pos[2];

    ssline >> packets[i].dir[0] >> packets[i].dir[1] >> packets[i].dir[2];

    ssline >> packets[i].tdecay;

    ssline >> packets[i].e_cmf >> packets[i].e_rf >> packets[i].nu_cmf >> packets[i].nu_rf;

    int escape_type = 0;
    ssline >> escape_type >> packets[i].escape_time;
    packets[i].escape_type = static_cast<enum packet_type>(escape_type);

    ssline >> packets[i].emissiontype >> packets[i].trueemissiontype;

    ssline >> packets[i].em_pos[0] >> packets[i].em_pos[1] >> packets[i].em_pos[2];

    ssline >> packets[i].absorptiontype >> packets[i].absorptionfreq >> packets[i].nscatterings;

    ssline >> packets[i].em_time;

    ssline >> packets[i].stokes[0] >> packets[i].stokes[1] >> packets[i].stokes[2];

    int int_originated_from_particlenotgamma = 0;
    ssline >> int_originated_from_particlenotgamma;
    packets[i].originated_from_particlenotgamma = (int_originated_from_particlenotgamma != 0);

    ssline >> packets[i].trueem_pos[0] >> packets[i].trueem_pos[1] >> packets[i].trueem_pos[2];

    ssline >> packets[i].trueem_time;

    ssline >> packets[i].pellet_nucindex;
  }

  if (packets_read < globals::npkts) {
    printout(
        "ERROR: Read failed after packet %d (expecting %d packets). Recompile exspec with the correct number of "
        "packets. Run (wc -l < packets00_0000.out) to count them.\n",
        packets_read, globals::npkts);
    std::abort();
  }
  return packets;
}
