#include "packet.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <format>
#include <iostream>
#include <ranges>
#include <span>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "artisoptions.h"
#include "constants.h"
#include "decay.h"
#include "globals.h"
#include "grid.h"
#include "input.h"
#include "mpi_logging.h"
#include "random.h"
#include "sn3d.h"
#include "vectors.h"

namespace {

constexpr auto get_packets_text_header() -> std::string {
  std::string header;
  header =
      "#number where type_id"
      " posx posy posz"
      " dirx diry dirz"
      " tdecay e_cmf e_rf nu_cmf nu_rf"
      " escape_type_id escape_time emissiontype trueemissiontype"
      " em_posx em_posy em_posz"
      " absorption_type absorption_freq nscatterings em_time";
  if constexpr (POL_ON) {
    header += " stokes1 stokes2 stokes3";
  }
  header +=
      " originated_from_particlenotgamma"
      " trueem_posx trueem_posy trueem_posz trueem_time"
      " pellet_nucindex pellet_decaytype";
  return header;
}

// Place pellet n with energy e_cmf_per_packet in cell m
void place_pellet(const double e_cmf_per_packet, const std::span<const double> en_cumulative, const int pktnumber,
                  Packet& pkt, const std::span<const double> energy_per_mass_nonemptymgi_decaypath) {
  const auto etot_simtime = en_cumulative.back();
  const double targetval = rng_uniform() * etot_simtime;

  // First choose a position for the pellet. In the cell.
  // first i such that en_cumulative[i] > targetval
  const int cellindex = static_cast<int>(std::ranges::upper_bound(en_cumulative, targetval) - en_cumulative.begin());
  assert_always(cellindex < grid::ngrid);

  pkt = Packet{};
  pkt.where = cellindex;
  pkt.number = pktnumber;  // record the packets number for debugging
  pkt.prop_time = globals::tmin;
  pkt.originated_from_particlenotgamma = false;

  pkt.pos = grid::get_propcell_random_position_tmin(cellindex);

  // ensure that the random position was inside the cell we selected
  assert_always(grid::get_cellindex_from_pos(pkt.pos, pkt.prop_time) == cellindex);

  const auto nonemptymgi = grid::get_propcell_nonemptymgi(cellindex);

  decay::setup_radioactive_pellet(e_cmf_per_packet, nonemptymgi, pkt, energy_per_mass_nonemptymgi_decaypath);

  // initial e_rf is probably never needed (e_rf is set at pellet decay time), but we
  // might as well give it a correct value since this code is fast and runs only once

  // pellet packet is moving with the homologous flow, so dir is proportional to pos
  pkt.dir = vec_norm(pkt.pos);  // assign dir = pos / vec_len(pos)
  const double dopplerfactor = calculate_doppler_nucmf_on_nurf(pkt.pos, pkt.dir, pkt.prop_time);
  pkt.e_rf = pkt.e_cmf / dopplerfactor;

  pkt.trueemissiontype = EMTYPE_NOTSET;
}

}  // anonymous namespace

void packet_init(std::span<Packet> packets)
// Subroutine that initialises the packets if we start a new simulation.
{
  MPI_Barrier_allranks();
  printlnlog("UNIFORM_PELLET_ENERGIES is {}", (UNIFORM_PELLET_ENERGIES ? "true" : "false"));

  printlnlog("INITIAL_PACKETS_ON is {}", (INITIAL_PACKETS_ON ? "on" : "off"));

  // The total number of pellets that we want to start with is just
  // npkts. The total energy of the pellets is given by etot.
  const double etot_tmodel_tinf = decay::get_global_etot_tmodel_tinf();

  printlnlog("etot {:g} (t_model to t_inf)", etot_tmodel_tinf);

  printlnlog("e_cmf per packet (t_model to t_inf) {:g} erg", etot_tmodel_tinf / MPKTS);

  const auto energy_per_mass_nonemptymgi_decaypath = decay::get_energy_per_mass_nonemptymgi_decaypath();

  // Need to get a normalisation factor
  auto en_cumulative = std::vector<double>(grid::ngrid);

  double etot_simtime = 0.;
  for (int propcellindex = 0; propcellindex < grid::ngrid; propcellindex++) {
    const int mgi = grid::get_propcell_modelgridindex(propcellindex);
    // some grid cells are empty
    if (mgi >= 0) {
      const auto nonemptymgi = grid::get_nonemptymgi_of_mgi(mgi);
      const double decay_en_per_mass =
          decay::get_modelcell_simtime_endecay_per_mass(nonemptymgi, energy_per_mass_nonemptymgi_decaypath.span());
      const auto initial_en_per_mass =
          (INITIAL_PACKETS_ON && USE_MODEL_INITIAL_ENERGY) ? grid::get_initenergyq(mgi) : 0.;

      etot_simtime += grid::get_propcell_volume_tmin(propcellindex) * grid::get_rho_tmin(mgi) *
                      (decay_en_per_mass + initial_en_per_mass);
    }
    en_cumulative[propcellindex] = etot_simtime;
  }
  assert_always(etot_simtime > 0);

  constexpr auto strtimelow{INITIAL_PACKETS_ON ? "tmodel" : "tmin"};
  printlnlog("etot ({} to tmax) {} erg", strtimelow, etot_simtime);

  // So energy per pellet is:
  const double e_cmf_per_packet = etot_simtime / MPKTS;
  printlnlog("e_cmf per packet ({} to tmax) {} erg", strtimelow, e_cmf_per_packet);

  // Now place the pellets in the ejecta and decide at what time they will decay.

  printlnlog("Placing pellets...");
  const auto allpktindices = std::ranges::iota_view{0, static_cast<int>(std::ssize(packets))};
  std::ranges::for_each(allpktindices, [&, e_cmf_per_packet](const int n) {
    place_pellet(e_cmf_per_packet, en_cumulative, n, packets[n], energy_per_mass_nonemptymgi_decaypath.span());
  });

  double e_cmf_total =
      std::ranges::fold_left(packets, 0., [](const double sum, const Packet& p) { return sum + p.e_cmf; });
  const double e_ratio = etot_simtime / e_cmf_total;
  printlnlog("packet energy sum {:g} should be {:g} normalisation factor: {:g}", e_cmf_total, etot_simtime, e_ratio);
  assert_always(std::isfinite(e_cmf_total));
  e_cmf_total *= e_ratio;
  for (auto& pkt : packets) {
    pkt.e_cmf *= e_ratio;
    pkt.e_rf *= e_ratio;
  }
  printlnlog("total energy that will be freed during simulation time: {:g} erg", e_cmf_total);
}

// read packets*.out text format file
auto read_text_packets(const std::string& filename) -> std::vector<Packet> {
  printlnlog("Reading {}", filename);
  auto packets_file = fstream_required(filename, std::ios::in);

  std::istringstream ssline;
  std::string line;
  std::vector<Packet> packets;
  std::getline(packets_file, line);  // read header line to make sure it matches
  assert_always(line == get_packets_text_header());

  packets.reserve(MPKTS);
  while (get_noncommentline(packets_file, line)) {
    ssline.clear();
    ssline.str(line);

    packets.emplace_back();
    Packet& pkt = packets.back();

    int pkt_type_in = 0;
    ssline >> pkt.number >> pkt.where >> pkt_type_in;
    pkt.type = static_cast<enum packet_type>(pkt_type_in);

    ssline >> pkt.pos[0] >> pkt.pos[1] >> pkt.pos[2];
    ssline >> pkt.dir[0] >> pkt.dir[1] >> pkt.dir[2];
    ssline >> pkt.tdecay;
    ssline >> pkt.e_cmf >> pkt.e_rf >> pkt.nu_cmf >> pkt.nu_rf;

    int escape_type = 0;
    ssline >> escape_type >> pkt.escape_time;
    pkt.escape_type = static_cast<enum packet_type>(escape_type);

    ssline >> pkt.emissiontype >> pkt.trueemissiontype;
    ssline >> pkt.em_pos[0] >> pkt.em_pos[1] >> pkt.em_pos[2];
    ssline >> pkt.absorptiontype >> pkt.absorptionfreq >> pkt.nscatterings;
    ssline >> pkt.em_time;

    if constexpr (POL_ON) {
      ssline >> pkt.stokes[0] >> pkt.stokes[1] >> pkt.stokes[2];
    }

    int int_originated_from_particlenotgamma = 0;
    ssline >> int_originated_from_particlenotgamma;
    pkt.originated_from_particlenotgamma = (int_originated_from_particlenotgamma != 0);

    ssline >> pkt.trueem_pos[0] >> pkt.trueem_pos[1] >> pkt.trueem_pos[2];
    ssline >> pkt.trueem_time;
    ssline >> pkt.pellet_nucindex;
    ssline >> pkt.pellet_decaytype;
  }

  if (std::ssize(packets) < MPKTS) {
    printlnlog("  found {} out of a possible {} packets.", std::ssize(packets), MPKTS);
  }
  packets.shrink_to_fit();
  return packets;
}

void write_text_packets(const std::string& filename, const std::span<const Packet> packets) {
  auto packets_file = fstream_required(filename, std::ios::out | std::ios::trunc);
  packets_file << get_packets_text_header() << '\n';

  for (const auto& pkt : packets) {
    if (!KEEP_ESCAPED_GAMMAS && pkt.type == TYPE_ESCAPE && pkt.escape_type == TYPE_GAMMA) {
      continue;
    }
    packets_file << pkt.number << ' ' << pkt.where << ' ' << std::to_underlying(pkt.type) << ' ';
    packets_file << pkt.pos[0] << ' ' << pkt.pos[1] << ' ' << pkt.pos[2] << ' ';
    packets_file << pkt.dir[0] << ' ' << pkt.dir[1] << ' ' << pkt.dir[2] << ' ';
    packets_file << pkt.tdecay << ' ';
    packets_file << pkt.e_cmf << ' ' << pkt.e_rf << ' ' << pkt.nu_cmf << ' ' << pkt.nu_rf << ' ';
    packets_file << std::to_underlying(pkt.escape_type) << ' ' << pkt.escape_time << ' ';
    packets_file << pkt.emissiontype << ' ' << pkt.trueemissiontype << ' ';
    packets_file << pkt.em_pos[0] << ' ' << pkt.em_pos[1] << ' ' << pkt.em_pos[2];
    packets_file << ' ' << pkt.absorptiontype << ' ' << pkt.absorptionfreq << ' ' << pkt.nscatterings << ' '
                 << pkt.em_time;
    if constexpr (POL_ON) {
      packets_file << ' ' << pkt.stokes[0] << ' ' << pkt.stokes[1] << ' ' << pkt.stokes[2];
    }
    packets_file << ' ' << static_cast<int>(pkt.originated_from_particlenotgamma);
    packets_file << ' ' << pkt.trueem_pos[0] << ' ' << pkt.trueem_pos[1] << ' ' << pkt.trueem_pos[2];
    packets_file << ' ' << pkt.trueem_time << ' ' << pkt.pellet_nucindex << ' ' << pkt.pellet_decaytype;
    packets_file << '\n';
  }
}

void read_temp_packetsfile(const int timestep, const int my_rank, std::vector<Packet>& packets) {
  // read binary packets file
  const auto filename = std::format("packets_{:04d}_ts{:d}.tmp", my_rank, timestep);

  printlnlog("Reading {}", filename);
  auto packets_file = fopen_required_uniqueptr(filename, "rb");
  std::int64_t packet_count_in_file = 0;
  assert_always(std::fread(&packet_count_in_file, sizeof(std::int64_t), 1, packets_file.get()) == 1);
  assert_always(packet_count_in_file > 0);
  assert_always(packet_count_in_file <= MPKTS);
  resize_exactly(packets, packet_count_in_file);
  assert_always(std::fread(packets.data(), sizeof(Packet), packet_count_in_file, packets_file.get()) ==
                static_cast<size_t>(packet_count_in_file));
  printlnlog("done");
}

void write_temp_packetsfile(const int timestep, const int my_rank, const std::span<const Packet> packets) {
  // write packets binary file (and retry if the write fails)
  const auto filename = std::format("packets_{:04d}_ts{:d}.tmp", my_rank, timestep);

  bool write_success = false;
  while (!write_success) {
    printlog("timestep {}: writing {}...", timestep, filename);
    FILE* packets_file = fopen(filename.c_str(), "wb");
    if (packets_file == nullptr) {
      printlnlog("ERROR: Could not open file '{}' for mode 'wb'.", filename);
      write_success = false;
    } else {
      auto packet_count = static_cast<std::int64_t>(std::ssize(packets));
      // write number of packets as header
      write_success = (std::fwrite(&packet_count, sizeof(std::int64_t), 1, packets_file) == 1);
      write_success = write_success &&
                      (std::fwrite(packets.data(), sizeof(Packet), packets.size(), packets_file) == packets.size());
      if (!write_success) {
        printlnlog("fwrite() FAILED! will retry...");
      }

      fclose(packets_file);
    }

    if (write_success) {
      printlnlog("done");
    }
  }
}
