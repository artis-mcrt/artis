#include "sn3d.h"
#include "grid_init.h"
#include "decay.h"
#include "packet_init.h"
#include "vectors.h"


static void place_pellet(const double e0, const int cellindex, const int pktnumber, PKT *pkt_ptr)
/// This subroutine places pellet n with energy e0 in cell m
{
  /// First choose a position for the pellet. In the cell.
  /// n is the index of the packet. m is the index for the grid cell.
  pkt_ptr->where = cellindex;
  pkt_ptr->number = pktnumber;  ///record the packets number for debugging
  pkt_ptr->originated_from_positron = false;

  if (grid_type == GRID_SPHERICAL1D)
  {
    const double zrand3 = gsl_rng_uniform(rng);
    const double r_inner = get_cellcoordmin(cellindex, 0);
    const double r_outer = get_cellcoordmin(cellindex, 0) + wid_init(cellindex);
    const double radius = pow(zrand3 * pow(r_inner, 3) + (1. - zrand3) * pow(r_outer, 3), 1/3.);
    // assert(radius >= r_inner);
    // assert(radius <= r_outer);

    get_rand_isotropic_unitvec(pkt_ptr->pos);
    vec_scale(pkt_ptr->pos, radius);
  }
  else
  {
    for (int axis = 0; axis < 3; axis++)
    {
      const double zrand = gsl_rng_uniform_pos(rng);
      pkt_ptr->pos[axis] = get_cellcoordmin(cellindex, axis) + (zrand * wid_init(0));
    }
  }

  const int mgi = cell[cellindex].modelgridindex;

  double cumulative_decay_energy_per_mass[DECAYPATH_COUNT];
  for (int i = 0; i < DECAYPATH_COUNT; i++)
  {
    const double lower_sum = ((i > 0) ? cumulative_decay_energy_per_mass[i - 1] : 0);
    cumulative_decay_energy_per_mass[i] = lower_sum + get_simtime_endecay_per_ejectamass(mgi, i);
  }

  const double zrand_chain = gsl_rng_uniform(rng) * cumulative_decay_energy_per_mass[DECAYPATH_COUNT - 1];
  enum decaypathways decaypath = DECAYPATH_COUNT;
  for (int i = 0; i < DECAYPATH_COUNT; i++)
  {
    if (zrand_chain <= cumulative_decay_energy_per_mass[i])
    {
      decaypath = i;
      break;
    }
  }
  assert(decaypath != DECAYPATH_COUNT) // Failed to select pellet

  pkt_ptr->tdecay = -1.; // ensure we enter the following loop

  #ifdef NO_INITIAL_PACKETS
  const double tdecaymin = tmin;
  #else
  const double tdecaymin = 0.; // allow decays before the first timestep
  #endif

  if (UNIFORM_PELLET_ENERGIES)
  {
    pkt_ptr->tdecay = sample_decaytime(decaypath, tdecaymin, tmax);
    pkt_ptr->e_cmf = e0;
  }
  else
  {
    // uniform decay time distribution (scale the packet energies instead)
    const double zrand = gsl_rng_uniform(rng);
    pkt_ptr->tdecay = zrand * tdecaymin + (1. - zrand) * tmax;
    pkt_ptr->e_cmf = get_decay_power_density(decaypath, mgi, pkt_ptr->tdecay);
  }

  bool from_positron;
  pkt_ptr->type = get_decay_pellet_type(decaypath, &from_positron); // set the packet tdecay and type
  pkt_ptr->originated_from_positron = from_positron;

  /// Now assign the energy to the pellet.
  const double dopplerfactor = doppler_packetpos(pkt_ptr, tmin);
  pkt_ptr->e_rf = pkt_ptr->e_cmf / dopplerfactor;
  pkt_ptr->trueemissiontype = -1;
}


static void setup_packets(int pktnumberoffset, PKT *pkt)
/// Subroutine that initialises the packets if we start a new simulation.
{
  float cont[MGRID + 1];

  /// The total number of pellets that we want to start with is just
  /// npkts. The total energy of the pellets is given by etot.
  const double etot_tinf = (
    (nucdecayenergy(NUCLIDE_NI56) + nucdecayenergy(NUCLIDE_CO56)) * totmassradionuclide[NUCLIDE_NI56] / nucmass(NUCLIDE_NI56) +
    nucdecayenergy(NUCLIDE_CO56) * totmassradionuclide[NUCLIDE_CO56] / nucmass(NUCLIDE_CO56) +
    (nucdecayenergy(NUCLIDE_NI57) + nucdecayenergy(NUCLIDE_CO57)) * totmassradionuclide[NUCLIDE_NI57] / nucmass(NUCLIDE_NI57) +
    nucdecayenergy(NUCLIDE_CO57) * totmassradionuclide[NUCLIDE_CO57] / nucmass(NUCLIDE_CO57) +
    (nucdecayenergy(NUCLIDE_V48) + nucdecayenergy(NUCLIDE_CR48)) * totmassradionuclide[NUCLIDE_CR48] / nucmass(NUCLIDE_CR48) +
    (nucdecayenergy(NUCLIDE_FE52) + nucdecayenergy(NUCLIDE_MN52)) * totmassradionuclide[NUCLIDE_FE52] / nucmass(NUCLIDE_FE52));
  printout("etot %g (t_0 ot t_inf)\n", etot_tinf);
  printout("decayenergy(NI56), decayenergy(CO56), decayenergy_gamma(CO56): %g, %g, %g\n", nucdecayenergy(NUCLIDE_NI56) / MEV, nucdecayenergy(NUCLIDE_CO56) / MEV, nucdecayenergygamma(NUCLIDE_CO56) / MEV);
  printout("decayenergy(NI57), decayenergy_gamma(NI57), nucdecayenergy(CO57): %g, %g, %g\n", nucdecayenergy(NUCLIDE_NI57) / MEV, nucdecayenergygamma(NUCLIDE_NI57) / MEV, nucdecayenergy(NUCLIDE_CO57) / MEV);
  printout("decayenergy(CR48), decayenergy(V48): %g %g\n", nucdecayenergy(NUCLIDE_CR48) / MEV, nucdecayenergy(NUCLIDE_V48) / MEV);
  printout("decayenergy(FE52), decayenergy(MN52): %g %g\n", nucdecayenergy(NUCLIDE_FE52) / MEV, nucdecayenergy(NUCLIDE_MN52) / MEV);

  const double e0_tinf = etot_tinf / npkts / n_out_it / n_middle_it;
  printout("packet e0 (t_0 to t_inf) %g erg\n", e0_tinf);

  // scale up the radioactive abundances to account for the missing masses in the model cells that are not associated with any propagation cells
  for (int iso = 0; iso < RADIONUCLIDE_COUNT; iso++)
  {
    if (totmassradionuclide[iso] <= 0)
      continue;
    double totmassradionuclide_actual = 0.;
    for (int mgi = 0; mgi < npts_model; mgi++)
    {
      totmassradionuclide_actual += get_modelinitradioabund(mgi, iso) * get_rhoinit(mgi) * vol_init_modelcell(mgi);
    }
    const double ratio = totmassradionuclide[iso] / totmassradionuclide_actual;
    if (totmassradionuclide_actual <= 0.)
      continue;
    // printout("nuclide %d ratio %g\n", iso, ratio);
    for (int mgi = 0; mgi < npts_model; mgi++)
    {
      if (get_numassociatedcells(mgi) > 0)
      {
        const double prev_abund = get_modelinitradioabund(mgi, iso);
        const double new_abund = prev_abund * ratio;
        set_modelinitradioabund(mgi, iso, new_abund);
      }
    }
  }

  double modelcell_decay_energy_density[npts_model];
  for (int mgi = 0; mgi < npts_model; mgi++)
  {
    modelcell_decay_energy_density[mgi] = 0.;
    for (int i = 0; i < DECAYPATH_COUNT; i++)
    {
      modelcell_decay_energy_density[mgi] += get_rhoinit(mgi) * get_simtime_endecay_per_ejectamass(mgi, i) * MH;
    }
  }

  // Need to get a normalisation factor.
  float norm = 0.0;
  for (int m = 0; m < ngrid; m++)
  {
    cont[m] = norm;
    const int mgi = cell[m].modelgridindex;

    norm += vol_init_gridcell(m) * modelcell_decay_energy_density[mgi];
  }
  cont[ngrid] = norm;

  const double etot = norm / MH;
  /// So energy per pellet is
  const double e0 = etot / npkts / n_out_it / n_middle_it;
  printout("packet e0 (in time range) %g\n", e0);

  printout("etot %g erg (in time range)\n", etot);

  // Now place the pellets in the ejecta and decide at what time they will decay.

  if (npkts > MPKTS)
  {
    printout("Too many packets. Abort.\n");
    abort();
  }

  printout("Placing pellets...\n");
  for (int n = 0; n < npkts; n++)
  {
    // Get random number.
    int mabove = ngrid;
    int mbelow = 0;
    double zrand = gsl_rng_uniform(rng);

    while (mabove != (mbelow + 1))
    {
      int m;

      if (mabove != (mbelow + 2))
        m = (mabove + mbelow) / 2;
      else
        m = mbelow + 1;

      if (cont[m] > (zrand * norm))
        mabove = m;
      else
        mbelow = m;
    }

    if (cont[mbelow] > (zrand * norm))
    {
      printout("mbelow %d cont[mbelow] %g zrand*norm %g\n", mbelow, cont[mbelow], zrand*norm);
      abort();
    }
    if ((cont[mabove] < (zrand * norm)) && (mabove != ngrid))
    {
      printout("mabove %d cont[mabove] %g zrand*norm %g\n", mabove, cont[mabove], zrand*norm);
      abort();
    }

    const int cellindex = mbelow;
    //printout("chosen cell %d (%d, %g, %g)\n", m, ngrid, zrand, norm);
    //abort();
    /*
    m=0;
    double runtot = 0.0;
    while (runtot < (zrand))
    {
      grid_ptr = &cell[m];
      runtot += grid_ptr->rho_init * f56ni(grid_ptr) * vol_init() / norm;
      m++;
    }
    m = m - 1;
    */

    assert(cellindex < ngrid);

    place_pellet(e0, cellindex, n + pktnumberoffset, &pkt[n]);
  }

  double e_cmf_total = 0.;
  for (int n = 0; n < npkts; n++)
  {
    pkt[n].interactions = 0;
    e_cmf_total += pkt[n].e_cmf;
  }
  const double e_ratio = etot / e_cmf_total;
  printout("packet energy sum %g should be %g normalisation factor: %g\n", e_cmf_total, etot, e_ratio);
  e_cmf_total *= e_ratio;
  for (int n = 0; n < npkts; n++)
  {
    pkt[n].e_cmf *= e_ratio;
    pkt[n].e_rf *= e_ratio;
  }
  printout("radioactive energy which will be freed during simulation time: %g erg\n", e_cmf_total);
}


void packet_init(int middle_iteration, int my_rank, PKT *pkt)
{
  printout("mem_usage: packets occupy %.1f MB\n", MPKTS * (sizeof(PKT *) + sizeof(PKT)) / 1024. / 1024.);
  if (simulation_continued_from_saved)
    return;

  const int pktnumberoffset = middle_iteration * npkts;
  setup_packets(pktnumberoffset, pkt);

  /* Consistency check to debug read/write
  PKT testpkt[MPKTS];
  int n;
  if (i > 0)
  {
  sprintf(filename,"packets%d_%d.tmp",i-1,my_rank);
  packets_file = fopen_required(filename, "rb");
  fread(&testpkt, sizeof(PKT), npkts, packets_file);
  fclose(packets_file);

  for (n=0; n<npkts; n++)
  {
  if (testpkt[n].number-pkt[n].number != 0) printout("fatal 1!!! %d, %d\n",testpkt[n].number,pkt[n].number);
  if (testpkt[n].where-pkt[n].where != 0) printout("fatal 2!!!\n");      /// The grid cell that the packet is in.
  if (testpkt[n].type-pkt[n].type != 0) printout("fatal 3!!!\n");       /// Identifies the type of packet (k-, r-, etc.)
  if (testpkt[n].pos[0]-pkt[n].pos[0] != 0) printout("fatal 4!!!\n");  /// Position of the packet (x,y,z).
  if (testpkt[n].pos[1]-pkt[n].pos[1] != 0) printout("fatal 5!!!\n");  /// Position of the packet (x,y,z).
  if (testpkt[n].pos[2]-pkt[n].pos[2] != 0) printout("fatal 6 !!!\n");  /// Position of the packet (x,y,z).
  if (testpkt[n].dir[0]-pkt[n].dir[0] != 0) printout("fatal 7!!!\n");  /// Direction of propagation. (x,y,z). Always a unit vector.
  if (testpkt[n].dir[1]-pkt[n].dir[1] != 0) printout("fatal 8!!!\n");  /// Direction of propagation. (x,y,z). Always a unit vector.
  if (testpkt[n].dir[2]-pkt[n].dir[2] != 0) printout("fatal 9!!!\n");  /// Direction of propagation. (x,y,z). Always a unit vector.
  if (testpkt[n].last_cross-pkt[n].last_cross != 0) printout("fatal 10!!!\n"); /// To avoid rounding errors on cell crossing.
  if (testpkt[n].tdecay-pkt[n].tdecay != 0) printout("fatal 11!!!\n");  /// Time at which pellet decays.
  if (testpkt[n].e_cmf-pkt[n].e_cmf != 0) printout("fatal 12!!!\n");   /// The energy the packet carries in the co-moving frame.
  if (testpkt[n].e_rf-pkt[n].e_rf != 0) printout("fatal 13!!!\n");    /// The energy the packet carries in the rest frame.
  if (testpkt[n].nu_cmf-pkt[n].nu_cmf != 0) printout("fatal 14!!!\n");  /// The frequency in the co-moving frame.
  if (testpkt[n].nu_rf-pkt[n].nu_rf != 0) printout("fatal 15!!!\n");   /// The frequency in the rest frame.
  if (testpkt[n].escape_type-pkt[n].escape_type != 0) printout("fatal 16!!!\n"); /// Flag to tell us in which form it escaped from the grid.
  if (testpkt[n].escape_time-pkt[n].escape_time != 0) printout("fatal 17!!!\n");
  if (testpkt[n].scat_count-pkt[n].scat_count != 0) printout("fatal 18!!!\n");  /// WHAT'S THAT???
  if (testpkt[n].next_trans-pkt[n].next_trans != 0) printout("fatal 19!!!\n");
  if (testpkt[n].interactions-pkt[n].interactions != 0) printout("fatal 20!!!\n");/// debug: number of interactions the packet undergone
  if (testpkt[n].last_event-pkt[n].last_event != 0) printout("fatal 21!!!\n");  /// debug: stores information about the packets history
}
}
  */
}


void write_packets(char filename[], PKT *pkt)
{
  FILE *packets_file = fopen_required(filename, "w");
  for (int i = 0; i < npkts; i++)
  {
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
    fprintf(packets_file, "%d ", pkt[i].escape_time);
    fprintf(packets_file, "%d ", pkt[i].scat_count);
    fprintf(packets_file, "%d ", pkt[i].next_trans);
    fprintf(packets_file, "%d ", pkt[i].interactions);
    fprintf(packets_file, "%d ", pkt[i].last_event);
    fprintf(packets_file, "%d ", pkt[i].emissiontype);
    fprintf(packets_file, "%d ", pkt[i].trueemissiontype);
    fprintf(packets_file, "%lg %lg %lg ", pkt[i].em_pos[0], pkt[i].em_pos[1], pkt[i].em_pos[2]);
    fprintf(packets_file, "%d ", pkt[i].absorptiontype);
    fprintf(packets_file, "%lg ", pkt[i].absorptionfreq);
    fprintf(packets_file, "%d ", pkt[i].nscatterings);
    fprintf(packets_file, "%d ", pkt[i].em_time);
    fprintf(packets_file, "%lg %lg %lg ", pkt[i].absorptiondir[0], pkt[i].absorptiondir[1], pkt[i].absorptiondir[2]);
    fprintf(packets_file, "%lg %lg %lg ", pkt[i].stokes[0], pkt[i].stokes[1], pkt[i].stokes[2]);
    fprintf(packets_file, "%lg %lg %lg ", pkt[i].pol_dir[0], pkt[i].pol_dir[1], pkt[i].pol_dir[2]);
    fprintf(packets_file, "%d ", pkt[i].originated_from_positron);
    fprintf(packets_file, "%g ", pkt[i].trueemissionvelocity);
    fprintf(packets_file, "\n");
  }
  fclose(packets_file);
}


void read_packets(char filename[], PKT *pkt)
{
  FILE *packets_file = fopen_required(filename, "r");
  char *line = (char *) malloc(sizeof(char) * 4096);

  int packets_read = 0;
  while (!feof(packets_file))
  {
    if (line != fgets(line, 4096, packets_file))
      break;

    packets_read++;
    const int i = packets_read - 1;

    if (i > npkts - 1)
    {
      printout("ERROR: More data found beyond packet %d (expecting %d packets). Recompile exspec with the correct number of packets. Run (wc -l < packets00_0000.out) to count them.\n",
               packets_read, npkts);
      abort();
    }

    char *linepos = line;
    int offset = 0;

    int pkt_type_in;
    sscanf(linepos, "%d %d %d%n", &pkt[i].number, &pkt[i].where, &pkt_type_in, &offset);
    linepos += offset;
    pkt[i].type = (enum packet_type) pkt_type_in;

    sscanf(linepos, "%lg %lg %lg%n", &pkt[i].pos[0], &pkt[i].pos[1], &pkt[i].pos[2], &offset);
    linepos += offset;

    sscanf(linepos, "%lg %lg %lg%n", &pkt[i].dir[0], &pkt[i].dir[1], &pkt[i].dir[2], &offset);
    linepos += offset;

    int last_cross_in;
    sscanf(linepos, "%d%n", &last_cross_in, &offset);
    linepos += offset;
    pkt[i].last_cross = (enum cell_boundary) last_cross_in;

    sscanf(linepos, "%lg%n", &pkt[i].tdecay, &offset);
    linepos += offset;

    sscanf(linepos, "%lg %lg %lg %lg%n", &pkt[i].e_cmf, &pkt[i].e_rf, &pkt[i].nu_cmf, &pkt[i].nu_rf, &offset);
    linepos += offset;

    int escape_type;
    sscanf(linepos, "%d %d %d%n", &escape_type, &pkt[i].escape_time, &pkt[i].scat_count, &offset);
    linepos += offset;
    pkt[i].escape_type = (enum packet_type) escape_type;

    sscanf(linepos, "%d %d %d%n", &pkt[i].next_trans, &pkt[i].interactions, &pkt[i].last_event, &offset);
    linepos += offset;

    sscanf(linepos, "%d %d%n", &pkt[i].emissiontype, &pkt[i].trueemissiontype, &offset);
    linepos += offset;

    sscanf(linepos, "%lg %lg %lg%n", &pkt[i].em_pos[0], &pkt[i].em_pos[1], &pkt[i].em_pos[2], &offset);
    linepos += offset;

    sscanf(linepos, "%d %lg %d%n", &pkt[i].absorptiontype, &pkt[i].absorptionfreq, &pkt[i].nscatterings, &offset);
    linepos += offset;

    sscanf(linepos, "%d%n", &pkt[i].em_time, &offset);
    linepos += offset;

    sscanf(linepos, "%lg %lg %lg%n", &pkt[i].absorptiondir[0], &pkt[i].absorptiondir[1], &pkt[i].absorptiondir[2], &offset);
    linepos += offset;

    sscanf(linepos, "%lg %lg %lg%n", &pkt[i].stokes[0], &pkt[i].stokes[1], &pkt[i].stokes[2], &offset);
    linepos += offset;

    sscanf(linepos, "%lg %lg %lg%n", &pkt[i].pol_dir[0], &pkt[i].pol_dir[1], &pkt[i].pol_dir[2], &offset);
    linepos += offset;

    int int_originated_from_positron;
    sscanf(linepos, "%d%n", &int_originated_from_positron, &offset);
    linepos += offset;
    pkt[i].originated_from_positron = (int_originated_from_positron != 0);

    sscanf(linepos, "%g%n", &pkt[i].trueemissionvelocity, &offset);
    linepos += offset;

  }

  if (packets_read < npkts)
  {
    printout("ERROR: Read failed after packet %d (expecting %d packets). Recompile exspec with the correct number of packets. Run (wc -l < packets00_0000.out) to count them.\n",
             packets_read, npkts);
    abort();
  }
  free(line);
  fclose(packets_file);
}
