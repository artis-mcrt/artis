#include "sn3d.h"
#include "grid_init.h"
#include "packet_init.h"
#include "vectors.h"


static double uniform_ni56(const int cellindex)
{
  const double r_on_rmax = get_cellradialpos(cellindex) / rmax;
  const double m_r = pow(r_on_rmax, 3) * mtot / MSUN; //this is the mass enclosed up to radius r in units of the total eject mass

  if (m_r < 0.5)
    return 1.0;
  else if (m_r < 0.75)
    return (1.0 - ((m_r - 0.5) * 4));
  else
    return 0.0;
}


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
  const double f56ni = (model_type == RHO_UNIFORM) ? uniform_ni56(cellindex) : get_modelinitradioabund(mgi, NUCLIDE_NI56);

  // first choose which of the decay chains to sample
  const int nchains = 10;
  double prob_chain[nchains];
  prob_chain[0] = f56ni / MNI56 * E56NI;
  prob_chain[1] = f56ni / MNI56 * E56CO;
  prob_chain[2] = get_modelinitradioabund(mgi, NUCLIDE_FE52) / MFE52 * E52FE;
  prob_chain[3] = get_modelinitradioabund(mgi, NUCLIDE_FE52) / MFE52 * E52MN;
  prob_chain[4] = get_modelinitradioabund(mgi, NUCLIDE_CR48) / MCR48 * E48CR;
  prob_chain[5] = get_modelinitradioabund(mgi, NUCLIDE_CR48) / MCR48 * E48V;
  prob_chain[6] = get_modelinitradioabund(mgi, NUCLIDE_CO56) / MCO56 * E56CO;
  prob_chain[7] = get_modelinitradioabund(mgi, NUCLIDE_NI57) / MNI57 * E57NI;
  prob_chain[8] = get_modelinitradioabund(mgi, NUCLIDE_NI57) / MNI57 * E57CO;
  prob_chain[9] = get_modelinitradioabund(mgi, NUCLIDE_CO57) / MCO57 * E57CO;

  double prob_sum = 0.;
  for (int i = 0; i < nchains; i++)
    prob_sum += prob_chain[i];

  const double zrand_chain = gsl_rng_uniform(rng) * prob_sum;
  double prob_accumulated = 0.;
  int selected_chain = -1;
  for (int i = 0; i < nchains; i++)
  {
    prob_accumulated += prob_chain[i];
    if (zrand_chain <= prob_accumulated)
    {
      selected_chain = i;
      break;
    }
  }
  if (selected_chain < 0)
  {
    printout("Failed to select pellet\n");
    abort();
  }

  double zrand = gsl_rng_uniform(rng);
  double zrand2;
  switch (selected_chain)
  {
    case 0:  // Ni56 pellet
      pkt_ptr->type = TYPE_56NI_PELLET;
      pkt_ptr->tdecay = -T56NI * log(zrand);
      break;

    case 1: // Ni56 -> Co56 pellet
      if (zrand < E56CO_GAMMA / E56CO)
      {
        pkt_ptr->type = TYPE_56CO_PELLET;
      }
      else
      {
        pkt_ptr->type = TYPE_56CO_POSITRON_PELLET;
        pkt_ptr->originated_from_positron = true;
      }

      zrand = gsl_rng_uniform(rng);
      zrand2 = gsl_rng_uniform(rng);
      pkt_ptr->tdecay = (-T56NI * log(zrand)) + (-T56CO * log(zrand2));
      break;

    case 2: // Fe52 pellet
      pkt_ptr->type = TYPE_52FE_PELLET;
      pkt_ptr->tdecay = -T52FE * log(zrand);
      break;

    case 3: // Fe52 -> Mn52 pellet
      pkt_ptr->type = TYPE_52MN_PELLET;
      zrand2 = gsl_rng_uniform(rng);
      pkt_ptr->tdecay = (-T52FE * log(zrand)) + (-T52MN * log(zrand2));
      break;

    case 4: // Cr48 pellet
      pkt_ptr->type = TYPE_48CR_PELLET;
      pkt_ptr->tdecay = -T48CR * log(zrand);
      break;

    case 5: // Cr48 -> V48 pellet
      pkt_ptr->type = TYPE_48V_PELLET;
      zrand2 = gsl_rng_uniform(rng);
      pkt_ptr->tdecay = (-T48CR * log(zrand)) + (-T48V * log(zrand2));
      break;

    case 6: // Co56 pellet
      /// Now it is a 56Co pellet, choose whether it becomes a positron
      if (zrand < E56CO_GAMMA / E56CO)
      {
        pkt_ptr->type = TYPE_56CO_PELLET;
      }
      else
      {
        pkt_ptr->type = TYPE_56CO_POSITRON_PELLET;
        pkt_ptr->originated_from_positron = true;
      }

      zrand = gsl_rng_uniform(rng);
      pkt_ptr->tdecay = -T56CO * log(zrand);
      break;

    case 7: // Ni57 pellet
      if (zrand < E57NI_GAMMA / E57NI)
      {
        pkt_ptr->type = TYPE_57NI_PELLET;
      }
      else
      {
        pkt_ptr->type = TYPE_57NI_POSITRON_PELLET;
        pkt_ptr->originated_from_positron = true;
      }

      zrand = gsl_rng_uniform(rng);
      pkt_ptr->tdecay = -T57NI * log(zrand);
      break;

    case 8: // Ni57 -> Co57 pellet
      pkt_ptr->type = TYPE_57CO_PELLET;
      zrand2 = gsl_rng_uniform(rng);
      pkt_ptr->tdecay = (-T57NI * log(zrand)) + (-T57CO * log(zrand2));
      break;

    case 9: // Co57 pellet
      pkt_ptr->type = TYPE_57CO_PELLET;
      pkt_ptr->tdecay = -T57CO * log(zrand);
      break;

    default:
      printout("Problem selecting pellet type\n");
      abort();
  }

  /// Now assign the energy to the pellet.
  pkt_ptr->e_cmf = e0;
  pkt_ptr->trueemissiontype = -1;
}


static void setup_packets(int pktnumberoffset, PKT *pkt)
/// Subroutine that initialises the packets if we start a new simulation.
{
  float cont[MGRID + 1];

  /// The total number of pellets that we want to start with is just
  /// npkts. The total energy of the pellets is given by etot.
  const double etot = (
    (E56NI + E56CO) * totmassradionuclide[NUCLIDE_NI56] / MNI56 +
    E56CO * totmassradionuclide[NUCLIDE_CO56] / MCO56 +
    (E57NI + E57CO) * totmassradionuclide[NUCLIDE_NI57] / MNI57 +
    E57CO * totmassradionuclide[NUCLIDE_CO57] / MCO57 +
    (E48V + E48CR) * totmassradionuclide[NUCLIDE_CR48] / MCR48 +
    (E52FE + E52MN) * totmassradionuclide[NUCLIDE_FE52] / MFE52);
  printout("etot %g\n", etot);
  printout("E56NI, E56CO, E56CO_GAMMA: %g, %g, %g\n", E56NI / MEV, E56CO / MEV, E56CO_GAMMA / MEV);
  printout("E57NI, E57NI_GAMMA, E57CO: %g, %g, %g\n", E57NI / MEV, E57NI_GAMMA / MEV, E57CO / MEV);
  printout("E48CR, E48V: %g %g\n", E48CR / MEV, E48V / MEV);
  printout("E52FE, E52MN: %g %g\n", E52FE / MEV, E52MN / MEV);

  /// So energy per pellet is
  const double e0 = etot / npkts / n_out_it / n_middle_it;
  printout("e0 %g\n", e0);

  /* Now place the pellets in the ejecta and decide at what time
  they will decay. */

  /* Need to get a normalisation factor. */
  float norm = 0.0;
  for (int m = 0; m < ngrid; m++)
  {
    cont[m] = norm;
    //printf("%g %g %g\n", (f56ni(grid_ptr)*(E56NI + E56CO)/MNI56),(f52fe(grid_ptr)*(E52FE + E52MN)/MFE52),(f48cr(grid_ptr)*(E48V + E48CR)/MCR48));
    const int mgi = cell[m].modelgridindex;
    const double f56ni = (model_type == RHO_UNIFORM) ? uniform_ni56(m) : get_modelinitradioabund(mgi, NUCLIDE_NI56);

    norm += get_rhoinit(mgi) * vol_init_gridcell(m) *
              ((f56ni * (E56NI + E56CO) / 56.)
               + (get_modelinitradioabund(mgi, NUCLIDE_CO56) * (E56CO) / 56.)
               + (get_modelinitradioabund(mgi, NUCLIDE_NI57) * (E57NI + E57CO) / 57.)
               + (get_modelinitradioabund(mgi, NUCLIDE_CO57) * (E57CO) / 57.)
               + (get_modelinitradioabund(mgi, NUCLIDE_FE52) * (E52FE + E52MN) / 52.)
               + (get_modelinitradioabund(mgi, NUCLIDE_CR48) * (E48V + E48CR) / 48.));
  }
  cont[ngrid] = norm;

  if (npkts > MPKTS)
  {
    printout("Too many packets. Abort.\n");
    abort();
  }

  int n = 0;
  int packet_reset = 0;
  printout("Placing pellets...\n");
  while (n < npkts)
  {
    /// Get random number.
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

    if (cellindex >= ngrid)
    {
      printout("Failed to place pellet. Abort.\n");
      abort();
    }

    /// Pellet is in cell m, pointer is grid_ptr.
    //  printout("Pellet in cell %d\n",m);
    place_pellet(e0, cellindex, n + pktnumberoffset, &pkt[n]);

    #ifdef NO_INITIAL_PACKETS
    if (pkt[n].tdecay < tmax && pkt[n].tdecay > tmin)
    #else
    if (pkt[n].tdecay < tmax)
    #endif
      n++;
    else
      packet_reset++;
  }

  /// Some fraction of the packets we reassigned because they were not going
  /// to activate in the time of interest so need to renormalise energies
  /// to account for this.
  for (int n = 0; n < npkts; n++)
  {
    pkt[n].e_cmf = pkt[n].e_cmf * npkts / (npkts + packet_reset);
    pkt[n].interactions = 0;
  }
  printout("radioactive energy which will be freed during simulation time %g\n", etot * npkts / (npkts + packet_reset));
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
  char *line = malloc(sizeof(char) * 4096);

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
    pkt[i].type = pkt_type_in;

    sscanf(linepos, "%lg %lg %lg%n", &pkt[i].pos[0], &pkt[i].pos[1], &pkt[i].pos[2], &offset);
    linepos += offset;

    sscanf(linepos, "%lg %lg %lg%n", &pkt[i].dir[0], &pkt[i].dir[1], &pkt[i].dir[2], &offset);
    linepos += offset;

    int last_cross_in;
    sscanf(linepos, "%d%n", &last_cross_in, &offset);
    linepos += offset;
    pkt[i].last_cross = last_cross_in;

    sscanf(linepos, "%lg%n", &pkt[i].tdecay, &offset);
    linepos += offset;

    sscanf(linepos, "%lg %lg %lg %lg%n", &pkt[i].e_cmf, &pkt[i].e_rf, &pkt[i].nu_cmf, &pkt[i].nu_rf, &offset);
    linepos += offset;

    int escape_type;
    sscanf(linepos, "%d %d %d%n", &escape_type, &pkt[i].escape_time, &pkt[i].scat_count, &offset);
    linepos += offset;
    pkt[i].escape_type = escape_type;

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
