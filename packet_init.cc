#include "sn3d.h"
#include "grid.h"
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
  pkt_ptr->prop_time = globals::tmin;
  pkt_ptr->originated_from_positron = false;

  if (globals::grid_type == GRID_SPHERICAL1D)
  {
    const double zrand3 = gsl_rng_uniform(rng);
    const double r_inner = get_cellcoordmin(cellindex, 0);
    const double r_outer = get_cellcoordmin(cellindex, 0) + wid_init(cellindex);
    const double radius = pow(zrand3 * pow(r_inner, 3) + (1. - zrand3) * pow(r_outer, 3), 1/3.);
    // assert_always(radius >= r_inner);
    // assert_always(radius <= r_outer);

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

  const int mgi = get_cell_modelgridindex(cellindex);

  decay::setup_radioactive_pellet(e0, mgi, pkt_ptr);

  // initial e_rf is probably never needed (e_rf is set at pellet decay time), but we
  // might as well give it a correct value since this code is fast and runs only once

  // pellet packet is moving with the homologous flow, so dir is proportional to pos
  vec_norm(pkt_ptr->pos, pkt_ptr->dir);  // assign dir = pos / vec_len(pos)
  const double dopplerfactor = doppler_packetpos(pkt_ptr);
  pkt_ptr->e_rf = pkt_ptr->e_cmf / dopplerfactor;

  pkt_ptr->prop_time = globals::tmin;
  pkt_ptr->trueemissiontype = -1;
}


void packet_init(int middle_iteration, int my_rank, PKT *pkt)
/// Subroutine that initialises the packets if we start a new simulation.
{
  printout("UNIFORM_PELLET_ENERGIES is %s\n", (UNIFORM_PELLET_ENERGIES ? "true" : "false"));

  const int pktnumberoffset = middle_iteration * globals::npkts;
  double cont[globals::ngrid + 1];

  /// The total number of pellets that we want to start with is just
  /// npkts. The total energy of the pellets is given by etot.
  const double etot_tinf = decay::get_global_etot_t0_tinf();

  printout("etot %g (t_0 to t_inf)\n", etot_tinf);

  const double e0_tinf = etot_tinf / globals::npkts / globals::n_out_it / globals::n_middle_it;
  printout("packet e0 (t_0 to t_inf) %g erg\n", e0_tinf);

  decay::setup_chain_energy_per_mass();

  // Need to get a normalisation factor.
  double norm = 0.0;
  for (int m = 0; m < globals::ngrid; m++)
  {
    cont[m] = norm;
    const int mgi = get_cell_modelgridindex(m);

    norm += vol_init_gridcell(m) * get_rhoinit(mgi) * decay::get_modelcell_endecay_per_mass(mgi);
  }
  assert_always(norm > 0);
  cont[globals::ngrid] = norm;

  const double etot = norm;
  /// So energy per pellet is
  const double e0 = etot / globals::npkts / globals::n_out_it / globals::n_middle_it;
  printout("packet e0 (in time range) %g erg\n", e0);

  printout("etot %g erg (in time range) erg\n", etot);

  // Now place the pellets in the ejecta and decide at what time they will decay.

  if (globals::npkts > MPKTS)
  {
    printout("Too many packets. Abort.\n");
    abort();
  }

  printout("Placing pellets...\n");
  for (int n = 0; n < globals::npkts; n++)
  {
    // Get random number.
    int mabove = globals::ngrid;
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
    if ((cont[mabove] < (zrand * norm)) && (mabove != globals::ngrid))
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
      grid_ptr = &globals::cell[m];
      runtot += grid_ptr->rho_init * f56ni(grid_ptr) * vol_init() / norm;
      m++;
    }
    m = m - 1;
    */

    assert_always(cellindex < globals::ngrid);

    place_pellet(e0, cellindex, n + pktnumberoffset, &pkt[n]);
  }

  decay::free_chain_energy_per_mass(); // will no longer be needed after packets are set up

  double e_cmf_total = 0.;
  for (int n = 0; n < globals::npkts; n++)
  {
    pkt[n].interactions = 0;
    e_cmf_total += pkt[n].e_cmf;
  }
  const double e_ratio = etot / e_cmf_total;
  printout("packet energy sum %g should be %g normalisation factor: %g\n", e_cmf_total, etot, e_ratio);
  assert_always(std::isfinite(e_cmf_total));
  e_cmf_total *= e_ratio;
  for (int n = 0; n < globals::npkts; n++)
  {
    pkt[n].e_cmf *= e_ratio;
    pkt[n].e_rf *= e_ratio;
  }
  printout("radioactive energy which will be freed during simulation time: %g erg\n", e_cmf_total);
}


void write_packets(char filename[], PKT *pkt)
{
  FILE *packets_file = fopen_required(filename, "w");
  for (int i = 0; i < globals::npkts; i++)
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
    fprintf(packets_file, "%d ", pkt[i].trueem_time);
    fprintf(packets_file, "%d ", pkt[i].pellet_nucindex);
    fprintf(packets_file, "\n");
  }
  fclose(packets_file);
}


void read_temp_packetsfile(const int timestep, const int my_rank, PKT *const pkt)
{
  char filename[100];
  sprintf(filename, "packets_%.4d_ts%d.tmp", my_rank, timestep);

  printout("Reading %s...", filename);
  FILE *packets_file = fopen_required(filename, "rb");
  fread(pkt, sizeof(PKT), globals::npkts, packets_file);
  //read_packets(packets_file);
  fclose(packets_file);
  printout("done\n");
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

    if (i > globals::npkts - 1)
    {
      printout("ERROR: More data found beyond packet %d (expecting %d packets). Recompile exspec with the correct number of packets. Run (wc -l < packets00_0000.out) to count them.\n",
               packets_read, globals::npkts);
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

    sscanf(linepos, "%d%n", &pkt[i].trueem_time, &offset);
    linepos += offset;

    sscanf(linepos, "%d%n", &pkt[i].pellet_nucindex, &offset);
    linepos += offset;
  }

  if (packets_read < globals::npkts)
  {
    printout("ERROR: Read failed after packet %d (expecting %d packets). Recompile exspec with the correct number of packets. Run (wc -l < packets00_0000.out) to count them.\n",
             packets_read, globals::npkts);
    abort();
  }
  free(line);
  fclose(packets_file);
}
