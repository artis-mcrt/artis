#include <stdio.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_vector_double.h>
#include "sn3d.h"
#include "atomic.h"
#include "grid_init.h"
#include "ltepop.h"
#include "macroatom.h"
#include "nonthermal.h"
#include "nltepop.h"
#include "ratecoeff.h"
#include "update_grid.h"

static FILE *nlte_file;


static inline
int get_nlte_vector_index(const int element, const int ion, const int level)
// this is the index for the NLTE solver that is handling all ions of a single element
// This is NOT an index into modelgrid[modelgridindex].nlte_pops that contains all elements
{
  // have to convert from nlte_pops index to nlte_vector index
  // the difference is that nlte vectors apply to a single element and include ground states
  // The (+ ion) term accounts for the ground state population indicies that are not counted in the NLTE array
  const int gs_index = elements[element].ions[ion].first_nlte - elements[element].ions[0].first_nlte + ion;

  // add in level or superlevel number
  const int level_index = gs_index + (is_nlte(element, ion, level) ? level : (get_nlevels_nlte(element, ion) + 1));

  return level_index;
}


static void get_ion_level_of_nlte_vector_index(const int index, const int element, int *ion, int *level)
{
  // this could easily be optimized if need be
  for (int dion = 0; dion < get_nions(element); dion++)
  {
    for (int dlevel = 0; dlevel < get_nlevels(element, dion); dlevel++)
    {
      if (get_nlte_vector_index(element, dion, dlevel) == index)
      {
        *ion = dion;
        *level = dlevel;
        return;
      }
    }
  }
}


static void eliminate_nlte_matrix_rowcol(
  const int index,
  const int gs_index,
  gsl_matrix *restrict rate_matrix,
  gsl_vector *restrict balance_vector)
{
  const gsl_matrix rate_matrix_var = *rate_matrix;

  const int colcount = rate_matrix_var.size2;
  for (int column = 0; column < colcount; column++)
    gsl_matrix_set(rate_matrix, index, column, 0.0);

  const int rowcount = rate_matrix_var.size1;
  for (int row = 1; row < rowcount; row++)
    gsl_matrix_set(rate_matrix, row, index, 0.0);

  gsl_matrix_set(rate_matrix, index, gs_index, -1.0);
  gsl_matrix_set(rate_matrix, index, index, 1.0);
  gsl_vector_set(balance_vector, index, 0.0);
}


static void filter_nlte_matrix(
  const int element,
  gsl_matrix *restrict rate_matrix,
  gsl_vector *restrict balance_vector,
  const gsl_vector *restrict pop_norm_factor_vec)
// find rows and columns that barely interaction with other levels, and effectively
// removing them by zeroing their interactions and setting their departure
// coeff to 1.0
{
  const gsl_matrix rate_matrix_var = *rate_matrix;
  const int nlte_dimension = rate_matrix_var.size1;
  for (int index = 0; index < nlte_dimension; index++)
  {
    double row_max = 0.0;
    for (int column = 0; column < nlte_dimension; column++)
    {
      const double element_value = fabs(gsl_matrix_get(rate_matrix, index, column));
      if (element_value > row_max)
        row_max = element_value;
    }
    double col_max = 0.0;
    for (int row = 1; row < nlte_dimension; row++) //skip the normalisation row 0
    {
      const double element_value = fabs(gsl_matrix_get(rate_matrix, row, index));
      if (element_value > col_max)
        col_max = element_value;
    }
    int ion = -1;
    int level = -1;
    get_ion_level_of_nlte_vector_index(index, element, &ion, &level);
    // printout("index%4d (ion_stage%2d level%4d) row_max %.1e col_max %.1e ",
    //          index,get_ionstage(element,ion),level,row_max,col_max);

    if ((row_max < 1e-100) || (col_max < 1e-100))
    {
      if (level == 0)
      {
        // printout("(Would eliminate but it's a ground state, so keeping it)");
        //printout("(Would eliminate but it's a ground state, so forcing pop=MINPOP=%g)",MINPOP);
        //gsl_vector_set(balance_vector, index, MINPOP / get_vector_get(pop_norm_factor_vec, index));
        //printout("(Eliminating this ground state)");
      }
      else
      {
        double gs_index = get_nlte_vector_index(element,ion, 0);
        eliminate_nlte_matrix_rowcol(index,gs_index, rate_matrix, balance_vector);
        // printout("(forcing LTE population)");
      }
    }
    // printout("\n");
  }
}


static double get_total_rate(
  const int index_selected, const gsl_matrix *restrict rate_matrix, const gsl_vector *restrict popvec,
  const bool into_level, const bool only_levels_below, const bool only_levels_above)
{
  double total_rate = 0.;
  assert(!only_levels_below || !only_levels_above);

  if (into_level)
  {
    // find rate into selected level
    gsl_vector_const_view row_view = gsl_matrix_const_row(rate_matrix, index_selected);
    const gsl_vector *rates_vec = &row_view.vector;

    // multiply incoming rate coefficients by their corresponding populations to get rates
    if (!only_levels_above) // add levels below
    {
      for (int index = 0; index < index_selected; index++)
      {
        total_rate += gsl_vector_get(rates_vec, index) * gsl_vector_get(popvec, index);
      }
    }

    if (!only_levels_below) // add levels above
    {
      for (unsigned int index = index_selected + 1; index < rate_matrix->size1; index++)
      {
        total_rate += gsl_vector_get(rates_vec, index) * gsl_vector_get(popvec, index);
      }
    }
  }
  else
  {
    // find rate out of selected level
    gsl_vector_const_view col_view = gsl_matrix_const_column(rate_matrix, index_selected);
    const gsl_vector *rates_vec = &col_view.vector;

    // multiply outgoing rate coefficients by the population of the selected level to get rates

    if (!only_levels_above) // add levels below
    {
      for (int index = 0; index < index_selected; index++)
      {
        total_rate += gsl_vector_get(rates_vec, index);
      }
    }

    if (!only_levels_below) // add levels above
    {
      for (unsigned int index = index_selected + 1; index < rate_matrix->size2; index++)
      {
        total_rate += gsl_vector_get(rates_vec, index);
      }
    }

    total_rate *= gsl_vector_get(popvec, index_selected);
  }

  return total_rate;
}


static double get_total_rate_in(
  const int index_to, const gsl_matrix *restrict rate_matrix, const gsl_vector *restrict popvec)
{
  return get_total_rate(index_to, rate_matrix, popvec, true, false, false);
}


static double get_total_rate_out(
  const int index_from, const gsl_matrix *restrict rate_matrix, const gsl_vector *restrict popvec)
{
  return get_total_rate(index_from, rate_matrix, popvec, false, false, false);
}



static void print_level_rates_summary(
  const int element, const int selected_ion, const int selected_level,
  const gsl_vector *restrict popvec,
  const gsl_matrix *rate_matrix_rad_bb,
  const gsl_matrix *rate_matrix_coll_bb,
  const gsl_matrix *rate_matrix_ntcoll_bb,
  const gsl_matrix *rate_matrix_rad_bf,
  const gsl_matrix *rate_matrix_coll_bf,
  const gsl_matrix *rate_matrix_ntcoll_bf)
{
  const int selected_index = get_nlte_vector_index(element, selected_ion, selected_level);

  for (int i = 0; i <= 3; i++)
  {
    // rates in from below, in from above, out to below, out to above
    if (i == 0)
    {
      const int nlevels_nlte = get_nlevels_nlte(element, selected_ion);
      if (ion_has_superlevel(element, selected_ion) && (selected_level == nlevels_nlte + 1))
      {
        printout("      superlevel  ");
      }
      else
      {
        printout("    level+1%5d  ", selected_level + 1);
      }
      printout(" %9.1e ", gsl_vector_get(popvec, selected_index));
    }
    else
    {
      printout("                             ");
    }

    const bool into_level = (i <= 1);
    const bool only_levels_below = i % 2;
    const bool only_levels_above = !only_levels_below;

    const double rad_bb_total = get_total_rate(selected_index, rate_matrix_rad_bb, popvec, into_level, only_levels_below, only_levels_above);
    const double coll_bb_total = get_total_rate(selected_index, rate_matrix_coll_bb, popvec, into_level, only_levels_below, only_levels_above);
    const double ntcoll_bb_total = get_total_rate(selected_index, rate_matrix_ntcoll_bb, popvec, into_level, only_levels_below, only_levels_above);
    const double rad_bf_total = get_total_rate(selected_index, rate_matrix_rad_bf, popvec, into_level, only_levels_below, only_levels_above);
    const double coll_bf_total = get_total_rate(selected_index, rate_matrix_coll_bf, popvec, into_level, only_levels_below, only_levels_above);
    const double ntcoll_bf_total = get_total_rate(selected_index, rate_matrix_ntcoll_bf, popvec, into_level, only_levels_below, only_levels_above);

    if (into_level)
    {
      // into this level
      if (only_levels_below)
      {
        printout(" from below ");
      }
      else
      {
        printout(" from above ");
      }
    }
    else
    {
      // out of this level
      if (only_levels_below)
      {
        printout("   to below ");
      }
      else
      {
        printout("   to above ");
      }
    }

    printout("%9.1e %9.1e %9.1e %9.1e %9.1e %9.1e\n",
             rad_bb_total, coll_bb_total, ntcoll_bb_total, rad_bf_total, coll_bf_total, ntcoll_bf_total);
  }
}


static void print_element_rates_summary(
  const int element,
  const gsl_vector *restrict popvec,
  const gsl_matrix *rate_matrix_rad_bb,
  const gsl_matrix *rate_matrix_coll_bb,
  const gsl_matrix *rate_matrix_ntcoll_bb,
  const gsl_matrix *rate_matrix_rad_bf,
  const gsl_matrix *rate_matrix_coll_bf,
  const gsl_matrix *rate_matrix_ntcoll_bf)
{
  const int nions = get_nions(element);
  for (int ion = 0; ion < nions; ion++)
  {
    const int nlevels = get_nlevels(element, ion);
    const int nlevels_nlte = get_nlevels_nlte(element, ion);

    const int atomic_number = get_element(element);
    const int ionstage = get_ionstage(element, ion);

    const int max_printed_levels = ion_has_superlevel(element, ion) ? nlevels_nlte + 2 : nlevels_nlte + 1;

    for (int level = 0; (level < max_printed_levels) && (level < nlevels) && (level <= nlevels_nlte + 1); level++)
    {
      if (level == 0)
      {
        printout("  Z=%2d ion_stage %2d      pop       rates    bb_rad    bb_col  bb_ntcol    bf_rad    bf_col  bf_ntcol\n",
                 atomic_number, ionstage);
      }

      print_level_rates_summary(element, ion, level, popvec, rate_matrix_rad_bb, rate_matrix_coll_bb, rate_matrix_ntcoll_bb, rate_matrix_rad_bf, rate_matrix_coll_bf, rate_matrix_ntcoll_bf);
    }

    if (ion_has_superlevel(element, ion) && max_printed_levels < (nlevels_nlte + 1))
    {
      const int level_superlevel = nlevels_nlte + 1;

      print_level_rates_summary(element, ion, level_superlevel, popvec, rate_matrix_rad_bb, rate_matrix_coll_bb, rate_matrix_ntcoll_bb, rate_matrix_rad_bf, rate_matrix_coll_bf, rate_matrix_ntcoll_bf);
    }
  }
}


static void print_level_rates(
  const int element, const int selected_ion, const int selected_level,
  const gsl_vector *restrict popvec,
  const gsl_matrix *rate_matrix_rad_bb,
  const gsl_matrix *rate_matrix_coll_bb,
  const gsl_matrix *rate_matrix_ntcoll_bb,
  const gsl_matrix *rate_matrix_rad_bf,
  const gsl_matrix *rate_matrix_coll_bf,
  const gsl_matrix *rate_matrix_ntcoll_bf)
{
  // very detailed output of the NLTE processes for a particular levels

  if (element > nelements - 1 || selected_ion > get_nions(element) - 1 || selected_level > (get_nlevels_nlte(element, selected_ion) + (ion_has_superlevel(element, selected_ion) ? 1 : 0)))
  {
    printout("print_level_rates: invalid element/ion/level arguments\n");
    abort();
  }

  if (rate_matrix_rad_bb == rate_matrix_coll_bb)
  {
    printout("print_level_rates: rate_matrix_rad_bb == rate_matrix_coll_bb. check individual_process_matricies is off\n");
    abort();
  }

  const gsl_vector popvector = *popvec;
  const int nlte_dimension = popvector.size;
  const int atomic_number = get_element(element);
  const int selected_ionstage = get_ionstage(element, selected_ion);
  const int selected_index = get_nlte_vector_index(element, selected_ion, selected_level);
  const double pop_selectedlevel = gsl_vector_get(popvec, selected_index);
  printout("NLTE level diagnostics for Z=%d ion_stage %d level %d (vector_index %d): rates into and out of this level\n",
           atomic_number, selected_ionstage, selected_level, selected_index);

  const double rad_bb_in_total = get_total_rate_in(selected_index, rate_matrix_rad_bb, popvec);
  const double coll_bb_in_total = get_total_rate_in(selected_index, rate_matrix_coll_bb, popvec);
  const double ntcoll_bb_in_total = get_total_rate_in(selected_index, rate_matrix_ntcoll_bb, popvec);
  const double rad_bf_in_total = get_total_rate_in(selected_index, rate_matrix_rad_bf, popvec);
  const double coll_bf_in_total = get_total_rate_in(selected_index, rate_matrix_coll_bf, popvec);
  const double ntcoll_bf_in_total = get_total_rate_in(selected_index, rate_matrix_ntcoll_bf, popvec);
  const double total_rate_in = rad_bb_in_total + coll_bb_in_total + rad_bf_in_total + coll_bf_in_total + ntcoll_bf_in_total;
  printout("  TOTAL rates in:             rad_bb_in  %7.1e coll_bb_in  %7.1e ntcoll_bb_in  %7.1e rad_bf_in  %7.1e coll_bf_in  %7.1e ntcoll_bf_in  %7.1e\n",
           rad_bb_in_total, coll_bb_in_total, ntcoll_bb_in_total, rad_bf_in_total, coll_bf_in_total, ntcoll_bf_in_total);

  const double rad_bb_out_total = get_total_rate_out(selected_index, rate_matrix_rad_bb, popvec);
  const double coll_bb_out_total = get_total_rate_out(selected_index, rate_matrix_coll_bb, popvec);
  const double ntcoll_bb_out_total = get_total_rate_out(selected_index, rate_matrix_ntcoll_bb, popvec);
  const double rad_bf_out_total = get_total_rate_out(selected_index, rate_matrix_rad_bf, popvec);
  const double coll_bf_out_total = get_total_rate_out(selected_index, rate_matrix_coll_bf, popvec);
  const double ntcoll_bf_out_total = get_total_rate_out(selected_index, rate_matrix_ntcoll_bf, popvec);
  const double total_rate_out = rad_bb_out_total + coll_bb_out_total + rad_bf_out_total + coll_bf_out_total + ntcoll_bf_out_total;
  printout("  TOTAL rates out:            rad_bb_out %7.1e coll_bb_out %7.1e ntcoll_bb_out %7.1e rad_bf_out %7.1e coll_bf_out %7.1e ntcoll_bf_out %7.1e\n",
          rad_bb_out_total, coll_bb_out_total, ntcoll_bb_out_total, rad_bf_out_total, coll_bf_out_total, ntcoll_bf_out_total);

  for (int index = 0; index < nlte_dimension; index++)
  {
    if (index == selected_index)
      continue;
    int ion;
    int level;
    get_ion_level_of_nlte_vector_index(index, element, &ion, &level);
    const int ionstage = get_ionstage(element, ion);
    // in means populating the selected level, out means depopulating the selected level
    const double pop = gsl_vector_get(popvec, index);
    const double rad_bb_in = gsl_matrix_get(rate_matrix_rad_bb, selected_index, index) * pop;
    const double rad_bb_out = gsl_matrix_get(rate_matrix_rad_bb, index, selected_index) * pop_selectedlevel;
    const double coll_bb_in = gsl_matrix_get(rate_matrix_coll_bb, selected_index, index) * pop;
    const double coll_bb_out = gsl_matrix_get(rate_matrix_coll_bb, index, selected_index) * pop_selectedlevel;
    const double ntcoll_bb_in = gsl_matrix_get(rate_matrix_ntcoll_bb, selected_index, index) * pop;
    const double ntcoll_bb_out = gsl_matrix_get(rate_matrix_ntcoll_bb, index, selected_index) * pop_selectedlevel;
    const double rad_bf_in = gsl_matrix_get(rate_matrix_rad_bf, selected_index, index) * pop;
    const double rad_bf_out = gsl_matrix_get(rate_matrix_rad_bf, index, selected_index) * pop_selectedlevel;
    const double coll_bf_in = gsl_matrix_get(rate_matrix_coll_bf, selected_index, index) * pop;
    const double coll_bf_out = gsl_matrix_get(rate_matrix_coll_bf, index, selected_index) * pop_selectedlevel;
    const double ntcoll_bf_in = gsl_matrix_get(rate_matrix_ntcoll_bf, selected_index, index) * pop;
    const double ntcoll_bf_out = gsl_matrix_get(rate_matrix_ntcoll_bf, index, selected_index) * pop_selectedlevel;

    const bool nonzero_rate_in = (fabs(rad_bb_in) > 0. || fabs(coll_bb_in) > 0. || fabs(ntcoll_bb_in) > 0. || fabs(rad_bf_in) > 0. || fabs(coll_bf_in) > 0. || fabs(ntcoll_bf_in) > 0.);
    const bool nonzero_rate_out = (fabs(rad_bb_out) > 0. || fabs(coll_bb_out) > 0. || fabs(ntcoll_bb_out) > 0. ||fabs(rad_bf_out) > 0. || fabs(coll_bf_out) > 0. || fabs(ntcoll_bf_out) > 0.);
    if (nonzero_rate_in || nonzero_rate_out)
    {
      const double epsilon_trans = fabs(epsilon(element, ion, level) - epsilon(element, selected_ion, selected_level));
      const double nu_trans = epsilon_trans / H;
      const double lambda = 1e8 * CLIGHT / nu_trans; // should be in Angstroms
      const double level_rate_in = rad_bb_in + coll_bb_in + ntcoll_bb_in + rad_bf_in + coll_bf_in + ntcoll_bf_in;
      const double level_rate_out = rad_bb_out + coll_bb_out + ntcoll_bb_in + rad_bf_out + coll_bf_out + ntcoll_bf_out;
      const int level_percent_in = round(level_rate_in / total_rate_in * 100);
      const int level_percent_out = round(level_rate_out / total_rate_out * 100);

      printout("  ionstage %d level %4d (%3d%% of in)  rad_bb_in  %7.1e coll_bb_in  %7.1e ntcoll_bb_in  %7.1e rad_bf_in  %7.1e coll_bf_in  %7.1e ntcoll_bf_in  %7.1e lambda %5.0f\n",
               ionstage, level, level_percent_in, rad_bb_in, coll_bb_in, ntcoll_bb_in, rad_bf_in, coll_bf_in, ntcoll_bf_in, lambda);
      printout("  ionstage %d level %4d (%3d%% of out) rad_bb_out %7.1e coll_bb_out %7.1e ntcoll_bb_out %7.1e rad_bf_out %7.1e coll_bf_out %7.1e ntcoll_bf_out %7.1e lambda %5.0f\n",
               ionstage, level, level_percent_out, rad_bb_out, coll_bb_out, ntcoll_bb_out, rad_bf_out, coll_bf_out, ntcoll_bf_out, lambda);
    }
  }
  printout("\n");
}


static void nltepop_reset_element(const int modelgridindex, const int element)
{
  const int nions = get_nions(element);
  for (int ion = 0; ion < nions; ion++)
  {
    const int nlte_start = elements[element].ions[ion].first_nlte;
    const int nlevels_nlte = get_nlevels_nlte(element,ion);
    for (int level = 1; level < nlevels_nlte; level++)
    {
      modelgrid[modelgridindex].nlte_pops[nlte_start + level - 1] = -1.0; // flag to indicate no useful data
    }

    if (ion_has_superlevel(element, ion))
    {
      modelgrid[modelgridindex].nlte_pops[nlte_start + nlevels_nlte] = -1.0;
    }
  }
}


static int get_element_nlte_dimension_and_slpartfunc(
  const int modelgridindex, const int element, const int nions, double *restrict superlevel_partfunc)
{
  int nlte_dimension = 0;
  for (int ion = 0; ion < nions; ion++)
  {
    superlevel_partfunc[ion] = 0.;
    const int nlevels_nlte = get_nlevels_nlte(element, ion);

    //this is the total number of nlte_levels (i.e. the size of the
    //storage). Our rate matrix will need to be of this dimension +2: the
    //ground state, the "super level".
    //If there's no super level needed then we only need +1
    if (ion_has_superlevel(element, ion))
    {
      nlte_dimension += nlevels_nlte + 2;
      const int nlevels = get_nlevels(element, ion);
      for (int level = nlevels_nlte + 1; level < nlevels; level++)
      {
        superlevel_partfunc[ion] += superlevel_boltzmann(modelgridindex, element, ion, level);
      }
      // printout("  NLTE: including ion_stage %d, which contributes %d to the vector dimension (including superlevel with partfunc %g)\n",
      //          get_ionstage(element, ion), nlevels_nlte + 2, superlevel_partfunc[ion]);
    }
    else
    {
      nlte_dimension += nlevels_nlte + 1;
      // printout("  NLTE: including ion_stage %d, which contributes %d to the vector dimension (no super level)\n",
      //          get_ionstage(element, ion), nlevels_nlte + 1);
    }
  }

  return nlte_dimension;
}


static void nltepop_matrix_add_boundbound(const int modelgridindex, const int element, const int ion,
                                          const double t_mid, double *restrict s_renorm,
                                          gsl_matrix *rate_matrix_rad_bb,
                                          gsl_matrix *rate_matrix_coll_bb,
                                          gsl_matrix *rate_matrix_ntcoll_bb)
{
  const float T_e = get_Te(modelgridindex);
  const float nne = get_nne(modelgridindex);
  const int nlevels = get_nlevels(element, ion);
  // const int Z = get_element(element);
  // const int ionstage = get_ionstage(element, ion);
  for (int level = 0; level < nlevels; level++)
  {
    const int level_index = get_nlte_vector_index(element, ion, level);
    const double epsilon_level = epsilon(element, ion, level);

    // de-excitation
    const int ndowntrans = get_ndowntrans(element, ion, level);
    for (int i = 0; i < ndowntrans; i++)
    {
      const int lineindex = elements[element].ions[ion].levels[level].downtrans_lineindicies[i];
      const int lower = linelist[lineindex].lowerlevelindex;

      const double epsilon_trans = epsilon_level - epsilon(element, ion, lower);

      const double R = rad_deexcitation_ratecoeff(modelgridindex, element, ion, level, lower, epsilon_trans, lineindex, t_mid) * s_renorm[level];
      const double C = col_deexcitation_ratecoeff(T_e, nne, epsilon_trans, lineindex) * s_renorm[level];

      const int upper_index = level_index;
      const int lower_index = get_nlte_vector_index(element,ion,lower);

      *gsl_matrix_ptr(rate_matrix_rad_bb, upper_index, upper_index) -= R;
      *gsl_matrix_ptr(rate_matrix_rad_bb, lower_index, upper_index) += R;
      *gsl_matrix_ptr(rate_matrix_coll_bb, upper_index, upper_index) -= C;
      *gsl_matrix_ptr(rate_matrix_coll_bb, lower_index, upper_index) += C;
      if ((R < 0) || (C < 0))
        printout("  WARNING: Negative de-excitation rate from ion_stage %d level %d to level %d\n", get_ionstage(element, ion), level, lower);
    }

    // excitation
    const int nuptrans = get_nuptrans(element, ion, level);
    for (int i = 0; i < nuptrans; i++)
    {
      const int lineindex = elements[element].ions[ion].levels[level].uptrans_lineindicies[i];
      const int upper = linelist[lineindex].upperlevelindex;
      const double epsilon_trans = epsilon(element, ion, upper) - epsilon_level;

      const double R = rad_excitation_ratecoeff(modelgridindex, element, ion, level, upper, epsilon_trans, lineindex, t_mid) * s_renorm[level];
      const double C = col_excitation_ratecoeff(T_e, nne, lineindex, epsilon_trans) * s_renorm[level];
      // const double NTC = nt_excitation_ratecoeff(modelgridindex, element, ion, level, upper, epsilon_trans, lineindex) * s_renorm[level];
      const double NTC = 0.;

      // if ((Z == 26) && (ionstage == 1) && (level == 0) && (upper <= 5))
      // {
      //   const double tau_sobolev = get_tau_sobolev(modelgridindex, lineindex, t_mid);
      //   const double nu_trans = epsilon_trans / H;
      //   const double lambda = 1e8 * CLIGHT / nu_trans; // should be in Angstroms
      //   printout("Z=%d ionstage %d lower %d upper %d lambda %6.1fÅ tau_sobolev=%g einstein_A %g osc_strength %g coll_str %g\n",
      //            Z, ionstage, level, upper, lambda, tau_sobolev,
      //            linelist[lineindex].einstein_A, linelist[lineindex].osc_strength, linelist[lineindex].coll_str);
      // }

      const int lower_index = level_index;
      const int upper_index = get_nlte_vector_index(element,ion,upper);

      *gsl_matrix_ptr(rate_matrix_rad_bb, lower_index, lower_index) -= R;
      *gsl_matrix_ptr(rate_matrix_rad_bb, upper_index, lower_index) += R;
      *gsl_matrix_ptr(rate_matrix_coll_bb, lower_index, lower_index) -= C;
      *gsl_matrix_ptr(rate_matrix_coll_bb, upper_index, lower_index) += C;
      *gsl_matrix_ptr(rate_matrix_ntcoll_bb, lower_index, lower_index) -= NTC;
      *gsl_matrix_ptr(rate_matrix_ntcoll_bb, upper_index, lower_index) += NTC;
      if ((R < 0) || (C < 0))
        printout("  WARNING: Negative excitation rate from ion %d level %d to level %d\n", get_ionstage(element, ion), level, upper);
    }
  }
}


static void nltepop_matrix_add_ionisation(
  const int modelgridindex, const int element, const int ion,
  double *restrict s_renorm, gsl_matrix *rate_matrix_rad_bf, gsl_matrix *rate_matrix_coll_bf)
{
  assert(ion + 1 < get_nions(element)); // can't ionise the top ion
  const float T_e = get_Te(modelgridindex);
  const float nne = get_nne(modelgridindex);
  const int nionisinglevels = get_ionisinglevels(element, ion);
  const int maxrecombininglevel = get_maxrecombininglevel(element, ion + 1);

  for (int level = 0; level < nionisinglevels; level++)
  {
    const int level_index = get_nlte_vector_index(element, ion, level);

    // thermal collisional ionization, photoionisation and recombination processes
    const double epsilon_current = epsilon(element, ion, level);
    const int lower_index = level_index;

    for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element, ion, level); phixstargetindex++)
    {
      const int upper = get_phixsupperlevel(element, ion, level, phixstargetindex);
      const int upper_index = get_nlte_vector_index(element, ion + 1, upper);
      const double epsilon_trans = epsilon(element, ion + 1, upper) - epsilon_current;

      // ionization

      // the R part is slow!
      const double R_ionisation = get_corrphotoioncoeff(element, ion, level, phixstargetindex, modelgridindex);
      const double C_ionisation = col_ionization_ratecoeff(T_e, nne, element, ion, level, phixstargetindex, epsilon_trans);

      *gsl_matrix_ptr(rate_matrix_rad_bf, lower_index, lower_index) -= R_ionisation * s_renorm[level];
      *gsl_matrix_ptr(rate_matrix_rad_bf, upper_index, lower_index) += R_ionisation * s_renorm[level];
      *gsl_matrix_ptr(rate_matrix_coll_bf, lower_index, lower_index) -= C_ionisation * s_renorm[level];
      *gsl_matrix_ptr(rate_matrix_coll_bf, upper_index, lower_index) += C_ionisation * s_renorm[level];

      if ((R_ionisation < 0) || (C_ionisation < 0))
        printout("  WARNING: Negative ionization rate from ion_stage %d level %d phixstargetindex %d\n",
                 get_ionstage(element, ion), level, phixstargetindex);

      // recombination
      if (upper <= maxrecombininglevel) // we can skip this part if the functions below will return zero anyway
      {
        const double R_recomb = rad_recombination_ratecoeff(T_e, nne, element, ion + 1, upper, level, modelgridindex);
        const double C_recomb = col_recombination_ratecoeff(modelgridindex, element, ion + 1, upper, level, epsilon_trans);

        *gsl_matrix_ptr(rate_matrix_rad_bf, upper_index, upper_index) -= R_recomb * s_renorm[upper];
        *gsl_matrix_ptr(rate_matrix_rad_bf, lower_index, upper_index) += R_recomb * s_renorm[upper];
        *gsl_matrix_ptr(rate_matrix_coll_bf, upper_index, upper_index) -= C_recomb * s_renorm[upper];
        *gsl_matrix_ptr(rate_matrix_coll_bf, lower_index, upper_index) += C_recomb * s_renorm[upper];

        if ((R_recomb < 0) || (C_recomb < 0))
          printout("  WARNING: Negative recombination rate to ion_stage %d level %d phixstargetindex %d\n",
                   get_ionstage(element, ion), level, phixstargetindex);
     }
    }
  }
}


static void nltepop_matrix_add_nt_ionisation(
  const int modelgridindex, const int element, const int ion,
  const double *restrict s_renorm, gsl_matrix *restrict rate_matrix_ntcoll_bf)
{
  // collisional ionization by non-thermal electrons

  assert(ion + 1 < get_nions(element)); // can't ionise the top ion
  const double Y_nt = nt_ionization_ratecoeff(modelgridindex, element, ion);
  if (Y_nt < 0.)
  {
    printout("  WARNING: Negative NT_ionization rate from ion_stage %d\n", get_ionstage(element, ion));
  }

  const int nlevels = get_nlevels(element, ion);

  for (int upperion = ion + 1; upperion <= nt_ionisation_maxupperion(element, ion); upperion++)
  {
    const double Y_nt_thisupperion = Y_nt * nt_ionization_upperion_probability(modelgridindex, element, ion, upperion, false);

    if (Y_nt_thisupperion > 0.)
    {
      const int upper_groundstate_index = get_nlte_vector_index(element, upperion, 0);
      for (int level = 0; level < nlevels; level++)
      {
        const int lower_index = get_nlte_vector_index(element, ion, level);

        *gsl_matrix_ptr(rate_matrix_ntcoll_bf, lower_index, lower_index) -= Y_nt_thisupperion * s_renorm[level];
        *gsl_matrix_ptr(rate_matrix_ntcoll_bf, upper_groundstate_index, lower_index) += Y_nt_thisupperion * s_renorm[level];
      }
    }
  }
}


static void nltepop_matrix_normalise(
  const int modelgridindex, const int element,
  gsl_matrix *restrict rate_matrix,
  gsl_vector *restrict pop_norm_factor_vec)
{
  const unsigned int nlte_dimension = pop_norm_factor_vec->size;
  assert(pop_norm_factor_vec->size == nlte_dimension);
  assert(rate_matrix->size1 == nlte_dimension);
  assert(rate_matrix->size2 == nlte_dimension);

  // TODO: consider replacing normalisation by LTE populations with
  // GSL's gsl_linalg_balance_matrix(gsl_matrix * A, gsl_vector * D) function instead
  for (unsigned int column = 0; column < nlte_dimension; column++)
  {
    int ion, level;
    get_ion_level_of_nlte_vector_index(column, element, &ion, &level);

    gsl_vector_set(pop_norm_factor_vec, column, calculate_levelpop_lte(modelgridindex, element, ion, level));

    if ((level != 0) && (is_nlte(element, ion, level) == false))
    {
      // level is a superlevel, so add populations of higher levels to the norm factor
      for (int dummylevel = level + 1; dummylevel < get_nlevels(element, ion); dummylevel++)
      {
        if (is_nlte(element, ion, dummylevel) == false)
        {
          *gsl_vector_ptr(pop_norm_factor_vec, column) += calculate_levelpop_lte(modelgridindex, element, ion, dummylevel);
        }
      }
      // NOTE: above calculation is not always equal to the sum of LTE populations
      // since calculate_levelpop_lte imposes MINPOP minimum
      // printout("superlevel norm factor index %d is %g, partfunc is %g, partfunc*levelpop(SL)/g(SL) %g\n",
      //          column, gsl_vector_get(pop_norm_factor_vec, column), superlevel_partfunc[ion],
      //          superlevel_partfunc[ion] * calculate_levelpop_lte(modelgridindex,element,ion,level) / stat_weight(element,ion,level));
    }

    // apply the normalisation factor to this column in the rate_matrix
    gsl_vector_view column_view = gsl_matrix_column(rate_matrix, column);
    gsl_vector_scale(&column_view.vector, gsl_vector_get(pop_norm_factor_vec, column));
  }
}


static void set_element_pops_lte(const int modelgridindex, const int element)
{
  nltepop_reset_element(modelgridindex, element);

  // const int nions = get_nions(element);
  // for (int ion = 0; ion < nions; ion++)
  //   modelgrid[modelgridindex].composition[element].partfunct[ion] = calculate_partfunct(element, ion, modelgridindex);
  //
  // const float nne = get_nne(modelgridindex);
  // const double nnelement = get_abundance(modelgridindex, element) / elements[element].mass * get_rho(modelgridindex);
  // for (int ion = 0; ion < nions; ion++)
  // {
  //   double nnion;
  //   if (ion == 0)
  //     nnion = nnelement * ionfract(element, ion, modelgridindex, nne);
  //   else
  //     nnion = MINPOP;
  //
  //   modelgrid[modelgridindex].composition[element].groundlevelpop[ion] = (
  //     nnion * stat_weight(element,ion,0) / modelgrid[modelgridindex].composition[element].partfunct[ion]);
  //
  //   assert(isfinite(modelgrid[modelgridindex].composition[element].groundlevelpop[ion]));
  // }
}


static bool lumatrix_is_singular(const gsl_matrix *LU, const int element)
{
  size_t n = LU->size1;
  bool is_singular = false;

  for (size_t i = 0; i < n; i++)
  {
    const double u = gsl_matrix_get(LU, i, i);
    if (u == 0)
    {
      int ion = -1;
      int level = -1;
      get_ion_level_of_nlte_vector_index(i, element, &ion, &level);
      printout("NLTE disconnected level: Z=%d ionstage %d level+1 %d\n", get_element(element), get_ionstage(element, ion), level + 1);
      is_singular = true;
    }
  }

 return is_singular;
}


static bool nltepop_matrix_solve(
  const int element,
  const gsl_matrix *restrict rate_matrix,
  const gsl_vector *restrict balance_vector,
  gsl_vector *restrict popvec,
  const gsl_vector *restrict pop_normfactor_vec)
// solve rate_matrix * x = balance_vector,
// then popvec[i] = x[i] / pop_norm_factor_vec[i]
{
  bool completed_solution;
  const unsigned int nlte_dimension = balance_vector->size;
  assert(pop_normfactor_vec->size == nlte_dimension);
  assert(rate_matrix->size1 == nlte_dimension);
  assert(rate_matrix->size2 == nlte_dimension);

  gsl_vector *x = gsl_vector_alloc(nlte_dimension); // population solution vector (normalised)

  // make a copy of the rate matrix for the LU decomp
  gsl_matrix *rate_matrix_LU_decomp = gsl_matrix_alloc(nlte_dimension, nlte_dimension);
  gsl_matrix_memcpy(rate_matrix_LU_decomp, rate_matrix);

  gsl_permutation *p = gsl_permutation_alloc(nlte_dimension);

  int s; // sign of the transformation
  gsl_linalg_LU_decomp(rate_matrix_LU_decomp, p, &s);

  if (lumatrix_is_singular(rate_matrix_LU_decomp, element))
  {
    printout("ERROR: NLTE matrix is singular for element Z=%d!\n", get_element(element));
    // abort();
    completed_solution = false;
  }
  else
  {
    gsl_error_handler_t *previous_handler = gsl_set_error_handler(gsl_error_handler_printout);

    // solve matrix equation: rate_matrix * x = balance_vector for x (population vector)
    gsl_linalg_LU_solve(rate_matrix_LU_decomp, p, balance_vector, x);

    gsl_set_error_handler(previous_handler);

    //gsl_linalg_HH_solve (&m.matrix, &b.vector, x);

    const double TOLERANCE = 1e-40;

    gsl_vector *gsl_work_vector = gsl_vector_alloc(nlte_dimension);
    double error_best = -1.;
    gsl_vector *x_best = gsl_vector_alloc(nlte_dimension); //population solution vector with lowest error
    gsl_vector *residual_vector = gsl_vector_alloc(nlte_dimension);
    int iteration;
    for (iteration = 0; iteration < 10; iteration++)
    {
      if (iteration > 0)
        gsl_linalg_LU_refine(rate_matrix, rate_matrix_LU_decomp, p, balance_vector, x, gsl_work_vector);

      gsl_vector_memcpy(residual_vector, balance_vector);
      gsl_blas_dgemv(CblasNoTrans, 1.0, rate_matrix, x, -1.0, residual_vector); // calculate Ax - b = residual
      const double error = fabs(gsl_vector_get(residual_vector, gsl_blas_idamax(residual_vector))); // value of the largest absolute residual

      if (error < error_best || error_best < 0.)
      {
        gsl_vector_memcpy(x_best, x);
        error_best = error;
      }
      //printout("Linear algebra solver iteration %d has a maximum residual of %g\n",iteration,error);
      if (error < TOLERANCE)
      {
        break;
      }
    }
    if (error_best >= 0.)
    {
      // printout("  NLTE solver matrix LU_refine: After %d iterations, keeping solution vector with a max residual of %g\n",iteration,error_best);
      if (error_best > 1e-10)
        printout("  NLTE solver matrix LU_refine: After %d iterations, best solution vector has a max residual of %g (WARNING!)\n",iteration,error_best);
      gsl_vector_memcpy(x, x_best);
    }
    gsl_vector_free(x_best);
    gsl_vector_free(gsl_work_vector);

    // get the real populations using the x vector and the normalisation factors
    gsl_vector_memcpy(popvec, x);
    gsl_vector_mul(popvec, pop_normfactor_vec);
    // popvec will be used contains the real population densities

    for (unsigned int row = 0; row < nlte_dimension; row++)
    {
      double recovered_balance_vector_elem = 0.;
      gsl_vector_const_view row_view = gsl_matrix_const_row(rate_matrix, row);
      gsl_blas_ddot(&row_view.vector, x, &recovered_balance_vector_elem);

      int ion, level;
      get_ion_level_of_nlte_vector_index(row, element, &ion, &level);

      // printout("index %4d (ion_stage %d level%4d): residual %+.2e recovered balance: %+.2e normed pop %.2e pop %.2e departure ratio %.4f\n",
      //          row,get_ionstage(element,ion),level, gsl_vector_get(residual_vector,row),
      //          recovered_balance_vector_elem, gsl_vector_get(x,row),
      //          gsl_vector_get(popvec, row),
      //          gsl_vector_get(x, row) / gsl_vector_get(x,get_nlte_vector_index(element,ion,0)));

      if (gsl_vector_get(popvec, row) < 0.0)
      {
        printout("  WARNING: NLTE solver gave negative population to index %d (Z=%d ion_stage %d level %d), pop = %g. Replacing with LTE pop of %g\n",
                 row, get_element(element), get_ionstage(element, ion), level, gsl_vector_get(x, row) * gsl_vector_get(pop_normfactor_vec, row), gsl_vector_get(pop_normfactor_vec, row));
        gsl_vector_set(popvec, row, gsl_vector_get(pop_normfactor_vec, row));
      }
    }

    gsl_vector_free(residual_vector);
    completed_solution = true;
  }

  gsl_vector_free(x);
  gsl_matrix_free(rate_matrix_LU_decomp);
  gsl_permutation_free(p);

  return completed_solution;
}


void solve_nlte_pops_element(const int element, const int modelgridindex, const int timestep, const int nlte_iter)
// solves the statistical balance equations to find NLTE level populations for all ions of an element
// (ionisation balance follows from this too)
{
  const int atomic_number = get_element(element);

  if (get_abundance(modelgridindex, element) <= 0.)
  {
    //abundance of this element is zero, so do not store any NLTE populations
    printout("Not solving for NLTE populations in cell %d at timestep %d for element Z=%d due to zero abundance\n",
             modelgridindex, timestep, atomic_number);

    nltepop_reset_element(modelgridindex, element);
    return;
  }

  // can save memory by using a combined rate matrix at the cost of diagnostic information
  const bool individual_process_matricies = true;

  const double t_mid = time_step[timestep].mid;
  const int nions = get_nions(element);

  printout("Solving for NLTE populations in cell %d at timestep %d iteration %d for element Z=%d (mass fraction %.2e, population %.2e)\n",
           modelgridindex, timestep, nlte_iter, atomic_number, get_abundance(modelgridindex, element),
           get_abundance(modelgridindex, element) / elements[element].mass * get_rho(modelgridindex));

  // LTE test, make sure binned radfield is off
  //set_TR(modelgridindex,3000);
  //set_W(modelgridindex,1.0);
  //printout("T_E %g T_R was %g, setting to 3000 \n",get_Te(modelgridindex),get_TR(modelgridindex));

  double superlevel_partfunc[nions]; // space is allocated for every ion, even if it does not have a superlevel
  const int nlte_dimension = get_element_nlte_dimension_and_slpartfunc(modelgridindex, element, nions, superlevel_partfunc);

  // printout("NLTE: the vector dimension is %d", nlte_dimension);

  gsl_matrix *rate_matrix = gsl_matrix_calloc(nlte_dimension, nlte_dimension);
  gsl_matrix *rate_matrix_rad_bb;
  gsl_matrix *rate_matrix_coll_bb;
  gsl_matrix *rate_matrix_ntcoll_bb;
  gsl_matrix *rate_matrix_rad_bf;
  gsl_matrix *rate_matrix_coll_bf;
  gsl_matrix *rate_matrix_ntcoll_bf;
  if (individual_process_matricies)
  {
    rate_matrix_rad_bb = gsl_matrix_calloc(nlte_dimension, nlte_dimension);
    rate_matrix_coll_bb = gsl_matrix_calloc(nlte_dimension, nlte_dimension);
    rate_matrix_ntcoll_bb = gsl_matrix_calloc(nlte_dimension, nlte_dimension);
    rate_matrix_rad_bf = gsl_matrix_calloc(nlte_dimension, nlte_dimension);
    rate_matrix_coll_bf = gsl_matrix_calloc(nlte_dimension, nlte_dimension);
    rate_matrix_ntcoll_bf = gsl_matrix_calloc(nlte_dimension, nlte_dimension);
  }
  else
  {
    // alias the single matrix accounting for all processes
    rate_matrix_rad_bb = rate_matrix;
    rate_matrix_coll_bb = rate_matrix;
    rate_matrix_ntcoll_bb = rate_matrix;
    rate_matrix_rad_bf = rate_matrix;
    rate_matrix_coll_bf = rate_matrix;
    rate_matrix_ntcoll_bf = rate_matrix;
  }

  gsl_vector *const balance_vector = gsl_vector_calloc(nlte_dimension);

  if (!balance_vector)
  {
    printout("Cannot allocate NLTE rate matrix/balance vector memory.\n");
    abort();
  }

  // printout("  Adding rates for ion stages:");
  for (int ion = 0; ion < nions; ion++)
  {
    // const int ionstage = get_ionstage(element, ion);
    // printout(" %d", ionstage);

    const int nlevels = get_nlevels(element, ion);
    const int nlevels_nlte = get_nlevels_nlte(element, ion);
    double s_renorm[nlevels];
    for (int level = 0; level <= nlevels_nlte; level++)
    {
      s_renorm[level] = 1.0;
    }

    for (int level = (nlevels_nlte + 1); level < nlevels; level++)
    {
      s_renorm[level] = superlevel_boltzmann(modelgridindex, element, ion, level) / superlevel_partfunc[ion];
    }

    nltepop_matrix_add_boundbound(
      modelgridindex, element, ion, t_mid, s_renorm, rate_matrix_rad_bb, rate_matrix_coll_bb, rate_matrix_ntcoll_bb);

    if (ion < nions - 1)
    {
      // this is the slowest component
      nltepop_matrix_add_ionisation(modelgridindex, element, ion, s_renorm, rate_matrix_rad_bf, rate_matrix_coll_bf);
      if (NT_ON)
        nltepop_matrix_add_nt_ionisation(modelgridindex, element, ion, s_renorm, rate_matrix_ntcoll_bf);
    }
  }
  // printout("\n");

  if (individual_process_matricies)
  {
    // sum the matricies for each transition process to get a total rate matrix
    gsl_matrix_add(rate_matrix, rate_matrix_rad_bb);
    gsl_matrix_add(rate_matrix, rate_matrix_coll_bb);
    gsl_matrix_add(rate_matrix, rate_matrix_ntcoll_bb);
    gsl_matrix_add(rate_matrix, rate_matrix_rad_bf);
    gsl_matrix_add(rate_matrix, rate_matrix_coll_bf);
    gsl_matrix_add(rate_matrix, rate_matrix_ntcoll_bf);
  }

  // replace the first row of the matrix and balance vector with the normalisation
  // constraint on the total element population
  gsl_vector_view first_row_view = gsl_matrix_row(rate_matrix, 0);
  gsl_vector_set_all(&first_row_view.vector, 1.0);
  // set first balance vector entry to the element population (all other entries will be zero)
  const double element_population = get_abundance(modelgridindex, element) / elements[element].mass * get_rho(modelgridindex);
  gsl_vector_set(balance_vector, 0, element_population);

  // calculate the normalisation factors and apply them to the matrix
  // columns and balance vector elements
  gsl_vector *pop_norm_factor_vec = gsl_vector_calloc(nlte_dimension);
  nltepop_matrix_normalise(modelgridindex, element, rate_matrix, pop_norm_factor_vec);

  // printout("Rate matrix | balance vector:\n");
  // for (int row = 0; row < nlte_dimension; row++)
  // {
  //   for (int column = 0; column < nlte_dimension; column++)
  //   {
  //     char str[15];
  //     sprintf(str, "%+.1e ", gsl_matrix_get(rate_matrix, row, column));
  //     printout(str);
  //   }
  //   printout("| ");
  //   char str[15];
  //   sprintf(str, "%+.1e\n", gsl_vector_get(balance_vector, row));
  //   printout(str);
  // }
  // printout("\n");

  // eliminate barely-interacting levels from the NLTE matrix by removing
  // their interactions and setting their normalised populations (probably departure coeff) to 1.0
  // filter_nlte_matrix(element, rate_matrix, balance_vector, pop_norm_factor_vec);

  gsl_vector *popvec = gsl_vector_alloc(nlte_dimension); // the true population densities

  const bool matrix_solve_success = nltepop_matrix_solve(element, rate_matrix, balance_vector, popvec, pop_norm_factor_vec);

  if (!matrix_solve_success)
  {
    printout("WARNING: Can't solve for NLTE populations in cell %d at timestep %d for element Z=%d due to singular matrix. Attempting to use LTE solution instead\n",
             modelgridindex, timestep, atomic_number);
    set_element_pops_lte(modelgridindex, element);
  }
  else
  {
    // double ion_populations[nions];
    for (int ion = 0; ion < nions; ion++)
    {
      const int nlevels_nlte = get_nlevels_nlte(element, ion);
      const int index_gs = get_nlte_vector_index(element, ion, 0);
      // const int ion_stage = get_ionstage(element, ion);
      // printout("  [ion_stage %d]\n", ion_stage);
      //
      // printout("    For ion_stage %d, the ground state populations are %g (function) and %g (matrix result with normed pop %g, ltepopnormfactor %g)\n",get_ionstage(element,ion),
      //          get_groundlevelpop(modelgridindex, element, ion), gsl_vector_get(popvec, index_gs),
      //          gsl_vector_get(x, index_gs), gsl_vector_get(pop_norm_factor_vec, index_gs));

      // store the NLTE level populations
      const int nlte_start = elements[element].ions[ion].first_nlte;
      double solution_ion_pop = 0.0;
      for (int level = 1; level <= nlevels_nlte; level++)
      {
        const int index = get_nlte_vector_index(element, ion, level);
        modelgrid[modelgridindex].nlte_pops[nlte_start + level - 1] = gsl_vector_get(popvec, index) / get_rho(modelgridindex);
        solution_ion_pop += gsl_vector_get(popvec, index);
      }

      // store the superlevel population if there is one
      if (ion_has_superlevel(element, ion)) //a superlevel exists
      {
        const int index_sl = get_nlte_vector_index(element, ion, nlevels_nlte + 1);
        modelgrid[modelgridindex].nlte_pops[nlte_start + nlevels_nlte] = (
          gsl_vector_get(popvec, index_sl) / modelgrid[modelgridindex].rho / superlevel_partfunc[ion]);

        // printout("solve_nlte_pops_element: The Z=%d ionstage %d superlevel population is %g with rho %g and superlevel_partfunc %g Te %g scaled pop stored as %g\n", get_element(element), get_ionstage(element, ion), gsl_vector_get(popvec, index_sl), modelgrid[modelgridindex].rho, superlevel_partfunc[ion], get_Te(modelgridindex), modelgrid[modelgridindex].nlte_pops[nlte_start + nlevels_nlte]);
        // the stored population is already divided by the partfunc, so just multiply
        // it by the superlevel_boltzmann to get the population of a level in the SL

        solution_ion_pop += gsl_vector_get(popvec, index_sl);
      }
      // printout("    I had a ground level pop of %g, a part fn of %g and therefore an ion pop of %g\n",
      //          modelgrid[modelgridindex].composition[element].groundlevelpop[ion],
      //          modelgrid[modelgridindex].composition[element].partfunct[ion],
      //          (modelgrid[modelgridindex].composition[element].partfunct[ion] *
      //            modelgrid[modelgridindex].composition[element].groundlevelpop[ion]
      //            / stat_weight(element, ion, 0)));

      // ionstagepop here must be called before setting the new ground level population
      // printout("    For ion_stage %d the total population is %g, but was previously %g\n",
      //          ion_stage,solution_ion_pop,ionstagepop(modelgridindex, element, ion));

      // store the ground level population
      modelgrid[modelgridindex].composition[element].groundlevelpop[ion] = gsl_vector_get(popvec, index_gs);
      solution_ion_pop += gsl_vector_get(popvec, index_gs);

      precalculate_partfuncts(modelgridindex);

      // ion_populations[ion] = solution_ion_pop;
      // if (ion > 0)
      // {
      //   const double gspopratio = modelgrid[modelgridindex].composition[element].groundlevelpop[ion-1] / modelgrid[modelgridindex].composition[element].groundlevelpop[ion];
      //
      //   const double ionpot = epsilon(element,ion,0) - epsilon(element,ion-1,0);
      //   const float T_e = get_Te(modelgridindex);
      //   const double partfunct_ratio = modelgrid[modelgridindex].composition[element].partfunct[ion-1] / modelgrid[modelgridindex].composition[element].partfunct[ion];
      //   const double gs_g_ratio = stat_weight(element,ion-1,0) / stat_weight(element,ion,0);
      //   const double sbphi_gs = gs_g_ratio * SAHACONST * pow(T_e,-1.5) * exp(ionpot/KB/T_e) * get_nne(modelgridindex);
      //   const double solution_ion_pop_ratio = ion_populations[ion-1] / ion_populations[ion];
      //   const double sbphi = partfunct_ratio * SAHACONST * pow(T_e,-1.5) * exp(ionpot/KB/T_e) * get_nne(modelgridindex);
      //
      //   printout("    The ratio of groundlevel pops (ion %d)/(ion %d) is %g, Saha-Boltzmann value is %g ratio %g\n",
      //            get_ionstage(element,ion-1),ion_stage,gspopratio,sbphi_gs,gspopratio/sbphi_gs);
      //   printout("    The ratio of total pops (ion %d)/(ion %d) is %g, Saha-Boltzmann value is %g ratio %g\n",
      //            get_ionstage(element,ion-1),ion_stage,solution_ion_pop_ratio,sbphi,solution_ion_pop_ratio/sbphi);
      //   const float nne = get_nne(modelgridindex);
      //   // calculate_partfunct(element, ion, modelgridindex) * modelgrid[modelgridindex].composition[element].groundlevelpop[ion] / stat_weight(element,ion,0))
      //   printout("    The corresponding phi-factor is: %g\n", (solution_ion_pop_ratio / nne));
      //   printout("    The ionization solver gives phi(ion_stage=%d / %d) = %g\n",
      //            get_ionstage(element, ion-1), get_ionstage(element, ion), phi(element, ion-1, modelgridindex));
      // }


      //double nne = get_nne(modelgridindex);
      //printout("  From ion fract, the ion pop should be %g\n", ionfract(element, ion, modelgridindex, nne)*get_abundance(modelgridindex,element)/ elements[element].mass * get_rho(modelgridindex));
      //printout("  I think that the element population is: %g (from abundance %g and rho %g)\n", get_abundance(modelgridindex,element)/elements[element].mass*get_rho(modelgridindex), get_abundance(modelgridindex,element), get_rho(modelgridindex));
      //printout("  I currently think that the top ion is: %d\n", elements_uppermost_ion[tid][element]);
    }

    const double elem_pop_abundance = get_abundance(modelgridindex, element) / elements[element].mass * get_rho(modelgridindex);
    const double elem_pop_matrix = gsl_blas_dasum(popvec);
    const double elem_pop_error_percent = fabs((elem_pop_abundance / elem_pop_matrix) - 1) * 100;
    if (elem_pop_error_percent > 1.0)
    {
      printout("  WARNING: The Z=%d element population is: %g (from abundance) and %g (from matrix solution sum of level pops), error: %.1f%%. Forcing element pops to LTE.\n",
               atomic_number, elem_pop_abundance, elem_pop_matrix, elem_pop_error_percent);
      set_element_pops_lte(modelgridindex, element);
    }

    if (individual_process_matricies)
    {
      print_element_rates_summary(element, popvec, rate_matrix_rad_bb, rate_matrix_coll_bb, rate_matrix_ntcoll_bb, rate_matrix_rad_bf, rate_matrix_coll_bf, rate_matrix_ntcoll_bf);
    }

    // if (individual_process_matricies && (atomic_number == 26 && timestep % 2 == 0))
    // {
    //   const int ionstage = 2;
    //   const int ion = ionstage - get_ionstage(element, 0);
    //
    //   print_level_rates(element, ion, 0, popvec, rate_matrix_rad_bb, rate_matrix_coll_bb, rate_matrix_ntcoll_bb, rate_matrix_rad_bf, rate_matrix_coll_bf, rate_matrix_ntcoll_bf);
    //   print_level_rates(element, ion, 1, popvec, rate_matrix_rad_bb, rate_matrix_coll_bb, rate_matrix_ntcoll_bb, rate_matrix_rad_bf, rate_matrix_coll_bf, rate_matrix_ntcoll_bf);
    //   print_level_rates(element, ion, 28, popvec, rate_matrix_rad_bb, rate_matrix_coll_bb, rate_matrix_ntcoll_bb, rate_matrix_rad_bf, rate_matrix_coll_bf, rate_matrix_ntcoll_bf);
    //   print_level_rates(element, ion, 29, popvec, rate_matrix_rad_bb, rate_matrix_coll_bb, rate_matrix_ntcoll_bb, rate_matrix_rad_bf, rate_matrix_coll_bf, rate_matrix_ntcoll_bf);
    //   print_level_rates(element, ion, 30, popvec, rate_matrix_rad_bb, rate_matrix_coll_bb, rate_matrix_ntcoll_bb, rate_matrix_rad_bf, rate_matrix_coll_bf, rate_matrix_ntcoll_bf);
    //   print_level_rates(element, ion, 36, popvec, rate_matrix_rad_bb, rate_matrix_coll_bb, rate_matrix_ntcoll_bb, rate_matrix_rad_bf, rate_matrix_coll_bf, rate_matrix_ntcoll_bf);
    //   print_level_rates(element, ion, 37, popvec, rate_matrix_rad_bb, rate_matrix_coll_bb, rate_matrix_ntcoll_bb, rate_matrix_rad_bf, rate_matrix_coll_bf, rate_matrix_ntcoll_bf);
    //   print_level_rates(element, ion, 38, popvec, rate_matrix_rad_bb, rate_matrix_coll_bb, rate_matrix_ntcoll_bb, rate_matrix_rad_bf, rate_matrix_coll_bf, rate_matrix_ntcoll_bf);
    //   print_level_rates(element, ion, 39, popvec, rate_matrix_rad_bb, rate_matrix_coll_bb, rate_matrix_ntcoll_bb, rate_matrix_rad_bf, rate_matrix_coll_bf, rate_matrix_ntcoll_bf);
    //   print_level_rates(element, ion, 40, popvec, rate_matrix_rad_bb, rate_matrix_coll_bb, rate_matrix_ntcoll_bb, rate_matrix_rad_bf, rate_matrix_coll_bf, rate_matrix_ntcoll_bf);
    //   print_level_rates(element, ion, 41, popvec, rate_matrix_rad_bb, rate_matrix_coll_bb, rate_matrix_ntcoll_bb, rate_matrix_rad_bf, rate_matrix_coll_bf, rate_matrix_ntcoll_bf);
    //   const int slindex = get_nlevels_nlte(element, ion) + 1;
    //   print_level_rates(element, ion, slindex, popvec, rate_matrix_rad_bb, rate_matrix_coll_bb, rate_matrix_ntcoll_bb, rate_matrix_rad_bf, rate_matrix_coll_bf, rate_matrix_ntcoll_bf);
    // }
  }

  if (individual_process_matricies)
  {
    gsl_matrix_free(rate_matrix_rad_bb);
    gsl_matrix_free(rate_matrix_coll_bb);
    gsl_matrix_free(rate_matrix_ntcoll_bb);
    gsl_matrix_free(rate_matrix_rad_bf);
    gsl_matrix_free(rate_matrix_coll_bf);
    gsl_matrix_free(rate_matrix_ntcoll_bf);
  }

  gsl_vector_free(popvec);

  gsl_matrix_free(rate_matrix);
  gsl_vector_free(balance_vector);
  gsl_vector_free(pop_norm_factor_vec);
}


// this does single ion solving and will be deprecated at some point
double solve_nlte_pops_ion(int element, int ion, int modelgridindex, int timestep)
//solves for nlte correction factors to level populations for levels
{
  if (get_nlevels(element,ion) > 1)
  {
    int lower, level_use, lower_use, upper_use;

    double *rate_matrix;
    double *balance_vector;

    int ionisinglevels, ndowntrans, nuptrans, i;
    double epsilon_current;
    double R, C;
    double t_mid;
    double epsilon_trans;
    int lineindex;

    int nlte_size, nlte_start;

    double test_ratio;//,lag, check;
    int super_level;

    const float T_e = get_Te(modelgridindex);
    const float nne = get_nne(modelgridindex);

    printout("Solving for NLTE populations in cell %d. Doing element %d, ion %d. I think it's timestep %d\n", modelgridindex, element, ion, timestep);
    //printout("Current Te %g and ne %g\n",T_e, nne);

    int nlevels_nlte = get_nlevels_nlte(element,ion);
    t_mid = time_step[timestep].mid;

    //PKT dummy;
    //dummy.where = modelgridindex;

    if (nlevels_nlte == (get_nlevels(element, ion) - 1))
    {
      nlte_size = nlevels_nlte + 2;
      super_level = 0;
    }
    else
    {
      nlte_size = nlevels_nlte + 3;
      super_level = 1;
    }
    //that's the total number of nlte_levels (i.e. the size of the
    //storage). Our rate matrix will need to be of this dimension +3: the
    //ground state, the "super level" and the ground state of the ion
    //above. If there's no super level needed then we only need +2
    nlte_start = elements[element].ions[ion].first_nlte;

    //      if (get_groundlevelpop(modelgridindex,element,ion) > (1.2*MINPOP))
    if (get_abundance(modelgridindex, element) > 0.0)
    {
      if ((rate_matrix = calloc((nlte_size)*(nlte_size), sizeof(double))) == NULL)
      {
        printout("Cannot allocate NLTE rate matrix memory.\n");
      }

      if ((balance_vector = calloc((nlte_size), sizeof(double))) == NULL)
      {
        printout("Cannot allocate NLTE vector memory.\n");
      }

      //printf("rate %p balance %p NULL %p\n", rate_matrix, balance_vector, NULL);
      //printout("I think there are %d levels to deal with and managed to allocate memory.\n", nlte_size);

      double superlevel_partition = 0.0;
      for (int level = 1; level < get_nlevels(element,ion); level++)
      {
        if (is_nlte(element, ion, level) == 0)
        {
          superlevel_partition += superlevel_boltzmann(modelgridindex,element,ion,level);
        }
      }

      double upperion_partition = 0.0;
      if (ion < get_nions(element))
      {
        for (int level = 0; level < get_nlevels(element,ion+1); level++)
        {
          upperion_partition += calculate_exclevelpop(modelgridindex,element,ion+1,level);
        }
      }


      if (superlevel_partition > 0.0)
      {
        printout("I found a superlevel and have computed a partition function for its substates of %g.\n", superlevel_partition);
      }
      else
      {
        printout("I don't know about any super level for this case.\n");
      }


      for (int level = 0; level < get_nlevels(element,ion); level++)
      {
        epsilon_current = epsilon(element,ion,level);
        ndowntrans = get_ndowntrans(element, ion, level);
        nuptrans = get_nuptrans(element, ion, level);

        double s_renorm;

        // deexcitation
        for (i = 0; i < ndowntrans; i++)
        {
          epsilon_trans = epsilon_current - epsilon(element, ion, lower);
          lineindex = elements[element].ions[ion].levels[level].downtrans_lineindicies[i];
          lower = linelist[lineindex].lowerlevelindex;

          mastate[tid].element = element;
          mastate[tid].ion = ion;
          mastate[tid].level = level;

          R = rad_deexcitation_ratecoeff(modelgridindex,element,ion,level,lower,epsilon_trans,lineindex,t_mid);
          C = col_deexcitation_ratecoeff(T_e, nne, epsilon_trans, lineindex);

          s_renorm = 1.0;

          if ((level == 0) || (is_nlte(element, ion, level) == 1))
          {
            level_use = level;
          }
          else
          {
            level_use = nlevels_nlte + 1;
            s_renorm = superlevel_boltzmann(modelgridindex,element,ion,level) / superlevel_partition;
          }

          if ((lower == 0) || (is_nlte(element, ion, lower) == 1))
          {
            lower_use = lower;
          }
          else
          {
            lower_use = nlevels_nlte + 1;
          }

          rate_matrix[level_use*(nlte_size) + level_use] -= (R + C) * s_renorm;
          rate_matrix[lower_use*(nlte_size) + level_use] += (R + C) * s_renorm;
          //printout("First using %d of %d\n", level_use*(nlte_size) + level_use, nlte_size*nlte_size);
          //printout("Second using %d of %d\n", lower_use*(nlte_size) + level_use, nlte_size*nlte_size);
        }

        // excitation
        for (i = 0; i < nuptrans; i++)
        {
          lineindex = elements[element].ions[ion].levels[level].uptrans_lineindicies[i];
          const int upper = linelist[lineindex].upperlevelindex;
          epsilon_trans = epsilon(ion, level, upper) - epsilon_current;

          mastate[tid].element = element;
          mastate[tid].ion = ion;
          mastate[tid].level = level;

          R = rad_excitation_ratecoeff(modelgridindex, element, ion, level, upper, epsilon_trans, lineindex, t_mid);
          C = col_excitation_ratecoeff(T_e, nne, lineindex, epsilon_trans);

          if ((level == 0) || (is_nlte(element, ion, level) == 1))
          {
            level_use = level;
            s_renorm = 1.0;
          }
          else
          {
            level_use = nlevels_nlte+1;
            s_renorm = superlevel_boltzmann(modelgridindex,element,ion,level) / superlevel_partition;
          }
          if ((upper == 0) || (is_nlte(element, ion, upper) == 1))
          {
            upper_use = upper;
          }
          else
          {
            upper_use = nlevels_nlte+1;
          }
          rate_matrix[level_use*(nlte_size) + level_use] -= (R + C) * s_renorm;
          rate_matrix[upper_use*(nlte_size) + level_use] += (R + C) * s_renorm;
          //printout("Third using %d of %d\n", level_use*(nlte_size) + level_use, nlte_size*nlte_size);
          //printout("Fourth using %d of %d\n", upper_use*(nlte_size) + level_use, nlte_size*nlte_size);
        }

        if (NT_ON && ion < get_nions(element)-1)
        {
          double Y = nt_ionization_ratecoeff(modelgridindex,element,ion);

          if ((level == 0) || (is_nlte(element, ion, level) == 1))
          {
            level_use = level;
            s_renorm = 1.0;
          }
          else
          {
            level_use = nlevels_nlte + 1; //the super level
            s_renorm = superlevel_boltzmann(modelgridindex,element,ion,level) / superlevel_partition;
          }

          upper_use = nlte_size - 1; //the continuum
          rate_matrix[level_use*(nlte_size) + level_use] -= Y * s_renorm;
          rate_matrix[upper_use*(nlte_size) + level_use] += Y * s_renorm;
        }

        //now put in the photoionization/recombination processes
        ionisinglevels = get_ionisinglevels(element,ion);
        if (ion < get_nions(element)-1 && level < ionisinglevels)  //&& get_ionstage(element,ion) < get_element(element)+1)
        {
          s_renorm = 1.0;

          if ((level == 0) || (is_nlte(element, ion, level) == 1))
          {
            level_use = level;
          }
          else
          {
            level_use = nlevels_nlte+1; //the super level
            s_renorm = superlevel_boltzmann(modelgridindex,element,ion,level) / superlevel_partition;
          }

          // ionization
          mastate[tid].element = element;
          mastate[tid].ion = ion;
          mastate[tid].level = level;
          for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion,level); phixstargetindex++)
          {
            int upper = get_phixsupperlevel(element,ion,level,phixstargetindex);
            epsilon_trans = epsilon(element,ion+1,upper) - epsilon_current;
            R = get_corrphotoioncoeff(element,ion,level,phixstargetindex,modelgridindex);
            C = col_ionization_ratecoeff(T_e, nne, element, ion, level, phixstargetindex, epsilon_trans);

            upper_use = nlte_size - 1; //ion above

            rate_matrix[level_use*(nlte_size) + level_use] -= (R + C) * s_renorm;
            rate_matrix[upper_use*(nlte_size) + level_use] += (R + C) * s_renorm;
          }

          // recombination
          mastate[tid].element = element;
          mastate[tid].ion = ion+1;
          for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion,level); phixstargetindex++)
          {
            int upper = get_phixsupperlevel(element,ion,level,phixstargetindex);
            mastate[tid].level = upper;
            epsilon_trans = epsilon(element,ion+1,upper) - epsilon_current;
            R = rad_recombination_ratecoeff(T_e, nne, element, ion+1, upper, level, modelgridindex);
            //printout("rad recombination of element %d, ion %d, level %d, to lower level %d has rate %g (ne %g and Te %g)\n",element,ion,mastate[tid].level,level,R/nne,nne,T_e);
            //printout("%d %d %d %d %g %g %g \n",element,ion,mastate[tid].level,level,R/nne,nne,T_e);
            C = col_recombination_ratecoeff(modelgridindex, element, ion + 1, upper, level, epsilon_trans);
            //C=C*1.e-10;

            double upper_renorm = calculate_exclevelpop(modelgridindex,element,ion+1,upper) / upperion_partition;
            //TODO: would the line below be correct, or is it equal to the above line?
            //double upper_renorm = superlevel_boltzmann(modelgridindex,element,ion+1,upper) / upperion_partition;
            rate_matrix[upper_use*(nlte_size) + upper_use] -= (R + C) * upper_renorm;
            rate_matrix[level_use*(nlte_size) + upper_use] += (R + C) * upper_renorm;
          }

          //balance_vector[level_use] += -1. * get_groundlevelpop(modelgridindex,element,ion+1) * (R + C);
        }

      }

      //replace the first line of the matrix with the normalisation constraint
      rate_matrix[0] = 1.0;
      balance_vector[0] = ionstagepop(modelgridindex,element,ion);
      for (int level = 1; level < nlte_size; level++)
      {
        rate_matrix[level] = 1.0;
        balance_vector[level] = 0.0;
      }
      rate_matrix[nlte_size-1] = 0.0;

      for (int level = 0; level < nlte_size; level++)
      {
        for (level_use = 0; level_use < nlte_size; level_use++)
        {
          //printout("%g ", rate_matrix[level*nlte_size + level_use]);
          if (!isfinite(rate_matrix[level*nlte_size + level_use]))
          {
            printout("[fatal]: NLTE matrix with non-finite element: %d %d %g\n", level, level_use, rate_matrix[level*nlte_size + level_use]);
            printout("[fatal]: found when handling element %d and ion %d\n", element, ion);
            printout("[fatal]: the relevant ground state populations are %g and %g\n",get_groundlevelpop(modelgridindex,element,ion),get_groundlevelpop(modelgridindex,element,ion+1));
            abort();
          }
        }
        // printout("\n");
      }
      //printout("\n");

      for (int level = 0; level < nlte_size; level++)
      {
        //	  printout("%g ",balance_vector[level] );
        if (!isfinite(balance_vector[level]))
        {
          printout("[fatal]: NLTE balance with non-finite element: %d %g\n", level, balance_vector[level]);
          printout("[fatal]: found when handling element %d and ion %d\n", element, ion);
          printout("[fatal]: the relevant ground state populations are %g and %g\n",get_groundlevelpop(modelgridindex,element,ion), get_groundlevelpop(modelgridindex,element,ion+1));
          abort();
        }
      }
      // printout("\n");
      // printout("\n");

      //printf("rate %p balance %p NULL %p\n", rate_matrix, balance_vector, NULL);


      //next bunch of lines are based on a gsl example

      gsl_matrix_view m = gsl_matrix_view_array(rate_matrix, nlte_size, nlte_size);
      gsl_permutation *p = gsl_permutation_alloc(nlte_size);

      int s; //sign of the transformation
      gsl_linalg_LU_decomp(&m.matrix, p, &s);

      gsl_vector_view b = gsl_vector_view_array(balance_vector, nlte_size);
      gsl_vector *x = gsl_vector_alloc(nlte_size);
      gsl_linalg_LU_solve(&m.matrix, p, &b.vector, x); //solve matrix equation m * x = b for x (populations)
      //gsl_linalg_HH_solve (&m.matrix, &b.vector, x);

      //printf("after solving: rate %p balance %p NULL %p\n", rate_matrix, balance_vector, NULL);

      //printout("I did the GSL stuff - didnt' blow up at least!\n");

      printout("The ground state populations were %g, %g, %g and %g\n",
        get_groundlevelpop(modelgridindex,element,ion),gsl_vector_get(x,0),
        get_groundlevelpop(modelgridindex,element,ion+1),gsl_vector_get(x,nlte_size-1));

      //printout("The partition functions (and ratios to gs) were: %g (%g) %g (%g)\n", modelgrid[modelgridindex].composition[element].partfunct[ion], modelgrid[modelgridindex].composition[element].partfunct[ion]/stat_weight(element,ion,0),modelgrid[modelgridindex].composition[element].partfunct[ion+1], modelgrid[modelgridindex].composition[element].partfunct[ion+1]/stat_weight(element,ion+1,0));


      if ((get_groundlevelpop(modelgridindex,element,ion) > (1.1 * MINPOP)) && (gsl_vector_get(x,0) > (1.1 * MINPOP)))
      {
        test_ratio = get_groundlevelpop(modelgridindex,element,ion)/gsl_vector_get(x,0);
        if (test_ratio < 1)
        {
          test_ratio = 1. / test_ratio;
        }
      }
      else
      {
        test_ratio = 0.0;
      }

      double test_ratio_upper;
      if ((get_groundlevelpop(modelgridindex,element,ion+1) > (1.1 * MINPOP)) && (gsl_vector_get(x,nlte_size-1) > (1.1 * MINPOP)))
      {
        test_ratio_upper = get_groundlevelpop(modelgridindex,element,ion+1) * modelgrid[modelgridindex].composition[element].partfunct[ion+1]
          / stat_weight(element,ion+1,0) / gsl_vector_get(x,nlte_size-1);
        if (test_ratio_upper < 1)
        {
          test_ratio_upper = 1. / test_ratio_upper;
        }
      }
      else
      {
        test_ratio_upper = 0.0;
      }

      if (test_ratio_upper > test_ratio)
      {
        test_ratio = test_ratio_upper;
      }

      printout("The test ratios are %g %g. Passing %g.\n", get_groundlevelpop(modelgridindex,element,ion)/gsl_vector_get(x,0),
        get_groundlevelpop(modelgridindex,element,ion+1)*modelgrid[modelgridindex].composition[element].partfunct[ion+1]/stat_weight(element,ion+1,0)/gsl_vector_get(x,nlte_size-1),
        test_ratio);
      //	  printout("The test ratio is %g. Passing %g.\n", get_groundlevelpop(modelgridindex,element,ion+1)*modelgrid[modelgridindex].composition[element].partfunct[ion+1]/stat_weight(element,ion+1,0)/gsl_vector_get(x,nlte_size-1), test_ratio);


      //printout("The top five excited states were %g, %g, %g, %g and %g.\n",gsl_vector_get(x,nlevels_nlte-4),gsl_vector_get(x,nlevels_nlte-3),gsl_vector_get(x,nlevels_nlte-2),gsl_vector_get(x,nlevels_nlte-1),gsl_vector_get(x,nlevels_nlte));

      //printout("The first five excited states were %g, %g, %g, %g and %g.\n",gsl_vector_get(x,1),gsl_vector_get(x,2),gsl_vector_get(x,3),gsl_vector_get(x,4),gsl_vector_get(x,5));

      //lag = 0.0;

      //if ((check=(fabs((1.0 - test_ratio)/(1.0 + test_ratio))) > 0.05) && (check < 0.2))
      //{
      //   lag=1.0;
      //}

      // Write the NLTE level populations to the array
      for (int level = 1; level < nlevels_nlte+1; level++)
      {
        modelgrid[modelgridindex].nlte_pops[nlte_start+level-1] = gsl_vector_get(x,level)/get_rho(modelgridindex);
        //printout("I have interfered with index %d.\n", nlte_start+level-1);
        //modelgrid[modelgridindex].nlte_pops[nlte_start+level-1] = ((lag*modelgrid[modelgridindex].nlte_pops[nlte_start+level-1]) + gsl_vector_get(x,level))/(lag + 1.0)/get_rho(modelgridindex);
      }
      // If there is a superlevel then write that too

      if (super_level == 1)
      {
        //printout("I thought the super level was: %g\n", modelgrid[modelgridindex].nlte_pops[nlte_start+nlevels_nlte]);

        modelgrid[modelgridindex].nlte_pops[nlte_start+nlevels_nlte] = gsl_vector_get(x,nlevels_nlte+1) / modelgrid[modelgridindex].rho / superlevel_partition;
        //	      modelgrid[modelgridindex].nlte_pops[nlte_start+nlevels_nlte] = ((lag*modelgrid[modelgridindex].nlte_pops[nlte_start+nlevels_nlte]) + gsl_vector_get(x,nlevels_nlte+1))/(lag + 1.0)/modelgrid[modelgridindex].rho/superlevel_partition;

        //printout("Now I think it is: %g\n", modelgrid[modelgridindex].nlte_pops[nlte_start+nlevels_nlte]);
        //printout("I also interfered with index %d.\n", nlte_start+nlevels_nlte);
      }

      //printout("I had a ground level pop of %g, a part fn of %g and therefore an ion pop of %g\n", modelgrid[modelgridindex].composition[element].groundlevelpop[ion], modelgrid[modelgridindex].composition[element].partfunct[ion], modelgrid[modelgridindex].composition[element].partfunct[ion]*modelgrid[modelgridindex].composition[element].groundlevelpop[ion]/stat_weight(element,ion,0));

      modelgrid[modelgridindex].composition[element].groundlevelpop[ion] = gsl_vector_get(x,0);

      //printout("Now I have a ground level pop of %g, a part fn of %g and therefore an ion pop of %g\n", modelgrid[modelgridindex].composition[element].groundlevelpop[ion], calculate_partfunct(element, ion, modelgridindex), calculate_partfunct(element, ion, modelgridindex)*modelgrid[modelgridindex].composition[element].groundlevelpop[ion]/stat_weight(element,ion,0));

      //printout("For the ion above, I had ground state of  %g, a part fn of %g and therefore an ion pop of %g\n",modelgrid[modelgridindex].composition[element].groundlevelpop[ion+1], modelgrid[modelgridindex].composition[element].partfunct[ion+1], modelgrid[modelgridindex].composition[element].partfunct[ion+1]*modelgrid[modelgridindex].composition[element].groundlevelpop[ion+1]/stat_weight(element,ion+1,0));

      //printout("Now I get that it should have an ion pop of %g\n", gsl_vector_get(x,nlte_size-1));

      //printout("Accordingly, my phi-factor is: %g\n", calculate_partfunct(element, ion, modelgridindex)*modelgrid[modelgridindex].composition[element].groundlevelpop[ion]/stat_weight(element,ion,0)/gsl_vector_get(x,nlte_size-1)/nne);
      //printout("The ionization solver gives phi = %g\n", phi(element, ion, modelgridindex));

      //printout("From ion fract, the lower and upper ion pops should be %g and %g\n", ionfract(element, ion, modelgridindex, nne)*get_abundance(modelgridindex,element)/ elements[element].mass * get_rho(modelgridindex), ionfract(element, ion+1, modelgridindex, nne)*get_abundance(modelgridindex,element)/ elements[element].mass * get_rho(modelgridindex));
      //printout("I think that the element population is: %g (from %g and %g)\n", get_abundance(modelgridindex,element)/ elements[element].mass * get_rho(modelgridindex), get_abundance(modelgridindex,element), get_rho(modelgridindex));
      //printout("I currently think that the top ion is: %d\n", elements_uppermost_ion[tid][element]);

      gsl_permutation_free(p);
      //printout("I freed up p\n");
      gsl_vector_free(x);
      //printout("I freed up x\n");
      //printf("rate %p balane %p NULL %p\n", rate_matrix, balance_vector, NULL);

      free(rate_matrix);
      //printout("I freed up rate_matrix\n");
      free(balance_vector);
      //printout("I freed up balance_vector\n");
    }
    else
    {
      //STUFF FOR "NOT USING" CASE

      nltepop_reset_element(modelgridindex, element);
      test_ratio = 0.0;
    }
    return test_ratio;
  }
  else
  {
    return 0; //Case for ion with only one level
  }
}


double superlevel_boltzmann(const int modelgridindex, const int element, const int ion, const int level)
{
  const int superlevel_index = get_nlevels_nlte(element,ion) + 1;
  // const double T_exc = get_TJ(modelgridindex);
  const double T_exc = get_Te(modelgridindex);
  const double E_level = epsilon(element, ion, level);
  const double E_superlevel = epsilon(element, ion, superlevel_index);

  return stat_weight(element, ion, level) / stat_weight(element, ion, superlevel_index) * exp(- (E_level - E_superlevel) / KB / T_exc);
}


void nltepop_open_file(const int my_rank)
{
  char filename[100];
  sprintf(filename, "nlte_%.4d.out",my_rank);
  nlte_file = fopen_required(filename, "w");
  fprintf(nlte_file, "%8s %14s %2s %9s %5s %11s %11s %11s\n",
          "timestep", "modelgridindex", "Z", "ion_stage", "level", "n_LTE", "n_NLTE", "ion_popfrac");
}


void nltepop_close_file(void)
{
  fclose(nlte_file);
}


void nltepop_write_to_file(const int modelgridindex, const int timestep)
{
  if (initial_iteration) // NLTE solver hasn't been run yet
    return;

  // fprintf(nlte_file,"#timestep %d modelgridindex %d T_R %g T_e %g W %g T_J %g nne %g\n",
  //         timestep, n, get_TR(n), get_Te(n), get_W(n), get_TJ(n), get_nne(n));

  for (int element = 0; element < nelements; element++)
  {
    const int nions = get_nions(element);
    const int atomic_number = get_element(element);

    for (int ion = 0; ion < nions; ion++)
    {
      const int nlevels_nlte = get_nlevels_nlte(element, ion);
      const int ion_first_nlte = elements[element].ions[ion].first_nlte;
      const int ion_stage = get_ionstage(element, ion);

      const int nsuperlevels = ion_has_superlevel(element, ion) ? 1 : 0;

      for (int level = 0; level <= nlevels_nlte + nsuperlevels; level++)
      {
        double nnlevellte = calculate_levelpop_lte(modelgridindex, element, ion, level);
        double nnlevelnlte;

        // use "%8d %14d %2d %9d " for fixed width
        fprintf(nlte_file, "%d %d %d %d ", timestep, modelgridindex, atomic_number, ion_stage);
        if (level <= nlevels_nlte)
        {
          fprintf(nlte_file, "%d ", level);

          if (level == 0)
          {
            nnlevelnlte = get_groundlevelpop(modelgridindex, element, ion);
          }
          else
          {
            nnlevelnlte =  (
               modelgrid[modelgridindex].nlte_pops[ion_first_nlte + level - 1] * modelgrid[modelgridindex].rho);
          }
        }
        else
        {
          // superlevel, so add the populations of all other levels in the superlevel
          const double slpopfactor = (
            modelgrid[modelgridindex].nlte_pops[ion_first_nlte + nlevels_nlte] * modelgrid[modelgridindex].rho);

          nnlevellte = 0;
          double superlevel_partfunc = 0;
          fprintf(nlte_file, "%d ", -1);
          for (int level_sl = nlevels_nlte + 1; level_sl < get_nlevels(element, ion); level_sl++)
          {
            nnlevellte += calculate_levelpop_lte(modelgridindex, element, ion, level_sl);
            superlevel_partfunc += superlevel_boltzmann(modelgridindex, element, ion, level_sl);
          }

          nnlevelnlte = slpopfactor * superlevel_partfunc;

          // printout("nltepop_write_to_file: The Z=%d ionstage %d superlevel population is %g with rho %g and superlevel_partfunc %g Te %g scaled pop stored as %g\n", get_element(element), get_ionstage(element, ion), nnlevelnlte, modelgrid[modelgridindex].rho, superlevel_partfunc, get_Te(modelgridindex), modelgrid[modelgridindex].nlte_pops[ion_first_nlte + nlevels_nlte]);
        }

        const double ion_popfrac = nnlevelnlte / ionstagepop(modelgridindex, element, ion);
        fprintf(nlte_file, "%11.5e %11.5e %11.5e\n", nnlevellte, nnlevelnlte, ion_popfrac);
      }
    }
  }

  fflush(nlte_file);
}


void nltepop_write_restart_data(FILE *restart_file)
{
  if (!NLTE_POPS_ON)
    return;

  printout("NLTE populations, ");

  fprintf(restart_file, "%d\n", 75618527); // special number marking the beginning of nlte data
  fprintf(restart_file, "%d\n", total_nlte_levels);

  for (int modelgridindex = 0; modelgridindex < MMODELGRID; modelgridindex++)
  {
    if (get_numassociatedcells(modelgridindex) > 0)
    {
      fprintf(restart_file, "%d\n", modelgridindex);
      for (int element = 0; element < nelements; element++)
      {
        const int nions = get_nions(element);
        for (int ion = 0; ion < nions; ion++)
        {
          fprintf(restart_file, "%d %g %g %lg\n",
                  ion,
                  modelgrid[modelgridindex].composition[element].groundlevelpop[ion],
                  modelgrid[modelgridindex].composition[element].partfunct[ion],
                  modelgrid[modelgridindex].cooling[element].contrib[ion]);
        }
      }
      for (int nlteindex = 0; nlteindex < total_nlte_levels; nlteindex++)
      {
        fprintf(restart_file, "%lg ", modelgrid[modelgridindex].nlte_pops[nlteindex]);
      }
    }
  }
}


void nltepop_read_restart_data(FILE *restart_file)
{
  if (!NLTE_POPS_ON)
    return;

  printout("Reading restart data for NLTE populations\n");

  int code_check;
  fscanf(restart_file, "%d\n", &code_check);
  if (code_check != 75618527)
  {
    printout("ERROR: Beginning of NLTE restart data not found!");
    abort();
  }

  int total_nlte_levels_in;
  fscanf(restart_file, "%d\n", &total_nlte_levels_in);
  if (total_nlte_levels_in != total_nlte_levels)
  {
    printout("ERROR: Expected NLTE levels but found %d in restart file", total_nlte_levels, total_nlte_levels_in);
    abort();
  }

  for (int modelgridindex = 0; modelgridindex < MMODELGRID; modelgridindex++)
  {
    if (get_numassociatedcells(modelgridindex) > 0)
    {
      int mgi_in;
      fscanf(restart_file, "%d\n", &mgi_in);
      if (mgi_in != modelgridindex)
      {
        printout("ERROR: expected data for cell %d but found cell %d\n", modelgridindex, mgi_in);
        abort();
      }

      for (int element = 0; element < nelements; element++)
      {
        const int nions = get_nions(element);
        for (int ion = 0; ion < nions; ion++)
        {
          int ion_in;
          fscanf(restart_file, "%d %g %g %lg\n",
                 &ion_in,
                 &modelgrid[modelgridindex].composition[element].groundlevelpop[ion],
                 &modelgrid[modelgridindex].composition[element].partfunct[ion],
                 &modelgrid[modelgridindex].cooling[element].contrib[ion]);
          if (ion_in != ion)
          {
            printout("ERROR: expected data for ion %d but found ion %d\n", ion, ion_in);
            abort();
          }
        }
      }
      for (int nlteindex = 0; nlteindex < total_nlte_levels; nlteindex++)
      {
        fscanf(restart_file, "%lg ", &modelgrid[modelgridindex].nlte_pops[nlteindex]);
      }
    }
  }
}