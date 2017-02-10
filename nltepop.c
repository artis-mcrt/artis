#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>
#include "sn3d.h"
#include "atomic.h"
#include "grid_init.h"
#include "ltepop.h"
#include "macroatom.h"
#include "nonthermal.h"
#include "nltepop.h"
#include "ratecoeff.h"
#include "update_grid.h"


static int get_nlte_vector_index(int element_in, int ion_in, int level_in)
{
  int index = 0;
  for (int ion = 0; ion < ion_in; ion++)
  {
    const int nlevels_nlte = get_nlevels_nlte(element_in,ion);

    index += nlevels_nlte + 1;
    if (nlevels_nlte != (get_nlevels(element_in,ion) - 1))
      index++; // there's a superlevel here
  }
  const int nlevels_nlte = get_nlevels_nlte(element_in,ion_in);
  if (is_nlte(element_in, ion_in, level_in) == true)
    index += level_in;
  else
    index += nlevels_nlte + 1; //the index of the superlevel

  return index;
}


static void get_ion_level_of_nlte_vector_index(int index, int element, int *ion, int *level)
{
  for (int dion = 0; dion < get_nions(element); dion++)
  {
    for (int dlevel = 0; dlevel < get_nlevels(element,dion); dlevel++)
    {
      if (get_nlte_vector_index(element,dion,dlevel) == index)
      {
        *ion = dion;
        *level = dlevel;
        return;
      }
    }
  }
}


static void eliminate_nlte_matrix_rowcol(int index, int gs_index, gsl_matrix *const rate_matrix, gsl_vector *const balance_vector)
{
  const gsl_matrix rate_matrix_var = *rate_matrix;
  const int colcount = rate_matrix_var.size2;
  for (int column = 0; column < colcount; column++)
    gsl_matrix_set(rate_matrix,index,column,0.0);

  const int rowcount = rate_matrix_var.size1;
  for (int row = 1; row < rowcount; row++)
    gsl_matrix_set(rate_matrix,row,index,0.0);

  gsl_matrix_set(rate_matrix,index,gs_index,-1.0);
  gsl_matrix_set(rate_matrix,index,index,1.0);
  gsl_vector_set(balance_vector,index,0.0);
}


static void filter_nlte_matrix(int element, gsl_matrix *const rate_matrix, gsl_vector *const balance_vector, gsl_vector *const pop_norm_factors)
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
      const double element_value = fabs(gsl_matrix_get(rate_matrix,index,column));
      if (element_value > row_max)
        row_max = element_value;
    }
    double col_max = 0.0;
    for (int row = 1; row < nlte_dimension; row++) //skip the normalisation row 0
    {
      const double element_value = fabs(gsl_matrix_get(rate_matrix,row,index));
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
        double gs_index = get_nlte_vector_index(element,ion,0);
        eliminate_nlte_matrix_rowcol(index,gs_index,rate_matrix,balance_vector);
        // printout("(forcing LTE population)");
      }
    }
    // printout("\n");
  }
}


static double get_total_rate_in(int index_to, gsl_matrix *rate_matrix, gsl_vector *popvec)
{
  const gsl_matrix rate_matrix_var = *rate_matrix;
  gsl_vector *rates_in_vec = gsl_vector_alloc(rate_matrix_var.size1);
  gsl_vector_view row_view = gsl_matrix_row(rate_matrix, index_to);
  gsl_vector_memcpy(rates_in_vec, &row_view.vector);
  gsl_vector_set(rates_in_vec, index_to, 0.);
  gsl_vector_mul(rates_in_vec, popvec);
  const double total_rate_in = gsl_blas_dasum(rates_in_vec);
  gsl_vector_free(rates_in_vec);
  return total_rate_in;
}


static double get_total_rate_out(int index_from, gsl_matrix *rate_matrix, gsl_vector *popvec)
{
  const gsl_matrix rate_matrix_var = *rate_matrix;
  gsl_vector *rates_out_vec = gsl_vector_alloc(rate_matrix_var.size2);
  gsl_vector_view col_view = gsl_matrix_column(rate_matrix, index_from);
  gsl_vector_memcpy(rates_out_vec, &col_view.vector);
  gsl_vector_set(rates_out_vec, index_from, 0.);
  gsl_vector_scale(rates_out_vec, gsl_vector_get(popvec, index_from));
  const double total_rate_out = gsl_blas_dasum(rates_out_vec);
  gsl_vector_free(rates_out_vec);
  return total_rate_out;
}


static void print_level_rates(int modelgridindex, int element, int selected_ion, int selected_level,
                              gsl_vector *popvec, gsl_matrix *rate_matrix_rad_bb,
                              gsl_matrix *rate_matrix_coll_bb, gsl_matrix *rate_matrix_rad_bf,
                              gsl_matrix *rate_matrix_coll_bf, gsl_matrix *rate_matrix_ntcoll_bf)

{
  if (element > nelements - 1 || selected_ion > get_nions(element) - 1 || selected_level > get_nlevels_nlte(element, selected_ion) - 1)
    return;
  const gsl_vector popvector = *popvec;
  const int nlte_dimension = popvector.size;
  const int selected_ionstage = get_ionstage(element, selected_ion);
  const int selected_index = get_nlte_vector_index(element, selected_ion, selected_level);
  const double pop_selectedlevel = gsl_vector_get(popvec, selected_index);
  printout("NLTE level diagnostics for index %d corresponding to ion_stage %d level %d (rates into and out of this level)\n", selected_index, selected_ionstage, selected_level);

  const double rad_bb_in_total = get_total_rate_in(selected_index, rate_matrix_rad_bb, popvec);
  const double coll_bb_in_total = get_total_rate_in(selected_index, rate_matrix_coll_bb, popvec);
  const double rad_bf_in_total = get_total_rate_in(selected_index, rate_matrix_rad_bf, popvec);
  const double coll_bf_in_total = get_total_rate_in(selected_index, rate_matrix_coll_bf, popvec);
  const double ntcoll_bf_in_total = get_total_rate_in(selected_index, rate_matrix_ntcoll_bf, popvec);
  const double total_rate_in = rad_bb_in_total + coll_bb_in_total + rad_bf_in_total + coll_bf_in_total + ntcoll_bf_in_total;
  printout("  TOTAL rates in:             rad_bb_in  %8.2e coll_bb_in  %8.2e rad_bf_in  %8.2e coll_bf_in  %8.2e ntcoll_bf_in  %8.2e\n",
           rad_bb_in_total, coll_bb_in_total, rad_bf_in_total, coll_bf_in_total, ntcoll_bf_in_total);

  const double rad_bb_out_total = get_total_rate_out(selected_index, rate_matrix_rad_bb, popvec);
  const double coll_bb_out_total = get_total_rate_out(selected_index, rate_matrix_coll_bb, popvec);
  const double rad_bf_out_total = get_total_rate_out(selected_index, rate_matrix_rad_bf, popvec);
  const double coll_bf_out_total = get_total_rate_out(selected_index, rate_matrix_coll_bf, popvec);
  const double ntcoll_bf_out_total = get_total_rate_out(selected_index, rate_matrix_ntcoll_bf, popvec);
  const double total_rate_out = rad_bb_out_total + coll_bb_out_total + rad_bf_out_total + coll_bf_out_total + ntcoll_bf_out_total;
  printout("  TOTAL rates out:            rad_bb_out %8.2e coll_bb_out %8.2e rad_bf_out %8.2e coll_bf_out %8.2e ntcoll_bf_out %8.2e\n",
          rad_bb_out_total, coll_bb_out_total, rad_bf_out_total, coll_bf_out_total, ntcoll_bf_out_total);

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
    const double rad_bf_in = gsl_matrix_get(rate_matrix_rad_bf, selected_index, index) * pop;
    const double rad_bf_out = gsl_matrix_get(rate_matrix_rad_bf, index, selected_index) * pop_selectedlevel;
    const double coll_bf_in = gsl_matrix_get(rate_matrix_coll_bf, selected_index, index) * pop;
    const double coll_bf_out = gsl_matrix_get(rate_matrix_coll_bf, index, selected_index) * pop_selectedlevel;
    const double ntcoll_bf_in = gsl_matrix_get(rate_matrix_ntcoll_bf, selected_index, index) * pop;
    const double ntcoll_bf_out = gsl_matrix_get(rate_matrix_ntcoll_bf, index, selected_index) * pop_selectedlevel;

    const bool nonzero_rate_in = (fabs(rad_bb_in) > 0. || fabs(coll_bb_in) > 0. || fabs(rad_bf_in) > 0. || fabs(coll_bf_in) > 0. || fabs(ntcoll_bf_in) > 0.);
    const bool nonzero_rate_out = (fabs(rad_bb_out) > 0. || fabs(coll_bb_out) > 0. || fabs(rad_bf_out) > 0. || fabs(coll_bf_out) > 0. || fabs(ntcoll_bf_out) > 0.);
    if (nonzero_rate_in || nonzero_rate_out)
    {
      const double epsilon_trans = fabs(epsilon(element,ion,level) - epsilon(element,selected_ion,selected_level));
      const double nu_trans = epsilon_trans / H;
      const double lambda = 1e8 * CLIGHT / nu_trans; // should be in Angstroms
      const double level_rate_in = rad_bb_in + coll_bb_in + rad_bf_in + coll_bf_in + ntcoll_bf_in;
      const double level_rate_out = rad_bb_out + coll_bb_out + rad_bf_out + coll_bf_out + ntcoll_bf_out;
      const int level_percent_in = round(level_rate_in / total_rate_in * 100);
      const int level_percent_out = round(level_rate_out / total_rate_out * 100);

      printout("  ion_stage %d level %4d (%3d%% rate in)  rad_bb_in  %8.2e coll_bb_in  %8.2e rad_bf_in  %8.2e coll_bf_in  %8.2e ntcoll_bf_in  %8.2e lambda %5.1f\n",
               ionstage, level, level_percent_in, rad_bb_in, coll_bb_in, rad_bf_in, coll_bf_in, ntcoll_bf_in, lambda);
      printout("  ion_stage %d level %4d (%3d%% rate out) rad_bb_out %8.2e coll_bb_out %8.2e rad_bf_out %8.2e coll_bf_out %8.2e ntcoll_bf_out %8.2e lambda %5.1f\n",
               ionstage, level, level_percent_out, rad_bb_out, coll_bb_out, rad_bf_out, coll_bf_out, ntcoll_bf_out, lambda);
    }
  }
  printout("\n");
}


void solve_nlte_pops_element(int element, int modelgridindex, int timestep)
//solves for nlte correction factors to level populations for levels in all ions of an element
{
  const double t_mid = time_step[timestep].mid;
  const int nions = get_nions(element);
  const int atomic_number = get_element(element);

  if (get_abundance(modelgridindex, element) > 0.0)
  {
    printout("Solving for NLTE populations in cell %d at timestep %d for element %d (Z=%d)\n",modelgridindex,timestep,element,atomic_number);
    //set_TR(modelgridindex,3000); //TODO: remove after testing complete
    //set_W(modelgridindex,1.0); //TODO: remove after testing complete
    //printout("T_E %g T_R was %g, setting to 3000 \n",get_Te(modelgridindex),get_TR(modelgridindex));

    int nlte_dimension = 0;
    double *const superlevel_partfunc = calloc(nions,sizeof(double)); //space is allocated for every ion, whether or not it has a superlevel
    for (int ion = 0; ion < nions; ion++)
    {
      const int nlevels_nlte = get_nlevels_nlte(element,ion);

      //this is the total number of nlte_levels (i.e. the size of the
      //storage). Our rate matrix will need to be of this dimension +2: the
      //ground state, the "super level".
      //If there's no super level needed then we only need +1
      if (nlevels_nlte == (get_nlevels(element, ion) - 1))
      {
        nlte_dimension += nlevels_nlte + 1;
        printout("NLTE: setting up ion %d, which contributes %d to the vector dimension (no super level)\n",
                 ion, nlevels_nlte + 1);
      }
      else
      {
        for (int level = 1; level < get_nlevels(element,ion); level++)
        {
          if (is_nlte(element, ion, level) == false)
            superlevel_partfunc[ion] += superlevel_boltzmann(modelgridindex,element,ion,level);
        }
        nlte_dimension += nlevels_nlte + 2;
        printout("NLTE: setting up ion %d, which contributes %d to the vector dimension (including superlevel with partfunc %g)\n",
                 ion, nlevels_nlte + 2, superlevel_partfunc[ion]);
      }
    }

    printout("NLTE: the vector dimension is %d.\n", nlte_dimension);

    gsl_matrix *const rate_matrix_rad_bb = gsl_matrix_calloc(nlte_dimension, nlte_dimension);
    gsl_matrix *const rate_matrix_coll_bb = gsl_matrix_calloc(nlte_dimension, nlte_dimension);
    gsl_matrix *const rate_matrix_rad_bf = gsl_matrix_calloc(nlte_dimension, nlte_dimension);
    gsl_matrix *const rate_matrix_coll_bf = gsl_matrix_calloc(nlte_dimension, nlte_dimension);
    gsl_matrix *const rate_matrix_ntcoll_bf = gsl_matrix_calloc(nlte_dimension, nlte_dimension);
    gsl_vector *const balance_vector = gsl_vector_calloc(nlte_dimension);

    if (!balance_vector)
    {
      printout("Cannot allocate NLTE rate matrix/balance vector memory.\n");
      abort();
    }

    for (int ion = 0; ion < nions; ion++)
    {
      double Y_nt = 0.; // rate coeff for collisional ionization by non-thermal electrons
      int upper_groundstate_index = -1;
      if (NT_ON && ion < nions - 1)
      {
        upper_groundstate_index = get_nlte_vector_index(element, ion + 1, 0);
        Y_nt = nt_ionization_ratecoeff(modelgridindex, element, ion);
        if (Y_nt < 0.)
        {
          printout("WARNING: Negative NT_ionization rate from ion_stage %d\n", get_ionstage(element, ion));
        }
      }

      const int nlevels = get_nlevels(element, ion);
      for (int level = 0; level < nlevels; level++)
      {
        const double epsilon_current = epsilon(element, ion, level);
        const int nionisinglevels = get_bfcontinua(element, ion);

        double s_renorm;
        if ((level != 0) && (is_nlte(element, ion, level) == false))
          s_renorm = superlevel_boltzmann(modelgridindex, element, ion, level) / superlevel_partfunc[ion];
        else
          s_renorm = 1.0;

        const int level_index = get_nlte_vector_index(element, ion, level);

        // de-excitation
        const int ndowntrans = elements[element].ions[ion].levels[level].downtrans[0].targetlevel;
        for (int i = 1; i <= ndowntrans; i++)
        {
          const int lower = elements[element].ions[ion].levels[level].downtrans[i].targetlevel;
          const double epsilon_target = elements[element].ions[ion].levels[level].downtrans[i].epsilon;
          const int lineindex = elements[element].ions[ion].levels[level].downtrans[i].lineindex;
          const double epsilon_trans = epsilon_current - epsilon_target;

          const double R = rad_deexcitation_ratecoeff(modelgridindex, element, ion, level, lower, epsilon_trans, lineindex, t_mid);
          const double C = col_deexcitation_ratecoeff(modelgridindex, epsilon_trans, lineindex);

          const int upper_index = level_index;
          const int lower_index = get_nlte_vector_index(element,ion,lower);

          *gsl_matrix_ptr(rate_matrix_rad_bb, upper_index, upper_index) -= R * s_renorm;
          *gsl_matrix_ptr(rate_matrix_rad_bb, lower_index, upper_index) += R * s_renorm;
          *gsl_matrix_ptr(rate_matrix_coll_bb, upper_index, upper_index) -= C * s_renorm;
          *gsl_matrix_ptr(rate_matrix_coll_bb, lower_index, upper_index) += C * s_renorm;
          if (R + C < 0)
            printout("WARNING: Negative de-excitation rate from ion_stage %d level %d to level %d\n",get_ionstage(element, ion), level, lower);
        }

        // excitation
        const int nuptrans = elements[element].ions[ion].levels[level].uptrans[0].targetlevel;
        for (int i = 1; i <= nuptrans; i++)
        {
          const int upper = elements[element].ions[ion].levels[level].uptrans[i].targetlevel;
          const double epsilon_target = elements[element].ions[ion].levels[level].uptrans[i].epsilon;
          const int lineindex = elements[element].ions[ion].levels[level].uptrans[i].lineindex;
          const double epsilon_trans = epsilon_target - epsilon_current;

          const double R = rad_excitation_ratecoeff(modelgridindex, element, ion, level, upper, epsilon_trans, lineindex, t_mid);
          const double C = col_excitation_ratecoeff(modelgridindex, lineindex, epsilon_trans);

          // if ((element == 0) && (ion == 1) && ((level <= 5) || (level == 35)) && (upper >= 74) && (upper <= 77))
          // {
          //   const double tau_sobolev = get_tau_sobolev(modelgridindex, lineindex, t_mid);
          //   printout("timestep %d lower %d upper %d tau_sobolev=%g\n", timestep, level, upper, tau_sobolev);
          // }

          const int lower_index = level_index;
          const int upper_index = get_nlte_vector_index(element,ion,upper);

          *gsl_matrix_ptr(rate_matrix_rad_bb, lower_index, lower_index) -= R * s_renorm;
          *gsl_matrix_ptr(rate_matrix_rad_bb, upper_index, lower_index) += R * s_renorm;
          *gsl_matrix_ptr(rate_matrix_coll_bb, lower_index, lower_index) -= C * s_renorm;
          *gsl_matrix_ptr(rate_matrix_coll_bb, upper_index, lower_index) += C * s_renorm;
          if (R + C < 0)
            printout("WARNING: Negative excitation rate from ion %d level %d to level %d\n",get_ionstage(element,ion),level,upper);
        }

        // collisional ionization by non-thermal electrons
        if (Y_nt > 0.)
        {
          const int lower_index = level_index;

          *gsl_matrix_ptr(rate_matrix_ntcoll_bf, lower_index, lower_index) -= Y_nt * s_renorm;
          *gsl_matrix_ptr(rate_matrix_ntcoll_bf, upper_groundstate_index, lower_index) += Y_nt * s_renorm;
        }

        // thermal collisional ionization, photoionisation and recombination processes
        if (ion < nions-1 && level < nionisinglevels)  //&& get_ionstage(element,ion) < get_element(element)+1)
        {
          // ionization
          const int lower_index = level_index;

          for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion,level); phixstargetindex++)
          {
            const int upper = get_phixsupperlevel(element,ion,level,phixstargetindex);
            const int upper_index = get_nlte_vector_index(element,ion+1,upper);
            const double epsilon_trans = epsilon(element,ion+1,upper) - epsilon_current;
            const double R = get_corrphotoioncoeff(element, ion, level, phixstargetindex, modelgridindex);
            const double C = col_ionization_ratecoeff(modelgridindex, element, ion, level, phixstargetindex, epsilon_trans);

            *gsl_matrix_ptr(rate_matrix_rad_bf, lower_index, lower_index) -= R * s_renorm;
            *gsl_matrix_ptr(rate_matrix_rad_bf, upper_index, lower_index) += R * s_renorm;
            *gsl_matrix_ptr(rate_matrix_coll_bf, lower_index, lower_index) -= C * s_renorm;
            *gsl_matrix_ptr(rate_matrix_coll_bf, upper_index, lower_index) += C * s_renorm;

            if (R + C < 0)
              printout("WARNING: Negative ionization rate from ion_stage %d level %d phixstargetindex %d\n",get_ionstage(element,ion),level,phixstargetindex);
          }

          // recombination
          for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion,level); phixstargetindex++)
          {
            const int upper = get_phixsupperlevel(element,ion,level,phixstargetindex);
            const int upper_index = get_nlte_vector_index(element,ion+1,upper);
            const double epsilon_trans = epsilon(element,ion+1,upper) - epsilon_current;
            const double R = rad_recombination_ratecoeff(modelgridindex, element, ion + 1, upper, level);
            const double C = col_recombination_ratecoeff(modelgridindex, element, ion + 1, upper, level, epsilon_trans);

            double s_renorm_upper = 1.0;
            if ((upper != 0) && (is_nlte(element,ion+1,upper) == false))
              s_renorm_upper = superlevel_boltzmann(modelgridindex,element,ion+1,upper) / superlevel_partfunc[ion+1];

            *gsl_matrix_ptr(rate_matrix_rad_bf, upper_index, upper_index) -= R * s_renorm_upper;
            *gsl_matrix_ptr(rate_matrix_rad_bf, lower_index, upper_index) += R * s_renorm_upper;
            *gsl_matrix_ptr(rate_matrix_coll_bf, upper_index, upper_index) -= C * s_renorm_upper;
            *gsl_matrix_ptr(rate_matrix_coll_bf, lower_index, upper_index) += C * s_renorm_upper;
            if (R + C < 0)
              printout("WARNING: Negative recombination rate to ion_stage %d level %d phixstargetindex %d\n",get_ionstage(element,ion),level,phixstargetindex);
          }

        } // if (level is an ionising level)

      } // level loop
    } // ion loop

    // sum the matricies for each transition process to get a total rate matrix
    gsl_matrix *const rate_matrix = gsl_matrix_calloc(nlte_dimension, nlte_dimension);
    gsl_matrix_add(rate_matrix, rate_matrix_rad_bb);
    gsl_matrix_add(rate_matrix, rate_matrix_coll_bb);
    gsl_matrix_add(rate_matrix, rate_matrix_rad_bf);
    gsl_matrix_add(rate_matrix, rate_matrix_coll_bf);
    gsl_matrix_add(rate_matrix, rate_matrix_ntcoll_bf);

    // replace the first row of the matrix and balance vector with the normalisation
    // constraint on the total element population
    gsl_vector_view first_row_view = gsl_matrix_row(rate_matrix, 0);
    gsl_vector_set_all(&first_row_view.vector, 1.0);
    // set first balance vector entry to the element population (all other entries will be zero)
    double element_population = get_abundance(modelgridindex,element) / elements[element].mass * get_rho(modelgridindex);
    gsl_vector_set(balance_vector, 0, element_population);

    // calculate the normalisation factors and apply them to the matrix
    // columns and balance vector elements
    gsl_vector *pop_norm_factor_vec = gsl_vector_calloc(nlte_dimension);
    for (int column = 0; column < nlte_dimension; column++)
    {
      int ion, level;
      get_ion_level_of_nlte_vector_index(column,element,&ion,&level);

      gsl_vector_set(pop_norm_factor_vec, column, calculate_levelpop_lte(modelgridindex,element,ion,level));

      if ((level != 0) && (is_nlte(element,ion,level) == false))
      {
        // level is a superlevel, so add populations of higher levels to the norm factor
        for (int dummylevel = level + 1; dummylevel < get_nlevels(element,ion); dummylevel++)
        {
          if (is_nlte(element,ion,dummylevel) == false)
          {
            *gsl_vector_ptr(pop_norm_factor_vec, column) += calculate_levelpop_lte(modelgridindex,element,ion,dummylevel);
          }
        }
        // NOTE: above calculation is not always equal to the sum of LTE populations
        // since calculate_levelpop_lte imposes MINPOP minimum
        printout("superlevel norm factor index %d is %g, partfunc is %g, partfunc*levelpop(SL)/g(SL) %g\n",
                 column, gsl_vector_get(pop_norm_factor_vec, column), superlevel_partfunc[ion],
                 superlevel_partfunc[ion] * calculate_levelpop_lte(modelgridindex,element,ion,level) / stat_weight(element,ion,level));
      }

      // apply the normalisation factor to this column in the rate_matrix
      gsl_vector_view column_view = gsl_matrix_column(rate_matrix,column);
      gsl_vector_scale(&column_view.vector, gsl_vector_get(pop_norm_factor_vec, column));
    }

    /*printout("Rate matrix | balance vector:\n");
    for (int row = 0; row < nlte_dimension; row++)
    {
      for (int column = 0; column < nlte_dimension; column++)
      {
        char str[15];
        sprintf(str, "%+.1e ", rate_matrix[row * nlte_dimension + column]);
        printout(str);
      }
      printout("| ");
      char str[15];
      sprintf(str, "%+.1e\n", balance_vector[row]);
      printout(str);
    }
    printout("\n");*/

    // eliminate barely-interacting levels from the NLTE matrix by removing
    // their interactions and setting their normalised populations (probably departure coeff) to 1.0
    filter_nlte_matrix(element,rate_matrix,balance_vector,pop_norm_factor_vec);

    // make a copy of the rate matrix for the LU decomp
    gsl_matrix *rate_matrix_LU_decomp = gsl_matrix_alloc(nlte_dimension, nlte_dimension);
    gsl_matrix_memcpy(rate_matrix_LU_decomp, rate_matrix);

    gsl_permutation *p = gsl_permutation_alloc(nlte_dimension);

    int s; //sign of the transformation
    gsl_linalg_LU_decomp(rate_matrix_LU_decomp, p, &s);

    gsl_vector *x = gsl_vector_alloc(nlte_dimension); // population solution vector

    // solve matrix equation: rate_matrix * x = balance_vector for x (population vector)
    gsl_linalg_LU_solve(rate_matrix_LU_decomp, p, balance_vector, x);

    //gsl_linalg_HH_solve (&m.matrix, &b.vector, x);

    const double TOLERANCE = 1e-40;

    double error_best = -1.;
    gsl_vector *x_best = gsl_vector_alloc(nlte_dimension); //population solution vector with lowest error
    gsl_vector *gsl_work_vector = gsl_vector_alloc(nlte_dimension);
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
        gsl_vector_memcpy(x_best,x);
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
      printout("NLTE solver LU_refine: After %d iterations, keeping solution vector that had a max residual of %g\n",iteration,error_best);
      gsl_vector_memcpy(x,x_best);
    }
    gsl_vector_free(gsl_work_vector);
    gsl_vector_free(x_best);
    gsl_permutation_free(p);

    gsl_vector *popvec = gsl_vector_alloc(nlte_dimension); // the true population densities
    gsl_vector_memcpy(popvec, x);                          // are equal to the normed pops multiplied
    gsl_vector_mul(popvec, pop_norm_factor_vec);           // by the normalisation factors

    for (int row = 0; row < nlte_dimension; row++)
    {
      double recovered_balance_vector_elem = 0.;
      gsl_vector_view row_view = gsl_matrix_row(rate_matrix, row);
      gsl_blas_ddot(&row_view.vector, x, &recovered_balance_vector_elem);

      int ion, level;
      get_ion_level_of_nlte_vector_index(row,element,&ion,&level);

      // printout("index %4d (ion_stage %d level%4d): residual %+.2e recovered balance: %+.2e normed pop %.2e pop %.2e departure ratio %.4f\n",
      //          row,get_ionstage(element,ion),level, gsl_vector_get(residual_vector,row),
      //          recovered_balance_vector_elem, gsl_vector_get(x,row),
      //          gsl_vector_get(popvec, row),
      //          gsl_vector_get(x, row) / gsl_vector_get(x,get_nlte_vector_index(element,ion,0)));

      if (gsl_vector_get(popvec, row) < 0.0)
      {
        printout("WARNING: NLTE solver gave negative population to index %d (ion_stage %d level %d), pop = %g. Forcing departure coeff to 1.0\n",
                 row, get_ionstage(element,ion), level, gsl_vector_get(x,row) * gsl_vector_get(pop_norm_factor_vec,row));
        gsl_vector_set(popvec, row, 1.);
      }
    }

    gsl_vector_free(residual_vector);

    double *ion_populations = calloc(nions,sizeof(double));
    for (int ion = 0; ion < nions; ion++)
    {
      const int ion_stage = get_ionstage(element, ion);
      const int nlevels_nlte = get_nlevels_nlte(element, ion);
      const int index_gs = get_nlte_vector_index(element, ion, 0);
      printout("[ion_stage %d]\n", get_ionstage(element, ion));

      printout("  For ion_stage %d, the ground state populations are %g (function) and %g (matrix result with normed pop %g, ltepopnormfactor %g)\n",get_ionstage(element,ion),
               get_groundlevelpop(modelgridindex,element,ion),gsl_vector_get(popvec, index_gs),gsl_vector_get(x,index_gs),gsl_vector_get(pop_norm_factor_vec,index_gs));

      // store the NLTE level populations
      const int nlte_start = elements[element].ions[ion].first_nlte;
      double solution_ion_pop = 0.0;
      for (int level = 1; level <= nlevels_nlte; level++)
      {
        const int index = get_nlte_vector_index(element,ion,level);
        modelgrid[modelgridindex].nlte_pops[nlte_start + level - 1] = gsl_vector_get(popvec, index) / modelgrid[modelgridindex].rho;
        solution_ion_pop += gsl_vector_get(popvec, index);
      }

      // store the superlevel population if there is one
      if (nlevels_nlte != (get_nlevels(element,ion) - 1)) //a superlevel exists
      {
        const int index_sl = get_nlte_vector_index(element, ion, nlevels_nlte + 1);
        modelgrid[modelgridindex].nlte_pops[nlte_start+nlevels_nlte] = gsl_vector_get(popvec, index_sl) /
                                                                       modelgrid[modelgridindex].rho /
                                                                       superlevel_partfunc[ion];

        solution_ion_pop += gsl_vector_get(popvec, index_sl);
      }
      printout("  I had a ground level pop of %g, a part fn of %g and therefore an ion pop of %g\n",
               modelgrid[modelgridindex].composition[element].groundlevelpop[ion], modelgrid[modelgridindex].composition[element].partfunct[ion], modelgrid[modelgridindex].composition[element].partfunct[ion]*modelgrid[modelgridindex].composition[element].groundlevelpop[ion]/stat_weight(element,ion,0));

      //ionstagepop here must be called before setting the new ground level population
      printout("  For ion_stage %d the total population is %g, but was previously %g\n",
               ion_stage,solution_ion_pop,ionstagepop(modelgridindex,element,ion));

      // store the ground level population
      modelgrid[modelgridindex].composition[element].groundlevelpop[ion] = gsl_vector_get(popvec, index_gs);
      solution_ion_pop += gsl_vector_get(popvec, index_gs);

      precalculate_partfuncts(modelgridindex);

      ion_populations[ion] = solution_ion_pop;
      if (ion > 0)
      {
        double gspopratio = modelgrid[modelgridindex].composition[element].groundlevelpop[ion-1] / modelgrid[modelgridindex].composition[element].groundlevelpop[ion];

        double ionpot = epsilon(element,ion,0) - epsilon(element,ion-1,0);
        double T_e = get_Te(modelgridindex);
        double partfunct_ratio = modelgrid[modelgridindex].composition[element].partfunct[ion-1] / modelgrid[modelgridindex].composition[element].partfunct[ion];
        double gs_g_ratio = stat_weight(element,ion-1,0) / stat_weight(element,ion,0);
        double sbphi_gs = gs_g_ratio * SAHACONST * pow(T_e,-1.5) * exp(ionpot/KB/T_e) * get_nne(modelgridindex);
        double solution_ion_pop_ratio = ion_populations[ion-1] / ion_populations[ion];
        double sbphi = partfunct_ratio * SAHACONST * pow(T_e,-1.5) * exp(ionpot/KB/T_e) * get_nne(modelgridindex);

        printout("  The ratio of groundlevel pops (ion %d)/(ion %d) is %g, Saha-Boltzmann value is %g ratio %g\n",
                 get_ionstage(element,ion-1),ion_stage,gspopratio,sbphi_gs,gspopratio/sbphi_gs);
        printout("  The ratio of total pops (ion %d)/(ion %d) is %g, Saha-Boltzmann value is %g ratio %g\n",
                 get_ionstage(element,ion-1),ion_stage,solution_ion_pop_ratio,sbphi,solution_ion_pop_ratio/sbphi);
        const double nne = get_nne(modelgridindex);
        // calculate_partfunct(element, ion, modelgridindex) * modelgrid[modelgridindex].composition[element].groundlevelpop[ion] / stat_weight(element,ion,0))
        printout("  The corresponding phi-factor is: %g\n", (solution_ion_pop_ratio / nne));
        printout("  The ionization solver gives phi(ion_stage=%d / %d) = %g\n",
                 get_ionstage(element, ion-1), get_ionstage(element, ion), phi(element, ion-1, modelgridindex));
      }


      //double nne = get_nne(modelgridindex);
      //printout("From ion fract, the ion pop should be %g\n", ionfract(element, ion, modelgridindex, nne)*get_abundance(modelgridindex,element)/ elements[element].mass * get_rho(modelgridindex));
      //printout("I think that the element population is: %g (from abundance %g and rho %g)\n", get_abundance(modelgridindex,element)/elements[element].mass*get_rho(modelgridindex), get_abundance(modelgridindex,element), get_rho(modelgridindex));
      //printout("I currently think that the top ion is: %d\n", elements_uppermost_ion[tid][element]);
    }
    printout("The element population is: %g (from abundance) and %g (from matrix solution)\n",
             get_abundance(modelgridindex,element) / elements[element].mass * get_rho(modelgridindex), gsl_blas_dasum(popvec));

    if (atomic_number == 26 && timestep % 20 == 0)
    {
      // print_level_rates(modelgridindex, element, 0, 61, popvec, rate_matrix_rad_bb, rate_matrix_coll_bb, rate_matrix_rad_bf, rate_matrix_coll_bf, rate_matrix_ntcoll_bf);
      // print_level_rates(modelgridindex, element, 0, 62, popvec, rate_matrix_rad_bb, rate_matrix_coll_bb, rate_matrix_rad_bf, rate_matrix_coll_bf, rate_matrix_ntcoll_bf);
      print_level_rates(modelgridindex, element, 1, 35, popvec, rate_matrix_rad_bb, rate_matrix_coll_bb, rate_matrix_rad_bf, rate_matrix_coll_bf, rate_matrix_ntcoll_bf);
      print_level_rates(modelgridindex, element, 1, 75, popvec, rate_matrix_rad_bb, rate_matrix_coll_bb, rate_matrix_rad_bf, rate_matrix_coll_bf, rate_matrix_ntcoll_bf);
      print_level_rates(modelgridindex, element, 1, 76, popvec, rate_matrix_rad_bb, rate_matrix_coll_bb, rate_matrix_rad_bf, rate_matrix_coll_bf, rate_matrix_ntcoll_bf);
      // print_level_rates(modelgridindex, element, 2, 50, popvec, rate_matrix_rad_bb, rate_matrix_coll_bb, rate_matrix_rad_bf, rate_matrix_coll_bf, rate_matrix_ntcoll_bf);
      // print_level_rates(modelgridindex, element, 3, 50, popvec, rate_matrix_rad_bb, rate_matrix_coll_bb, rate_matrix_rad_bf, rate_matrix_coll_bf, rate_matrix_ntcoll_bf);
    }

    gsl_matrix_free(rate_matrix_rad_bb);
    gsl_matrix_free(rate_matrix_coll_bb);
    gsl_matrix_free(rate_matrix_rad_bf);
    gsl_matrix_free(rate_matrix_coll_bf);
    gsl_matrix_free(rate_matrix_ntcoll_bf);

    free(ion_populations);

    gsl_vector_free(x);
    gsl_vector_free(popvec);

    gsl_matrix_free(rate_matrix);
    gsl_matrix_free(rate_matrix_LU_decomp);
    gsl_vector_free(balance_vector);
    gsl_vector_free(pop_norm_factor_vec);
    free(superlevel_partfunc);
  }
  else
  {
    //abundance of this element is zero, so do not store any NLTE populations
    printout("Not solving for NLTE populations in cell %d at timestep %d for element %d due to zero abundance\n",
             modelgridindex,timestep,element);

    for (int ion = 0; ion < nions; ion++)
    {
      const int nlte_start = elements[element].ions[ion].first_nlte;
      const int nlevels_nlte = get_nlevels_nlte(element,ion);
      for (int level = 1; level < nlevels_nlte+1; level++)
      {
        modelgrid[modelgridindex].nlte_pops[nlte_start+level-1] = -1.0; // flag to indicate no useful data
      }
      if (nlevels_nlte != (get_nlevels(element,ion) - 1)) // a superlevel exists
      {
        modelgrid[modelgridindex].nlte_pops[nlte_start+nlevels_nlte] = -1.0;
      }
    }
  }
}


// this does single ion solving and will be deprecated at some point
double solve_nlte_pops(int element, int ion, int modelgridindex, int timestep)
//solves for nlte correction factors to level populations for levels
{
  int lower, level_use, lower_use, upper_use;

  double *rate_matrix;
  double *balance_vector;

  int ionisinglevels, ndowntrans, nuptrans, i;
  double epsilon_current;
  double R, C;
  double t_mid;
  double epsilon_target, epsilon_trans;
  int lineindex;
  double s_renorm;

  int nlte_size, nlte_start;

  double test_ratio;//,lag, check;
  double test_ratio_upper;

  int super_level;

  if (get_nlevels(element,ion) > 1)
  {
    //double nne = get_nne(modelgridindex);
    //double T_e = get_Te(modelgridindex);

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
        ndowntrans = elements[element].ions[ion].levels[level].downtrans[0].targetlevel;
        nuptrans = elements[element].ions[ion].levels[level].uptrans[0].targetlevel;

        // deexcitation
        for (i = 1; i <= ndowntrans; i++)
        {
          lower = elements[element].ions[ion].levels[level].downtrans[i].targetlevel;
          epsilon_target = elements[element].ions[ion].levels[level].downtrans[i].epsilon;
          //double statweight_target = elements[element].ions[ion].levels[level].downtrans[i].stat_weight;
          lineindex = elements[element].ions[ion].levels[level].downtrans[i].lineindex;
          epsilon_trans = epsilon_current - epsilon_target;

          mastate[tid].element = element;
          mastate[tid].ion = ion;
          mastate[tid].level = level;
          // mastate[tid].statweight = statweight;
          mastate[tid].nnlevel = 1.0;

          R = rad_deexcitation_ratecoeff(modelgridindex,element,ion,level,lower,epsilon_trans,lineindex,t_mid);
          C = col_deexcitation_ratecoeff(modelgridindex, epsilon_trans, lineindex);

          s_renorm = 1.0;

          if ((level == 0) || (is_nlte(element, ion, level) == 1))
          {
            level_use = level;
          }
          else
          {
            level_use = nlevels_nlte + 1;
            s_renorm = superlevel_boltzmann(modelgridindex,element,ion,level)/superlevel_partition;
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
        for (i = 1; i <= nuptrans; i++)
        {
          int upper = elements[element].ions[ion].levels[level].uptrans[i].targetlevel;
          epsilon_target = elements[element].ions[ion].levels[level].uptrans[i].epsilon;
          //double statweight_target = elements[element].ions[ion].levels[level].uptrans[i].stat_weight;
          lineindex = elements[element].ions[ion].levels[level].uptrans[i].lineindex;
          epsilon_trans = epsilon_target - epsilon_current;

          mastate[tid].element = element;
          mastate[tid].ion = ion;
          mastate[tid].level = level;
          // mastate[tid].statweight = statweight;
          mastate[tid].nnlevel = 1.0;

          R = rad_excitation_ratecoeff(modelgridindex, element, ion, level, upper, epsilon_trans, lineindex, t_mid);
          C = col_excitation_ratecoeff(modelgridindex,lineindex,epsilon_trans);

          if ((level == 0) || (is_nlte(element, ion, level) == 1))
          {
            level_use = level;
            s_renorm = 1.0;
          }
          else
          {
            level_use = nlevels_nlte+1;
            s_renorm = superlevel_boltzmann(modelgridindex,element,ion,level)/superlevel_partition;
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
        ionisinglevels = get_bfcontinua(element,ion);
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
          // mastate[tid].statweight = statweight;
          mastate[tid].nnlevel = 1.0;
          for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion,level); phixstargetindex++)
          {
            int upper = get_phixsupperlevel(element,ion,level,phixstargetindex);
            epsilon_trans = epsilon(element,ion+1,upper) - epsilon_current;
            R = get_corrphotoioncoeff(element,ion,level,phixstargetindex,modelgridindex);
            C = col_ionization_ratecoeff(modelgridindex, element, ion, level, phixstargetindex, epsilon_trans);

            upper_use = nlte_size - 1; //ion above

            rate_matrix[level_use*(nlte_size) + level_use] -= (R + C) * s_renorm;
            rate_matrix[upper_use*(nlte_size) + level_use] += (R + C) * s_renorm;
          }

          // recombination
          mastate[tid].element = element;
          mastate[tid].ion = ion+1;
          mastate[tid].nnlevel = 1.0;
          for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion,level); phixstargetindex++)
          {
            int upper = get_phixsupperlevel(element,ion,level,phixstargetindex);
            mastate[tid].level = upper;
            // mastate[tid].statweight = stat_weight(element,ion+1,upper);
            epsilon_trans = epsilon(element,ion+1,upper) - epsilon_current;
            R = rad_recombination_ratecoeff(modelgridindex, element, ion+1, upper, level);
            //printout("rad recombination of element %d, ion %d, level %d, to lower level %d has rate %g (ne %g and Te %g)\n",element,ion,mastate[tid].level,level,R/nne,nne,T_e);
            //printout("%d %d %d %d %g %g %g \n",element,ion,mastate[tid].level,level,R/nne,nne,T_e);
            C = col_recombination_ratecoeff(modelgridindex,element,ion+1,upper,level,epsilon_trans);
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
        modelgrid[modelgridindex].nlte_pops[nlte_start+level-1] = gsl_vector_get(x,level)/modelgrid[modelgridindex].rho;
        //printout("I have interfered with index %d.\n", nlte_start+level-1);
        //modelgrid[modelgridindex].nlte_pops[nlte_start+level-1] = ((lag*modelgrid[modelgridindex].nlte_pops[nlte_start+level-1]) + gsl_vector_get(x,level))/(lag + 1.0)/modelgrid[modelgridindex].rho;
      }
      // If there is a superlevel then write that too

      if (super_level == 1)
      {
        //printout("I thought the super level was: %g\n", modelgrid[modelgridindex].nlte_pops[nlte_start+nlevels_nlte]);

        modelgrid[modelgridindex].nlte_pops[nlte_start+nlevels_nlte] = gsl_vector_get(x,nlevels_nlte+1)/modelgrid[modelgridindex].rho/superlevel_partition;
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
      for (int level = 1; level < nlevels_nlte+1; level++)
      {
        modelgrid[modelgridindex].nlte_pops[nlte_start+level-1] = -1.0;///flag to indicate no useful data
      }
      if (super_level == 1)
      {
        modelgrid[modelgridindex].nlte_pops[nlte_start+nlevels_nlte] = -1.0;
      }
      test_ratio = 0.0;

    }
    return(test_ratio);
  }
  else
  {
    return 0; //Case for ion with only one level
  }
}


double superlevel_boltzmann(int modelgridindex, int element, int ion, int level)
{
  double T_exc = get_TJ(modelgridindex);
  double E_level = epsilon(element,ion,level);
  double E_superlevel = epsilon(element,ion,get_nlevels_nlte(element,ion)+1);

  return stat_weight(element,ion,level) * exp(-(E_level-E_superlevel)/KB/T_exc);
}


void nltepop_open_file(int my_rank)
{
  char filename[100];
  sprintf(filename,"nlte_%.4d.out",my_rank);
  if ((nlte_file = fopen(filename, "w")) == NULL)
  {
    printout("Cannot open %s.\n",filename);
    abort();
  }
}


void nltepop_write_to_file(int n, int timestep)
{
  if (initial_iteration) // NLTE solver hasn't been run yet
    return;

  fprintf(nlte_file,"timestep %d modelgridindex %d T_R %g T_e %g W %g T_J %g nne %g\n",timestep,n,get_TR(n),get_Te(n),get_W(n),get_TJ(n),get_nne(n));
  for (int nlte = 0; nlte < total_nlte_levels; nlte++)
  {
    int element = -1;
    int ion = -1;
    int level = -1;
    for (int dum_element = 0; dum_element < nelements; dum_element++)
    {
      const int nions = get_nions(dum_element);
      for (int dum_ion = 0; dum_ion < nions; dum_ion++)
      {
        const int ion_first_nlte = elements[dum_element].ions[dum_ion].first_nlte;
        //printout("NLTETEST: element %d ion %d first_nlte %d looking for %d\n",dum_element,dum_ion,ion_first_nlte,nlte);
        if (nlte >= ion_first_nlte)
        {
          element = dum_element;
          ion = dum_ion;
          level = nlte - ion_first_nlte + 1;
        }
      }
    }
    if (level == -1)
    {
      printout("FATAL: matching ion & level for NLTE index %d not found",nlte);
      abort();
    }
    else
    {
      //printout("found it in element %d ion %d level %d\n",element,ion,level);
    }
    double nnlevelnlte = modelgrid[n].nlte_pops[nlte] * modelgrid[n].rho;
    double E_level = epsilon(element,ion,level);
    double E_ground = epsilon(element,ion,0);
    double T_exc = get_TJ(n); //TODO: is this right? should be TJ?
    //double W = get_W(n);
    double W = 1.;
    double nnlevellte = get_groundlevelpop(n,element,ion) * W *
        stat_weight(element,ion,level)/stat_weight(element,ion,0) *
        exp(-(E_level-E_ground)/KB/T_exc);
    //double nnlevellte = calculate_levelpop_lte(n,element,ion,level); this function includes a minpop
    //if (ion == 1)
    {
      if (level == 1)
        fprintf(nlte_file,"nlte_index - element %d ion %d level 0 energy 0 nlte_pop - nnlevel_NLTE - nnlevel_LTE %g\n",
              element,ion,get_groundlevelpop(n,element,ion));

      fprintf(nlte_file,"nlte_index %d element %d ion %d ",nlte,element,ion);
      if (is_nlte(element, ion, level))
        fprintf(nlte_file,"level %d ",level);
      else
        fprintf(nlte_file,"level SL ");

      fprintf(nlte_file,"energy %g nlte_pop %g nnlevel_NLTE %g nnlevel_LTE %g\n",
              E_level-E_ground,modelgrid[n].nlte_pops[nlte],nnlevelnlte,nnlevellte);
    }
  }
  fprintf(nlte_file,"\n");
  //printout("I just wrote %g (really %g\n", modelgrid[0].nlte_pops[820] , modelgrid[0].nlte_pops[820]*modelgrid[0].rho);

  //fprintf(nlte_file,"%d %g %g %g %g ",n,get_TR(n),get_Te(n),get_W(n),get_TJ(n));
  // for (int dummy_element = 0; dummy_element < nelements; dummy_element++)
  // {
  //   int nions = get_nions(dummy_element);
  //   for (int dummy_ion = 0; dummy_ion < nions; dummy_ion++)
  //   {
  //     int nlte = get_nlevels_nlte(dummy_element, dummy_ion);
  //     if (nlte > 0)
  //     {
  //       for (int dummy_level = 1; dummy_level < (nlte+1); dummy_level++)
  //       {
  //         fprintf(nlte_file,"%g ", calculate_exclevelpop_old(n, dummy_element, dummy_ion, dummy_level));
  //       }
  //     }
  //   }
  // }
  //fprintf(nlte_file,"\n");
  fflush(nlte_file);
}