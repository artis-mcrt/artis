#include "sn3d.h"
#include "atomic.h"
#include "grid_init.h"
#include "ltepop.h"
#include "macroatom.h"
#include "nltepop.h"
#include "update_grid.h"

void nlte_pops_element(int element, int modelgridindex, int timestep)
//solves for nlte correction factors to level populations for levels in all ions of an element
{
  double t_mid = time_step[timestep].mid;

  printout("Solving for NLTE populations in cell %d at timestep %d for element %d\n",modelgridindex,timestep,element);
  printout("T_E %g T_R was %g, setting to 3000 \n",get_Te(modelgridindex),get_TR(modelgridindex));
  set_TR(modelgridindex,3000);

  int nions = get_nions(element);
  //int nions = 3; //TODO: remove, testing

  int nlte_dimension = 0;
  double *superlevel_partfunc = calloc(nions,sizeof(double)); //space is allocated for every ion, whether or not it has a superlevel
  for (int ion = 0; ion < nions; ion++)
  {
    int nlevels_nlte = get_nlevels_nlte(element,ion);

    //this is the total number of nlte_levels (i.e. the size of the
    //storage). Our rate matrix will need to be of this dimension +2: the
    //ground state, the "super level".
    //If there's no super level needed then we only need +1
    nlte_dimension += nlevels_nlte + 1;
    printout("NLTE: setting up ion %d, which contributes %d to the vector dimension ",ion,nlevels_nlte+1);
    if (nlevels_nlte == (get_nlevels(element, ion) - 1))
    {
      printout("(no super level)\n");
      superlevel_partfunc[ion] = 0.0;
    }
    else
    {
      nlte_dimension += 1;
      superlevel_partfunc[ion] = 0.0;
      for (int level = 1; level < get_nlevels(element,ion); level++)
      {
        if (is_nlte(element, ion, level) == 0)
          superlevel_partfunc[ion] += superlevel_boltzmann(modelgridindex,element,ion,level);
      }
      printout("(plus a superlevel with partfunc %g)\n",superlevel_partfunc[ion]);
    }
  }

  if (get_abundance(modelgridindex, element) > 0.0)
  {
    double *rate_matrix = calloc(nlte_dimension*nlte_dimension,sizeof(double));
    double *balance_vector = calloc(nlte_dimension,sizeof(double));

    if (!rate_matrix)
    {
      printout("Cannot allocate NLTE rate matrix memory.\n");
    }

    if (!balance_vector)
    {
      printout("Cannot allocate NLTE vector memory.\n");
    }

    printout("NLTE: the vector dimension is %d.\n", nlte_dimension);

    for (int ion = 0; ion < nions; ion++)
    {
      for (int level = 0; level < get_nlevels(element,ion); level++)
      {
        double statweight = stat_weight(element,ion,level);
        double epsilon_current = epsilon(element,ion,level);
        int ndowntrans = elements[element].ions[ion].levels[level].downtrans[0].targetlevel;
        int nuptrans = elements[element].ions[ion].levels[level].uptrans[0].targetlevel;

        //if (level < 10)
        //  printout("NLTE: the index of ion %d level %d is %d\n",ion,level,get_nlte_vector_index(element,ion,level));

        // deexcitation
        for (int i = 1; i <= ndowntrans; i++)
        {
          int lower = elements[element].ions[ion].levels[level].downtrans[i].targetlevel;
          double epsilon_target = elements[element].ions[ion].levels[level].downtrans[i].epsilon;
          int lineindex = elements[element].ions[ion].levels[level].downtrans[i].lineindex;
          double epsilon_trans = epsilon_current - epsilon_target;

          mastate[tid].element = element;
          mastate[tid].ion = ion;
          mastate[tid].level = level;
          mastate[tid].statweight = statweight;
          mastate[tid].nnlevel = 1.0;

          double R = nlte_matrix_rad_deexcitation(modelgridindex,lower,epsilon_trans,lineindex,t_mid);
          double C = col_deexcitation(modelgridindex,lower,epsilon_trans,lineindex);
          //double C = 0.0; //TODO: remove, testing only

          //TOOD: remove
          //if (level < 10 && lower < 10)
          //  printout("deexc: level %d, lower %d, R/C %g\n",level,lower,R/C);

          int upper_index = get_nlte_vector_index(element,ion,level);
          int lower_index = get_nlte_vector_index(element,ion,lower);

          double s_renorm = 1.0;

          if ((level != 0) && (is_nlte(element,ion,level) != 1))
            s_renorm = superlevel_boltzmann(modelgridindex,element,ion,level)/superlevel_partfunc[ion];

          rate_matrix[upper_index * nlte_dimension + upper_index] -= (R + C) * s_renorm;
          rate_matrix[lower_index * nlte_dimension + upper_index] += (R + C) * s_renorm;
          if (R + C < 0)
            printout("WARNING: Negative de-excitation rate from ion %d level %d to level %d\n",ion,level,lower);
        }

        // excitation
        for (int i = 1; i <= nuptrans; i++)
        {
          int upper = elements[element].ions[ion].levels[level].uptrans[i].targetlevel;
          double epsilon_target = elements[element].ions[ion].levels[level].uptrans[i].epsilon;
          int lineindex = elements[element].ions[ion].levels[level].uptrans[i].lineindex;
          double epsilon_trans = epsilon_target - epsilon_current;

          mastate[tid].element = element;
          mastate[tid].ion = ion;
          mastate[tid].level = level;
          mastate[tid].statweight = statweight;
          mastate[tid].nnlevel = 1.0;

          double R = nlte_matrix_rad_excitation(modelgridindex,upper,epsilon_trans,lineindex,t_mid); //CAUSES massive departures from LTE, probably broken
          double C = col_excitation(modelgridindex,upper,lineindex,epsilon_trans);
          //double C = 0.0; //TODO: remove, testing only

          int lower_index = get_nlte_vector_index(element,ion,level);
          int upper_index = get_nlte_vector_index(element,ion,upper);

          double s_renorm = 1.0;

          if ((level != 0) && (is_nlte(element,ion,level) != 1))
            s_renorm = superlevel_boltzmann(modelgridindex,element,ion,level)/superlevel_partfunc[ion];

          rate_matrix[lower_index * nlte_dimension + lower_index] -= (R + C) * s_renorm;
          rate_matrix[upper_index * nlte_dimension + lower_index] += (R + C) * s_renorm;
          if (R + C < 0)
            printout("WARNING: Negative excitation rate from ion %d level %d tolevel %d\n",ion,level,upper);
        }

        #ifdef NT_ON
          if (ion < nions-1)
          {
            double Y = nt_ionization_rate(modelgridindex,element,ion);
            Y = 0.0; // TODO: remove, testing only

            int lower_index = get_nlte_vector_index(element,ion,level);
            int upper_index = get_nlte_vector_index(element,ion+1,0);

            double s_renorm = 1.0;

            if ((level != 0) && (is_nlte(element,ion,level) != 1))
            {
              s_renorm = superlevel_boltzmann(modelgridindex,element,ion,level) / superlevel_partfunc[ion];
            }

            rate_matrix[lower_index * nlte_dimension + lower_index] -= Y * s_renorm;
            rate_matrix[upper_index * nlte_dimension + lower_index] += Y * s_renorm;
            if (Y < 0)
              printout("WARNING: Negative NT_ionization rate from ion %d level %d\n",ion,level);
          }
        #endif

        //now put in the photoionization/recombination processes
        int ionisinglevels = get_bfcontinua(element,ion);
        if (ion < nions-1 && level < ionisinglevels)  //&& get_ionstage(element,ion) < get_element(element)+1)
        {
          // ionization
          int lower_index = get_nlte_vector_index(element,ion,level);

          mastate[tid].element = element;
          mastate[tid].ion = ion;
          mastate[tid].level = level;
          mastate[tid].statweight = statweight;
          mastate[tid].nnlevel = 1.0;
          for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion,level); phixstargetindex++)
          {
            int upper = get_phixsupperlevel(element,ion,level,phixstargetindex);
            int upper_index = get_nlte_vector_index(element,ion+1,upper);
            double epsilon_trans = epsilon(element,ion+1,upper) - epsilon_current;
            //double R = photoionization(modelgridindex,phixstargetindex,epsilon_trans);
            double C = col_ionization(modelgridindex,phixstargetindex,epsilon_trans);
            double R = 0.0; //TODO: remove, testing only

            double s_renorm = 1.0;
            if ((level != 0) && (is_nlte(element,ion,level) != 1))
              s_renorm = superlevel_boltzmann(modelgridindex,element,ion,level)/superlevel_partfunc[ion];

            rate_matrix[lower_index * nlte_dimension + lower_index] -= (R + C) * s_renorm;
            rate_matrix[upper_index * nlte_dimension + lower_index] += (R + C) * s_renorm;
            if (R + C < 0)
              printout("WARNING: Negative ionization rate from ion %d level %d phixstargetindex %d\n",ion,level,phixstargetindex);
          }

          // recombination
          mastate[tid].element = element;
          mastate[tid].ion = ion+1;
          mastate[tid].nnlevel = 1.0;
          for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion,level); phixstargetindex++)
          {
            int upper = get_phixsupperlevel(element,ion,level,phixstargetindex);
            int upper_index = get_nlte_vector_index(element,ion+1,upper);
            mastate[tid].level = upper;
            mastate[tid].statweight = stat_weight(element,ion+1,upper);
            double epsilon_trans = epsilon(element,ion+1,upper) - epsilon_current;
            //double R = rad_recombination(modelgridindex,level,epsilon_trans);
            double C = col_recombination(modelgridindex,level,epsilon_trans);
            double R = 0.0; //TODO: remove, testing only

            //printout("rad recombination of element %d, ion %d, level %d, to lower level %d has rate %g (ne %g and Te %g)\n",element,ion,mastate[tid].level,level,R/nne,nne,T_e);
            //printout("%d %d %d %d %g %g %g \n",element,ion,mastate[tid].level,level,R/nne,nne,T_e);

            double s_renorm = 1.0;
            if ((upper != 0) && (is_nlte(element,ion+1,upper) != 1))
              s_renorm = superlevel_boltzmann(modelgridindex,element,ion+1,upper)/superlevel_partfunc[ion+1];

            rate_matrix[upper_index * nlte_dimension + upper_index] -= (R + C) * s_renorm;
            rate_matrix[lower_index * nlte_dimension + upper_index] += (R + C) * s_renorm;
            if (R + C < 0)
              printout("WARNING: Negative recombination rate to ion %d level %d phixstargetindex %d\n",ion,level,phixstargetindex);
          }

        } // if (level is an ionising level)

      } // level loop
    } // ion loop

    // replace the first row of the matrix and balance vector with the normalisation constraint
    // and apply normalisation factors
    double *pop_norm_factors = calloc(nlte_dimension,sizeof(double));
    for (int column = 0; column < nlte_dimension; column++)
    {
      int ion, level;
      get_ion_level_of_nlte_vector_index(column,element,&ion,&level);

      pop_norm_factors[column] = calculate_levelpop_lte(modelgridindex,element,ion,level);

      if ((level != 0) && (is_nlte(element,ion,level) != 1))
      { //level is a superlevel, so add populations of higher levels to the norm factor
        for (int dummylevel = level+1; dummylevel < get_nlevels(element,ion); dummylevel++)
        {
          if (is_nlte(element,ion,dummylevel) != 1)
            pop_norm_factors[column] += calculate_levelpop_lte(modelgridindex,element,ion,dummylevel);
        }
        //NOTE: above calulation is not always the sum of LTE populations,
        //since calculate_levelpop_lte imposes MINPOP minimum
        printout("superlevel norm factor index %d is %g, partfunc is %g, partfunc*SL_levelpop %g\n",
                 column,pop_norm_factors[column],superlevel_partfunc[ion],
                 superlevel_partfunc[ion]*calculate_levelpop_lte(modelgridindex,element,ion,level)/stat_weight(element,ion,level));
      }
      //printout("LEVELFIND: for index %d, the ion is %d, level is %d, normalising factor (lte_pop) is %g\n",column,ion,level,pop_norm_factors[column]);

      rate_matrix[column] = pop_norm_factors[column];
      balance_vector[column] = 0.0;

      for (int row = 1; row < nlte_dimension; row++)
      {
        rate_matrix[row * nlte_dimension + column] *= pop_norm_factors[column];
      }
    }

    //set first vector entry to the element population
    balance_vector[0] = get_abundance(modelgridindex,element) / elements[element].mass * get_rho(modelgridindex);

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
    // their interactions and setting their normalised population (probably departure coeff) to 1.0
    filter_nlte_matrix(element,nlte_dimension,rate_matrix,balance_vector);

    //make a copy of the rate matrix, since gsl routines will alter the original
    double *rate_matrix_copy = calloc(nlte_dimension*nlte_dimension,sizeof(double));
    for (int row = 0; row < nlte_dimension; row++)
    {
      for (int column = 0; column < nlte_dimension; column++)
        rate_matrix_copy[row * nlte_dimension + column] = rate_matrix[row * nlte_dimension + column];
    }

    gsl_matrix_view m = gsl_matrix_view_array(rate_matrix, nlte_dimension, nlte_dimension);
    gsl_permutation *p = gsl_permutation_alloc(nlte_dimension);

    int s; //sign of the transformation
    gsl_linalg_LU_decomp(&m.matrix, p, &s);

    gsl_vector_view b = gsl_vector_view_array(balance_vector, nlte_dimension);
    gsl_vector *x = gsl_vector_alloc(nlte_dimension);
    gsl_linalg_LU_solve(&m.matrix, p, &b.vector, x); //solve matrix equation m * x = b for x (populations)
    //gsl_linalg_HH_solve (&m.matrix, &b.vector, x);

    //printf("after solving: rate %p balance %p NULL %p\n", rate_matrix, balance_vector, NULL);
    for (int row = 0; row < nlte_dimension; row++)
    {
      double sum = 0.;
      for (int column = 0; column < nlte_dimension; column++)
      {
        sum += rate_matrix_copy[row * nlte_dimension + column] * gsl_vector_get(x,column);
      }

      int ion, level;
      get_ion_level_of_nlte_vector_index(row,element,&ion,&level);

      printout("row %4d (ion %d level%4d): Recovered balance: %+.2e, Pop vector: value %+.2e pop %.2e depature ratio %.2e\n",
               row,ion,level,sum,gsl_vector_get(x,row),gsl_vector_get(x,row)*pop_norm_factors[row],gsl_vector_get(x,row)/gsl_vector_get(x,get_nlte_vector_index(element,ion,0)));
      if (gsl_vector_get(x,row)*pop_norm_factors[row] < 0.0)
      {
        printout("WARNING: NLTE solver gave negative population to index %d (ion %d level %d), pop = %g\n",row,ion,level,gsl_vector_get(x,row)*pop_norm_factors[row]);
      }
    }

    double *ion_populations = calloc(nions,sizeof(double));
    double elempopmatrix = 0.0;
    for (int ion = 0; ion < nions; ion++)
    {
      int nlevels_nlte = get_nlevels_nlte(element,ion);
      int index_gs = get_nlte_vector_index(element,ion,0);
      printout("For ion %d, the ground state populations are %g (function) and %g (matrix result with departure coeff %g, ltepopnormfactor %g)\n",ion,
        get_groundlevelpop(modelgridindex,element,ion),pop_norm_factors[index_gs]*gsl_vector_get(x,index_gs),gsl_vector_get(x,index_gs),pop_norm_factors[index_gs]);

      // store the NLTE level populations
      int nlte_start = elements[element].ions[ion].first_nlte;
      double ionstagepopmatrix = 0.0;
      ionstagepopmatrix += gsl_vector_get(x,get_nlte_vector_index(element,ion,0));
      for (int level = 1; level <= nlevels_nlte; level++)
      {
        int index = get_nlte_vector_index(element,ion,level);
        modelgrid[modelgridindex].nlte_pops[nlte_start+level-1] = pop_norm_factors[index]*gsl_vector_get(x,index)/modelgrid[modelgridindex].rho;
        ionstagepopmatrix += gsl_vector_get(x,index) * pop_norm_factors[index];
        //printout("I have interfered with index %d.\n", nlte_start+level-1);
        //modelgrid[modelgridindex].nlte_pops[nlte_start+level-1] = ((lag*modelgrid[modelgridindex].nlte_pops[nlte_start+level-1]) + gsl_vector_get(x,level))/(lag + 1.0)/modelgrid[modelgridindex].rho;
      }

      // store the superlevel population if there is one
      if (nlevels_nlte != (get_nlevels(element,ion) - 1)) //a superlevel exists
      {
        int index_sl = get_nlte_vector_index(element,ion,nlevels_nlte+1);
        modelgrid[modelgridindex].nlte_pops[nlte_start+nlevels_nlte] = pop_norm_factors[index_sl]*gsl_vector_get(x,index_sl)/modelgrid[modelgridindex].rho/superlevel_partfunc[ion];

        ionstagepopmatrix += gsl_vector_get(x,index_sl) * pop_norm_factors[index_sl];
      }
      //printout("I had a ground level pop of %g, a part fn of %g and therefore an ion pop of %g\n", modelgrid[modelgridindex].composition[element].groundlevelpop[ion], modelgrid[modelgridindex].composition[element].partfunct[ion], modelgrid[modelgridindex].composition[element].partfunct[ion]*modelgrid[modelgridindex].composition[element].groundlevelpop[ion]/stat_weight(element,ion,0));

      // store the ground level population
      modelgrid[modelgridindex].composition[element].groundlevelpop[ion] = gsl_vector_get(x,index_gs) * pop_norm_factors[index_gs];
      ionstagepopmatrix += modelgrid[modelgridindex].composition[element].groundlevelpop[ion];

      ion_populations[ion] = ionstagepopmatrix;
      if (ion > 0)
      {
        double gspopratio = modelgrid[modelgridindex].composition[element].groundlevelpop[ion-1] / modelgrid[modelgridindex].composition[element].groundlevelpop[ion];

        double ionpot = epsilon(element,ion,0) - epsilon(element,ion-1,0);
        float T_e = get_Te(modelgridindex);
        double partfunct_ratio = modelgrid[modelgridindex].composition[element].partfunct[ion-1]/modelgrid[modelgridindex].composition[element].partfunct[ion];
        double gs_g_ratio = stat_weight(element,ion-1,0)/stat_weight(element,ion,0);
        double sbphi_gs = gs_g_ratio * SAHACONST * pow(T_e,-1.5) * exp(ionpot/KB/T_e) * get_nne(modelgridindex);
        double matrix_ion_pop_ratio = ion_populations[ion-1]/ion_populations[ion];
        double sbphi = partfunct_ratio * SAHACONST * pow(T_e,-1.5) * exp(ionpot/KB/T_e) * get_nne(modelgridindex);

        printout("IONBALANCE: the ratio of ion %d to ion %d groundlevel pops is %g, Saha-Boltzman value is %g ratio %g\n",ion-1,ion,gspopratio,sbphi_gs,gspopratio/sbphi_gs);
        printout("IONBALANCE: the ratio of ion %d to ion %d total pops is %g, Saha-Boltzman value is %g ratio %g\n",ion-1,ion,matrix_ion_pop_ratio,sbphi,matrix_ion_pop_ratio/sbphi);
      }

      printout("NLTE: for ion %d the total population is %g, but was previously %g\n",ion,ionstagepopmatrix,ionstagepop(modelgridindex,element,ion));

      //printout("Accordingly, my phi-factor is: %g\n", calculate_partfunct(element, ion, modelgridindex)*modelgrid[modelgridindex].composition[element].groundlevelpop[ion]/stat_weight(element,ion,0))/nne);
      //printout("The ionization solver gives phi = %g\n", phi(element, ion, modelgridindex));

      //double nne = get_nne(modelgridindex);
      //printout("From ion fract, the ion pop should be %g\n", ionfract(element, ion, modelgridindex, nne)*get_abundance(modelgridindex,element)/ elements[element].mass * get_rho(modelgridindex));
      //printout("I think that the element population is: %g (from abundance %g and rho %g)\n", get_abundance(modelgridindex,element)/elements[element].mass*get_rho(modelgridindex), get_abundance(modelgridindex,element), get_rho(modelgridindex));
      //printout("I currently think that the top ion is: %d\n", elements_uppermost_ion[tid][element]);
      elempopmatrix += ionstagepopmatrix;
    }
    printout("The element population is: %g (from abundance) and %g (from matrix solution)\n", get_abundance(modelgridindex,element)/elements[element].mass*get_rho(modelgridindex), elempopmatrix);

    free(ion_populations);

    gsl_permutation_free(p);
    gsl_vector_free(x);

    free(superlevel_partfunc);
    free(rate_matrix);
    free(rate_matrix_copy);
    free(balance_vector);
    free(pop_norm_factors);
  }
  else
  {
    //abundance of this element is zero, so do not store any NLTE populations
    for (int ion = 0; ion < nions; ion++)
    {
      int nlte_start = elements[element].ions[ion].first_nlte;
      int nlevels_nlte = get_nlevels_nlte(element,ion);
      for (int level = 1; level < nlevels_nlte+1; level++)
      {
        modelgrid[modelgridindex].nlte_pops[nlte_start+level-1] = -1.0;///flag to indicate no useful data
      }
      if (nlevels_nlte != (get_nlevels(element,ion) - 1)) //a superlevel exists
      {
        modelgrid[modelgridindex].nlte_pops[nlte_start+nlevels_nlte] = -1.0;
      }
    }
  }
}

int get_nlte_vector_index(int element_in, int ion_in, int level_in)
{
  int index = 0;
  for (int ion = 0; ion < ion_in; ion++)
  {
    int nlevels_nlte = get_nlevels_nlte(element_in,ion);

    index += nlevels_nlte + 1;
    if (nlevels_nlte != (get_nlevels(element_in,ion) - 1))
      index += 1; // there's a superlevel here
  }
  int nlevels_nlte = get_nlevels_nlte(element_in,ion_in);
  if (is_nlte(element_in, ion_in, level_in) == 1)
    index += level_in;
  else
    index += nlevels_nlte + 1; //the index of the superlevel

  return index;
}


void get_ion_level_of_nlte_vector_index(int index, int element, int *ion, int *level)
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


void filter_nlte_matrix(int element, int nlte_dimension, double *rate_matrix, double *balance_vector)
// find rows and columns that barely interaction with other levels, and effectively
// removing them by zeroing their interactions and setting their departure
// coeff to 1.0
{
  for (int index = 0; index < nlte_dimension; index++)
  {
    double row_max = 0.0;
    for (int column = 0; column < nlte_dimension; column++)
    {
      double element_value = rate_matrix[index * nlte_dimension + column];
      if (element_value > row_max)
        row_max = element_value;
    }
    double col_max = 0.0;
    for (int row = 0; row < nlte_dimension; row++)
    {
      double element_value = rate_matrix[row * nlte_dimension + index];
      if (element_value > col_max)
        col_max = element_value;
    }
    if ((row_max < MINPOP) || (col_max < MINPOP))
    {
      int ion = -1;
      int level = -1;
      get_ion_level_of_nlte_vector_index(index,element,&ion,&level);
      printout("Eliminating row and col with index %d (ion %d level %d) row_max %g col_max %g\n",index,ion,level,row_max,col_max);
      if (level == 0)
      {
        printout("(Except that it's a ground state, so leaving as is.)\n");
      }
      else
      {
        double gs_index = get_nlte_vector_index(element,ion,0);
        eliminate_nlte_matrix_rowcol(index,gs_index,nlte_dimension,rate_matrix,balance_vector);
      }
    }
  }
}


void eliminate_nlte_matrix_rowcol(int index, int gs_index, int nlte_dimension, double *rate_matrix, double *balance_vector)
{
  for (int column = 0; column < nlte_dimension; column++)
    rate_matrix[index * nlte_dimension + column] = 0.0;

  for (int row = 0; row < nlte_dimension; row++)
    rate_matrix[row * nlte_dimension + index] = 0.0;

  rate_matrix[index * nlte_dimension + gs_index] = -1.0;
  rate_matrix[index * nlte_dimension + index] = 1.0;
  balance_vector[index] = 0.0;
}


// this does single ion solving and will be deprecated at some point
double nlte_pops(int element, int ion, int modelgridindex, int timestep)
//solves for nlte correction factors to level populations for levels
{
  int lower, level_use, lower_use, upper_use;

  double *rate_matrix;
  double *balance_vector;

  int ionisinglevels, ndowntrans, nuptrans, i;
  double statweight, epsilon_current;
  double R, C, Y;
  double t_mid;
  double epsilon_target, statweight_target, epsilon_trans;
  int lineindex;
  double s_renorm;
  PKT dummy;

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

    dummy.where = modelgridindex;

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

      /*
      if (superlevel_partition > 0.0)
      {
        printout("I found a superlevel and have computed a partition function for its substates of %g.\n", superlevel_partition);
      }
      else
      {
        printout("I don't know about any super level for this case.\n");
      }
      */

      for (int level = 0; level < get_nlevels(element,ion); level++)
      {
        statweight = stat_weight(element,ion,level);
        epsilon_current = epsilon(element,ion,level);
        ndowntrans = elements[element].ions[ion].levels[level].downtrans[0].targetlevel;
        nuptrans = elements[element].ions[ion].levels[level].uptrans[0].targetlevel;

        // deexcitation
        for (i = 1; i <= ndowntrans; i++)
        {
          lower = elements[element].ions[ion].levels[level].downtrans[i].targetlevel;
          epsilon_target = elements[element].ions[ion].levels[level].downtrans[i].epsilon;
          statweight_target = elements[element].ions[ion].levels[level].downtrans[i].stat_weight;
          lineindex = elements[element].ions[ion].levels[level].downtrans[i].lineindex;
          epsilon_trans = epsilon_current - epsilon_target;

          mastate[tid].element = element;
          mastate[tid].ion = ion;
          mastate[tid].level = level;
          mastate[tid].statweight = statweight;
          mastate[tid].nnlevel = 1.0;

          R = rad_deexcitation(modelgridindex,lower,epsilon_trans,lineindex,t_mid);
          C = col_deexcitation(modelgridindex,lower,epsilon_trans,lineindex);

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
          statweight_target = elements[element].ions[ion].levels[level].uptrans[i].stat_weight;
          lineindex = elements[element].ions[ion].levels[level].uptrans[i].lineindex;
          epsilon_trans = epsilon_target - epsilon_current;

          mastate[tid].element = element;
          mastate[tid].ion = ion;
          mastate[tid].level = level;
          mastate[tid].statweight = statweight;
          mastate[tid].nnlevel = 1.0;

          R = rad_excitation(modelgridindex,upper,epsilon_trans,lineindex,t_mid);//,T_R,W);
          C = col_excitation(modelgridindex,upper,lineindex,epsilon_trans);

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

        #ifdef NT_ON
          if (ion < get_nions(element)-1)
          {
            Y = nt_ionization_rate(modelgridindex,element,ion);

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
        #endif

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
          mastate[tid].statweight = statweight;
          mastate[tid].nnlevel = 1.0;
          for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion,level); phixstargetindex++)
          {
            int upper = get_phixsupperlevel(element,ion,level,phixstargetindex);
            epsilon_trans = epsilon(element,ion+1,upper) - epsilon_current;
            R = photoionization(modelgridindex,phixstargetindex,epsilon_trans);
            C = col_ionization(modelgridindex,phixstargetindex,epsilon_trans);

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
            mastate[tid].statweight = stat_weight(element,ion+1,upper);
            epsilon_trans = epsilon(element,ion+1,upper) - epsilon_current;
            R = rad_recombination(modelgridindex,level,epsilon_trans);
            //printout("rad recombination of element %d, ion %d, level %d, to lower level %d has rate %g (ne %g and Te %g)\n",element,ion,mastate[tid].level,level,R/nne,nne,T_e);
            //printout("%d %d %d %d %g %g %g \n",element,ion,mastate[tid].level,level,R/nne,nne,T_e);
            C = col_recombination(modelgridindex,level,epsilon_trans);
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

      /*printout("Rate matrix | balance vector:\n");
      for (int row = 0; row < nlte_size; row++)
      {
        for (int column = 0; column < nlte_size; column++)
        {
          char str[15];
          sprintf(str, "%+.2e ", rate_matrix[row * nlte_size + column]);
          printout(str);
        }
        printout("| ");
        char str[15];
        sprintf(str, "%+.3e\n", balance_vector[row]);
        printout(str);
      }
      printout("\n");*/

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
            exit(0);
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
          exit(0);
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


      if ((get_groundlevelpop(modelgridindex,element,ion) > (1.1*MINPOP)) && (gsl_vector_get(x,0) > (1.1*MINPOP)))
      {
        test_ratio = get_groundlevelpop(modelgridindex,element,ion)/gsl_vector_get(x,0);
        if (test_ratio < 1)
        {
          test_ratio = 1./test_ratio;
        }
      }
      else
      {
        test_ratio = 0.0;
      }

      if ((get_groundlevelpop(modelgridindex,element,ion+1) > (1.1*MINPOP)) && (gsl_vector_get(x,nlte_size-1) > (1.1*MINPOP)))
      {
        test_ratio_upper = get_groundlevelpop(modelgridindex,element,ion+1) * modelgrid[modelgridindex].composition[element].partfunct[ion+1]
          / stat_weight(element,ion+1,0) / gsl_vector_get(x,nlte_size-1);
        if (test_ratio_upper < 1)
        {
          test_ratio_upper = 1./test_ratio_upper;
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


//*****************************************************************
double get_tot_nion(int modelgridindex)
{
  double result = 0.;
  for (int element = 0; element < nelements; element++)
  {
    result += modelgrid[modelgridindex].composition[element].abundance / elements[element].mass * get_rho(modelgridindex);

    //int nions = get_nions(element);
    //for (ion = 0; ion < nions; ion++)
    //{
    //  result += ionstagepop(modelgridindex,element,ion);
    //}
  }

  return result;
}

//*****************************************************************
double get_oneoverw(int element, int ion, int modelgridindex)
{
  // Routine to compute the work per ion pair for doing the NT ionization calculation. Makes use of EXTREMELY SIMPLE approximations - high energy limits only */

  // Work in terms of 1/W since this is actually what we want. It is given by sigma/(Latom + Lelec).
  // We are going to start by taking all the high energy limits and ignoring Lelec, so that the
  // denominator is extremely simplified. Need to get the mean Z value.

  double Zbar = 0.0;
  for (int ielement = 0; ielement < nelements; ielement++)
  {
    Zbar += modelgrid[modelgridindex].composition[ielement].abundance * elements[ielement].anumber;
  }
  //printout("cell %d has Zbar of %g\n", modelgridindex, Zbar);

  double Aconst = 1.33e-14 * EV * EV;
  double binding = get_mean_binding_energy(element, ion);
  double oneoverW = Aconst * binding / Zbar / (2*3.14159*pow(QE,4.0));
  //printout("For element %d ion %d I got W of %g (eV)\n", element, ion, 1./oneoverW/EV);

  return oneoverW;
}


//*****************************************************************
double get_mean_binding_energy(int element, int ion)
{
  int q[M_NT_SHELLS];
  double total;

  int ioncharge = get_ionstage(element,ion) - 1;
  int nbound = elements[element].anumber - ioncharge; //number of bound electrons

  if (nbound > 0)
  {
    double use1, use2, use3;
    for (int electron_loop = 0; electron_loop < M_NT_SHELLS; electron_loop++)
    {
      q[electron_loop] = 0;
    }

    for (int electron_loop = 0; electron_loop < nbound; electron_loop++)
    {
      if (q[0] < 2) //K 1s
      {
        q[0] += 1;
      }
      else if (q[1] < 2) //L1 2s
      {
        q[1] += 1;
      }
      else if (q[2] < 2) //L2 2p[1/2]
      {
        q[2] += 1;
      }
      else if (q[3] < 4) //L3 2p[3/2]
      {
        q[3] += 1;
      }
      else if (q[4] < 2) //M1 3s
      {
        q[4] += 1;
      }
      else if (q[5] < 2) //M2 3p[1/2]
      {
        q[5] += 1;
      }
      else if (q[6] < 4) //M3 3p[3/2]
      {
        q[6] += 1;
      }
      else if (ioncharge == 0)
      {
        if (q[9] < 2) //N1 4s
        {
          q[9] += 1;
        }
        else if (q[7] < 4) //M4 3d[3/2]
        {
          q[7] += 1;
        }
        else if (q[8] < 6) //M5 3d[5/2]
        {
          q[8] += 1;
        }
        else
        {
          printout("Going beyond the 4s shell in NT calculation. Abort!\n");
          exit(0);
        }
      }
      else if (ioncharge == 1)
      {
        if (q[9] < 1) // N1 4s
        {
          q[9] += 1;
        }
        else if (q[7] < 4) //M4 3d[3/2]
        {
          q[7] += 1;
        }
        else if (q[8] < 6) //M5 3d[5/2]
        {
          q[8] += 1;
        }
        else
        {
          printout("Going beyond the 4s shell in NT calculation. Abort!\n");
          exit(0);
        }
      }
      else if (ioncharge > 1)
      {
        if (q[7] < 4) //M4 3d[3/2]
        {
          q[7] += 1;
        }
        else if (q[8] < 6) //M5 3d[5/2]
        {
          q[8] += 1;
        }
        else
        {
          printout("Going beyond the 4s shell in NT calculation. Abort!\n");
          exit(0);
        }
      }
    }

    //      printout("For element %d ion %d I got q's of: %d %d %d %d %d %d %d %d %d %d\n", element, ion, q[0], q[1], q[2], q[3], q[4], q[5], q[6], q[7], q[8], q[9]);
    //printout("%g %g %g %g %g %g %g %g %g %g\n", electron_binding[elements[element].anumber-1][0], electron_binding[elements[element].anumber-1][1], electron_binding[elements[element].anumber-1][2],electron_binding[elements[element].anumber-1][3],electron_binding[elements[element].anumber-1][4],electron_binding[elements[element].anumber-1][5],electron_binding[elements[element].anumber-1][6],electron_binding[elements[element].anumber-1][7],electron_binding[elements[element].anumber-1][8],electron_binding[elements[element].anumber-1][9]);

    total = 0.0;
    for (int electron_loop = 0; electron_loop < M_NT_SHELLS; electron_loop++)
    {
      if ((use1 = q[electron_loop]) > 0)
      {
        if ((use2 = electron_binding[elements[element].anumber-1][electron_loop]) > 0)
        {
          if (use2 < (use3 = elements[element].ions[ion].ionpot))
          {
            total += use1/use3;
          }
          else
          {
            total += use1/use2;
          }
        }
        else
        {
          use2 = electron_binding[elements[element].anumber-1][electron_loop-1];
          if (use2 < (use3 = elements[element].ions[ion].ionpot))
          {
            total += use1/use3;
          }
          else
          {
            total += use1/use2;
          }
          //		  total += use1/electron_binding[elements[element].anumber-1][electron_loop-1];
          if (electron_loop != 8)
          {
            //For some reason in the Lotz data, this is no energy for the M5 shell before Ni. So if the complaint
            //is for 8 (corresponding to that shell) then just use the M4 value
            printout("Huh? I'm trying to use a binding energy when I have no data. element %d ion %d\n",element,ion);
            exit(0);
          }
        }
      }
      //printout("total %g\n", total);
    }

  }
  else
  {
    total = 0.0;
  }

  //printout("For element %d ion %d I got mean binding energy of %g (eV)\n", element, ion, 1./total/EV);

  return total;
}


//*****************************************************************
int read_binding_energies()
{
  FILE *binding;
  if ((binding = fopen("binding_energies.txt", "r")) == NULL)
  {
    printout("Cannot open binding_energies.txt.\n");
    exit(0);
  }

  int dum1, dum2;
  fscanf(binding, "%d %d", &dum1, &dum2); //dimensions of the table
  if ((dum1 != M_NT_SHELLS) || (dum2 != MAX_Z_BINDING))
  {
    printout("Wrong size for the binding energy tables!\n");
    exit(0);
  }

  for (int index1 = 0; index1 < dum2; index1++)
  {
    float dum[10];
    fscanf(binding, "%g %g %g %g %g %g %g %g %g %g", &dum[0],&dum[1], &dum[2],&dum[3],&dum[4],&dum[5],&dum[6],&dum[7],&dum[8],&dum[9]);
    for (int index2 = 0; index2 < 10; index2++)
    {
      electron_binding[index1][index2] = dum[index2]*EV;
    }
  }

  fclose(binding);
  return 0;
}


//*****************************************************************
double nt_ionization_rate(int modelgridindex, int element, int ion)
{
  double gammadeposition = rpkt_emiss[modelgridindex] * 1.e20 * 4. * PI;
  // Above is the gamma-ray bit. Below is *supposed* to be the kinetic energy of positrons created by 56Co and 48V. These formulae should be checked, however.
  double positroncobalt = (0.610*0.19*MEV) *
          (exp(-1.*time_step[nts_global].mid/TCOBALT) - exp(-1.*time_step[nts_global].mid/TNICKEL)) /
          (TCOBALT-TNICKEL) * modelgrid[modelgridindex].fni * get_rho(modelgridindex) / MNI56;
  double positron48v = (0.290*0.499*MEV) *
          (exp(-1.*time_step[nts_global].mid/T48V) - exp(-1.*time_step[nts_global].mid/T48CR)) /
          (T48V-T48CR) * modelgrid[modelgridindex].f48cr * get_rho(modelgridindex) / MCR48;

  //printout("nt_ionization_rate: element: %d, ion %d\n",element,ion);
  //printout("nt_ionization_rate: gammadep: %g, poscobalt %g pos48v %g\n",
  //  gammadeposition,positroncobalt,positron48v);

  // to get the non-thermal ionization rate we need to divide the energy deposited
  // per unit volume per unit time in the grid cell (sum of terms above)
  // by the total ion number density and the "work per ion pair"
  return (gammadeposition + positroncobalt + positron48v) /
    get_tot_nion(modelgridindex) * get_oneoverw(element, ion, modelgridindex);
}


//***************************************************************/
double superlevel_boltzmann(int modelgridindex, int element, int ion, int level)
{
  double T_exc = get_TJ(modelgridindex);
  double E_level = epsilon(element,ion,level);
  double E_superlevel = epsilon(element,ion,get_nlevels_nlte(element,ion)+1);

  return (stat_weight(element,ion,level) * exp(-(E_level-E_superlevel)/KB/T_exc));
}




///***************************************************************************/
double nlte_matrix_rad_deexcitation(int modelgridindex, int lower, double epsilon_trans, int lineindex, double t_current)
///radiative deexcitation rate: paperII 3.5.2
{
  int element = mastate[tid].element;
  int ion = mastate[tid].ion;
  int upper = mastate[tid].level;

  //int modelgridindex = cell[pkt_ptr->where].modelgridindex;

  #ifdef DEBUG_ON
    if (upper <= lower)
    {
      printout("[fatal] rad_deexcitation: tried to calculate upward transition ... abort\n");
      abort();
    }
  #endif

  //double statweight_target = statw_down(lineindex);
  double nu_trans = epsilon_trans/H;
  //A_ul = einstein_spontaneous_emission(element,ion,upper,lower);
  double A_ul = einstein_spontaneous_emission(lineindex);
  double B_ul = CLIGHTSQUAREDOVERTWOH / pow(nu_trans,3) * A_ul;
  //double B_lu = mastate[tid].statweight/statweight_target * B_ul;

  //double n_u = get_levelpop(modelgridindex,element,ion,upper);
  //double n_l = get_levelpop(modelgridindex,element,ion,lower);
  double n_u = calculate_exclevelpop(modelgridindex,element,ion,upper);
  double n_l = calculate_exclevelpop(modelgridindex,element,ion,lower); //TODO: do these two ways of calculating level pops give different results?
  double T_R = get_Te(modelgridindex);
  double W   = get_W(modelgridindex);
  //double tau_sobolev = (B_lu * n_l - B_ul * n_u) * HCLIGHTOVERFOURPI * t_current;

  //double beta = 1.0;
  double R = 0.0;

  if ((n_l < 1.1 * MINPOP) || (n_u < 1.1 * MINPOP))
  {
    //beta = 1.0;
    R = 0.0;
  }
  else
  {
    //beta = 1.0 / tau_sobolev * (1 - exp(-tau_sobolev));
    R = A_ul  + B_ul * W * TWOHOVERCLIGHTSQUARED * pow(nu_trans,3) * 1.0/(exp(HOVERKB*nu_trans/T_R) - 1);
  }


  #ifdef DEBUG_ON
    if (debuglevel == 2) printout("[debug] rad_rates_down: element, ion, upper, lower %d, %d, %d, %d\n",element,ion,upper,lower);
    //printout("[debug] rad_rates_down: tau_sobolev, beta %g, %g\n",tau_sobolev,beta);
    //if (debuglevel == 2) printout("[debug] rad_deexc: A_ul %g, tau_sobolev %g, n_u %g\n",A_ul,tau_sobolev,n_u);
    //if (debuglevel == 777) printout("[debug] rad_deexc: A_ul %g, tau_sobolev %g, n_u %g\n",A_ul,tau_sobolev,n_u);
    if (!isfinite(R))
    {
      printout("fatal a1: abort\n");
      abort();
    }
  #endif

  return R;
}



double nlte_matrix_rad_excitation(int modelgridindex, int upper, double epsilon_trans, int lineindex, double t_current)//, double T_R, double W)
///radiative excitation rate: paperII 3.5.2
{
  int element = mastate[tid].element;
  int ion = mastate[tid].ion;
  int lower = mastate[tid].level;

  //int modelgridindex = cell[pkt_ptr->where].modelgridindex;

  #ifdef DEBUG_ON
    if (upper <= lower)
    {
      printout("[fatal] rad_excitation: tried to calculate downward transition ... abort\n");
      abort();
    }
  #endif

  double statweight_target = statw_up(lineindex);

  double nu_trans = epsilon_trans/H;
  //A_ul = einstein_spontaneous_emission(element,ion,upper,lower);
  double A_ul = einstein_spontaneous_emission(lineindex);
  double B_ul = CLIGHTSQUAREDOVERTWOH / pow(nu_trans,3) * A_ul;
  double B_lu = statweight_target / mastate[tid].statweight * B_ul;

  //double n_u = get_levelpop(modelgridindex,element,ion,upper);
  //double n_l = get_levelpop(modelgridindex,element,ion,lower);
  double n_u = calculate_exclevelpop(modelgridindex,element,ion,upper);
  double n_l = calculate_exclevelpop(modelgridindex,element,ion,lower); //TODO: do these two ways of calculating level pops give different results?
  double T_R = get_Te(modelgridindex);
  double W   = get_W(modelgridindex);
  //n_u = n_l * W * g_ratio * exp(-epsilon_trans/KB/T_R);
  //double tau_sobolev = (B_lu * n_l - B_ul * n_u) * HCLIGHTOVERFOURPI * t_current;

  //double beta = 1.0;
  double R = 0.0;
  if ((n_l < 1.1 * MINPOP) || (n_u < 1.1 * MINPOP))
  {
    R = 0.0;
  }
  else
  {
    //beta = 1.0 / tau_sobolev * (1. - exp(-tau_sobolev));
    //printout("[check] rad_excitation: %g, %g, %g\n",1.0/tau_sobolev,exp(-tau_sobolev),1.0/tau_sobolev * (1. - exp(-tau_sobolev)));
    //n_u2 = calculate_levelpop_fromreflevel(pkt_ptr->where,element,ion,upper,lower,mastate[tid].nnlevel);
    //R = (B_lu*mastate[tid].nnlevel - B_ul * n_u2) * beta * radfield(nu_trans,pkt_ptr->where);
    R = B_lu * W * TWOHOVERCLIGHTSQUARED * pow(nu_trans,3) * 1.0/(exp(HOVERKB*nu_trans/T_R) - 1);
  }


  #ifdef DEBUG_ON
    //printout("tau_sobolev, beta: %g, %g\n",tau_sobolev,beta);
    if (debuglevel == 2) printout("[debug] rad_rates_up: element, ion, upper, lower, A_ul, n_u: %d, %d, %d, %d, %g, %g\n",element,ion,upper,lower,A_ul,n_l);
    //if (debuglevel == 2) printout("[debug] rad_exc: A_ul %g, tau_sobolev %g, n_u %g, n_l %g, radfield %g\n",A_ul,tau_sobolev,n_u,n_l,radfield(nu_trans,modelgridindex));
    //if (debuglevel == 777) printout("[debug] rad_exc: A_ul %g, tau_sobolev %g, n_u %g, n_l %g, radfield %g\n",A_ul,tau_sobolev,n_u,n_l,radfield(nu_trans,modelgridindex));
    if (!isfinite(R))
    {
      printout("[fatal] rad_excitation: abort\n");
      //printout("[fatal] rad_excitation: R %g, mastate[tid].nnlevel %g, B_lu %g, B_ul %g, n_u %g, n_l %g, beta %g, radfield %g,tau_sobolev %g, t_current %g\n",R,mastate[tid].nnlevel,B_lu,B_ul,n_u,n_l,beta,radfield(nu_trans,modelgridindex),tau_sobolev,t_current);
      //printout("[fatal] rad_excitation: %g, %g, %g\n",1.0/tau_sobolev,exp(-tau_sobolev),1.0/tau_sobolev * (1. - exp(-tau_sobolev)));
      abort();
    }
  #endif

  return R;
}
