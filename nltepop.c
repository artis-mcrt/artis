#include "sn3d.h"
#include "nltepop.h"

double nlte_pops(int element, int ion, int modelgridindex, int timestep)
//solves for nlte correction factors to level populations for levels
{
  int nlte_levels;
  int level, upper, lower, level_use, lower_use, upper_use;

  double* rate_matrix;
  double* balance_vector;

  int ionisinglevels, ndowntrans, nuptrans, i;
  double statweight, epsilon_current;
  double tmid;
  double superlevel_partition;
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
  //double get_nne(int modelgridindex);
  //double get_Te(int modelgridindex);
  double T_e, nne;

  if (get_nlevels(element,ion) > 1)
  {
    nne = get_nne(modelgridindex);
    T_e = get_Te(modelgridindex);

    printout("Solving for NLTE populations in cell %d. Doing element %d, ion %d. I think it's timestep %d\n", modelgridindex, element, ion, timestep);
    //printout("Current Te %g and ne %g\n",T_e, nne);

    nlte_levels = get_nlevels_nlte(element,ion);
    t_mid = time_step[timestep].mid;

    dummy.where = modelgridindex;

    if (nlte_levels == (get_nlevels(element, ion)-1))
    {
      nlte_size = nlte_levels+2;
      super_level = 0;
    }
    else
    {
      nlte_size = nlte_levels+3;
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

      //printf("rate %p balane %p NULL %p\n", rate_matrix, balance_vector, NULL);

      //printout("I think there are %d levels to deal with and managed to allocate memory.\n", nlte_size);

      superlevel_partition = 0.0;
      for (level = 1; level < get_nlevels(element,ion); level++)
      {
        if (is_nlte(element, ion, level) == 0)
        {
          superlevel_partition += superlevel_boltzmann(modelgridindex,element,ion,level);
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

      for (level = 0; level < nlte_size; level++)
      {
        balance_vector[level] = 0.0;
      }

      for (level = 0; level < get_nlevels(element,ion); level++)
      {
        statweight = stat_weight(element,ion,level);
        epsilon_current = epsilon(element,ion,level);
        ndowntrans = elements[element].ions[ion].levels[level].downtrans[0].targetlevel;
        nuptrans = elements[element].ions[ion].levels[level].uptrans[0].targetlevel;

        /* deexcitation */
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

          R = rad_deexcitation(&dummy,lower,epsilon_trans,statweight_target,lineindex,t_mid);
          C = col_deexcitation(modelgridindex,lower,epsilon_trans,statweight_target,lineindex);

          s_renorm = 1.0;

          if ((level == 0) || (is_nlte(element, ion, level) == 1))
          {
            level_use = level;
          }
          else
          {
            level_use = nlte_levels+1;
            s_renorm = superlevel_boltzmann(modelgridindex,element, ion, level)/superlevel_partition;
          }
          if ((lower == 0) || (is_nlte(element, ion, lower) == 1))
          {
            lower_use = lower;
          }
          else
          {
            lower_use = nlte_levels+1;
          }

          rate_matrix[level_use*(nlte_size) + level_use] -= (R + C) * s_renorm;
          rate_matrix[lower_use*(nlte_size) + level_use] += (R + C) * s_renorm;
          //printout("First using %d of %d\n", level_use*(nlte_size) + level_use, nlte_size*nlte_size);
          //printout("Second using %d of %d\n", lower_use*(nlte_size) + level_use, nlte_size*nlte_size);
        }

        /* excitation */
        for (i = 1; i <= nuptrans; i++)
        {
          upper = elements[element].ions[ion].levels[level].uptrans[i].targetlevel;
          epsilon_target = elements[element].ions[ion].levels[level].uptrans[i].epsilon;
          statweight_target = elements[element].ions[ion].levels[level].uptrans[i].stat_weight;
          lineindex = elements[element].ions[ion].levels[level].uptrans[i].lineindex;
          epsilon_trans = epsilon_target - epsilon_current;

          mastate[tid].element = element;
          mastate[tid].ion = ion;
          mastate[tid].level = level;
          mastate[tid].statweight = statweight;
          mastate[tid].nnlevel = 1.0;

          R = rad_excitation(&dummy,upper,epsilon_trans,statweight_target,lineindex,t_mid);//,T_R,W);
          C = col_excitation(modelgridindex,upper,lineindex,epsilon_trans);

          s_renorm = 1.0;

          if ((level == 0) || (is_nlte(element, ion, level) == 1))
          {
            level_use = level;
          }
          else
          {
            level_use = nlte_levels+1;
            s_renorm = superlevel_boltzmann(modelgridindex,element, ion, level)/superlevel_partition;
          }
          if ((upper == 0) || (is_nlte(element, ion, upper) == 1))
          {
            upper_use = upper;
          }
          else
          {
            upper_use = nlte_levels+1;
          }
          rate_matrix[level_use*(nlte_size) + level_use] -= (R + C) * s_renorm;
          rate_matrix[upper_use*(nlte_size) + level_use] += (R + C) * s_renorm;
          //printout("Third using %d of %d\n", level_use*(nlte_size) + level_use, nlte_size*nlte_size);
          //printout("Fourth using %d of %d\n", upper_use*(nlte_size) + level_use, nlte_size*nlte_size);

        }

        #ifdef NT_ON
          if (ion < get_nions(element)-1)
          {
            Y = nt_ionization_rate(modelgridindex, element, ion);
            s_renorm = 1.0;

            if ((level == 0) || (is_nlte(element, ion, level) == 1))
            {
              level_use = level;
            }
            else
            {
              level_use = nlte_levels+1; //the super level
              s_renorm = superlevel_boltzmann(modelgridindex,element, ion, level)/superlevel_partition;
            }

            upper_use = nlte_size-1; //the continuum
            rate_matrix[level_use*(nlte_size) + level_use] -= (Y)*s_renorm;
            rate_matrix[upper_use*(nlte_size) + level_use] += (Y)*s_renorm;
          }
        #endif

        //now put in the photoionization/recombination processes
        ionisinglevels = get_bfcontinua(element,ion);
        if (ion < get_nions(element)-1 && level < ionisinglevels)  //&& get_ionstage(element,ion) < get_element(element)+1)
        {
          mastate[tid].element = element;
          mastate[tid].ion = ion;
          mastate[tid].level = level;
          mastate[tid].statweight = statweight;
          mastate[tid].nnlevel = 1.0;

          R = 0.0;
          C = 0.0;
          for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion,level); phixstargetindex++)
          {
            upper = get_phixsupperlevel(element,ion,level,phixstargetindex);
            epsilon_trans = epsilon(element,ion+1,upper) - epsilon_current;
            R += photoionization(modelgridindex,phixstargetindex,epsilon_trans);
            C += col_ionization(modelgridindex,phixstargetindex,epsilon_trans);
          }

          s_renorm = 1.0;

          if ((level == 0) || (is_nlte(element, ion, level) == 1))
          {
            level_use = level;
          }
          else
          {
            level_use = nlte_levels+1; //the super level
            s_renorm = superlevel_boltzmann(modelgridindex,element,ion,level) / superlevel_partition;
          }

          upper_use = nlte_size-1; //the continuum

          rate_matrix[level_use*(nlte_size) + level_use] -= (R + C) * s_renorm;
          rate_matrix[upper_use*(nlte_size) + level_use] += (R + C) * s_renorm;

          s_renorm = 1.0;

          R = 0.0;
          C = 0.0;
          mastate[tid].element = element;
          mastate[tid].ion = ion+1;
          mastate[tid].nnlevel = 1.0;
          for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion,level); phixstargetindex++)
          {
            upper = get_phixsupperlevel(element,ion,level,phixstargetindex);
            mastate[tid].level = upper;
            mastate[tid].statweight = stat_weight(element,ion+1,upper);
            epsilon_trans = epsilon(element,ion+1,upper) - epsilon_current;
            R += rad_recombination(modelgridindex,level,epsilon_trans);
            //printout("rad recombination of element %d, ion %d, level %d, to lower level %d has rate %g (ne %g and Te %g)\n",element,ion,mastate[tid].level,level,R/nne,nne,T_e);
            printout("%d %d %d %d %g %g %g \n",element,ion,mastate[tid].level,level,R/nne,nne,T_e);
            C += col_recombination(modelgridindex,level,epsilon_trans);
            //C=C*1.e-10;
          }

          rate_matrix[upper_use*(nlte_size) + upper_use] -= (R + C) * s_renorm;
          rate_matrix[level_use*(nlte_size) + upper_use] += (R + C) * s_renorm;

          //balance_vector[level_use] += -1. * get_groundlevelpop(modelgridindex,element,ion+1) * (R + C);
        }

      }


      //printout("I've filled out a rate matrix. Probably you should check it at some point!\n");

      /*
      printout("Rate matrix:\n");
      for (level = 0; level < nlte_size; level++)
        {
          for (level_use = 0; level_use < nlte_size; level_use++)
      {
        printout("%g ", rate_matrix[level*nlte_size + level_use]);
      }
          printout("\n");
        }
      printout("\n");
      */

      /*replace the first line of the matrix with the normalisation constraint*/
      for (level = 1; level < nlte_size; level++)
      {
        rate_matrix[level] = 1.0;
        balance_vector[level] = 0.0;
      }
      balance_vector[0] = get_groundlevelpop(modelgridindex,element,ion)*modelgrid[modelgridindex].composition[element].partfunct[ion] / stat_weight(element,ion,0);
      rate_matrix[0] = 1.0;
      rate_matrix[nlte_size-1] = 0.0;

      for (level = 0; level < nlte_size; level++)
      {
        for (level_use = 0; level_use < nlte_size; level_use++)
        {
          //printout("%g ", rate_matrix[level*nlte_size + level_use]);
          if (!isfinite(rate_matrix[level*nlte_size + level_use]))
          {
            printout("[fatal]: NLTE matrix with non-finite element: %d %d %g\n", level, level_use, rate_matrix[level*nlte_size + level_use]);
            printout("[fatal]: found when handling element %d and ion %d\n", element, ion);
            printout("[fatal]: the relevant ground state populations are %g and %g\n",get_groundlevelpop(modelgridindex,element,ion), get_groundlevelpop(modelgridindex,element,ion+1));
            exit(0);
          }
        }
        // printout("\n");
      }
      //printout("\n");

      for (level = 0; level < nlte_size; level++)
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

      //printf("rate %p balane %p NULL %p\n", rate_matrix, balance_vector, NULL);


      /*next bunch of lines are based on a gsl example */

      gsl_matrix_view m
        = gsl_matrix_view_array (rate_matrix, nlte_size, nlte_size);

      gsl_vector_view b
        = gsl_vector_view_array (balance_vector, nlte_size);

      gsl_vector *x = gsl_vector_alloc (nlte_size);

      int s;

      gsl_permutation * p = gsl_permutation_alloc (nlte_size);

      gsl_linalg_LU_decomp (&m.matrix, p, &s);

      gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
      //gsl_linalg_HH_solve (&m.matrix, &b.vector, x);

      //printf("rate %p balane %p NULL %p\n", rate_matrix, balance_vector, NULL);

      //printout("I did the GSL stuff - didnt' blow up at least!\n");

      //printout("The ground state populations were %g, %g, %g and %g\n", get_groundlevelpop(modelgridindex,element,ion), gsl_vector_get(x,0),  get_groundlevelpop(modelgridindex,element,ion+1), gsl_vector_get(x,nlte_size-1));

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
        test_ratio=test_ratio_upper;
      }

      printout("The test ratios are %g %g. Passing %g.\n", get_groundlevelpop(modelgridindex,element,ion)/gsl_vector_get(x,0), get_groundlevelpop(modelgridindex,element,ion+1)*modelgrid[modelgridindex].composition[element].partfunct[ion+1]/stat_weight(element,ion+1,0)/gsl_vector_get(x,nlte_size-1), test_ratio);
      //	  printout("The test ratio is %g. Passing %g.\n", get_groundlevelpop(modelgridindex,element,ion+1)*modelgrid[modelgridindex].composition[element].partfunct[ion+1]/stat_weight(element,ion+1,0)/gsl_vector_get(x,nlte_size-1), test_ratio);


      //printout("The top five excited states were %g, %g, %g, %g and %g.\n",gsl_vector_get(x,nlte_levels-4),gsl_vector_get(x,nlte_levels-3),gsl_vector_get(x,nlte_levels-2),gsl_vector_get(x,nlte_levels-1),gsl_vector_get(x,nlte_levels));

      //printout("The first five excited states were %g, %g, %g, %g and %g.\n",gsl_vector_get(x,1),gsl_vector_get(x,2),gsl_vector_get(x,3),gsl_vector_get(x,4),gsl_vector_get(x,5));

      //lag = 0.0;

      //if ((check=(fabs((1.0 - test_ratio)/(1.0 + test_ratio))) > 0.05) && (check < 0.2))
      //{
      //   lag=1.0;
      //}


      /* Write the NLTE level populations to the array*/
      for (level=1; level < nlte_levels+1; level++)
      {
        modelgrid[modelgridindex].nlte_pops[nlte_start+level-1] = gsl_vector_get(x,level)/modelgrid[modelgridindex].rho;
        //printout("I have interfered with index %d.\n", nlte_start+level-1);
        //modelgrid[modelgridindex].nlte_pops[nlte_start+level-1] = ((lag*modelgrid[modelgridindex].nlte_pops[nlte_start+level-1]) + gsl_vector_get(x,level))/(lag + 1.0)/modelgrid[modelgridindex].rho;
      }
      /* If there is a superlevel then write that too*/

      if (super_level == 1)
      {
        //printout("I thought the super level was: %g\n", modelgrid[modelgridindex].nlte_pops[nlte_start+nlte_levels]);

        modelgrid[modelgridindex].nlte_pops[nlte_start+nlte_levels] = gsl_vector_get(x,nlte_levels+1)/modelgrid[modelgridindex].rho/superlevel_partition;
        //	      modelgrid[modelgridindex].nlte_pops[nlte_start+nlte_levels] = ((lag*modelgrid[modelgridindex].nlte_pops[nlte_start+nlte_levels]) + gsl_vector_get(x,nlte_levels+1))/(lag + 1.0)/modelgrid[modelgridindex].rho/superlevel_partition;

        //printout("Now I think it is: %g\n", modelgrid[modelgridindex].nlte_pops[nlte_start+nlte_levels]);
        //printout("I also interfered with index %d.\n", nlte_start+nlte_levels);
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

      gsl_permutation_free (p);
      //printout("I freed up p\n");
      gsl_vector_free (x);
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
      for (level = 1; level < nlte_levels+1; level++)
      {
        modelgrid[modelgridindex].nlte_pops[nlte_start+level-1] = -1.0;///flag to indicate no useful data
      }
      if (super_level == 1)
      {
        modelgrid[modelgridindex].nlte_pops[nlte_start+nlte_levels] = -1.0;
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


/*****************************************************************/
double get_tot_nion(int modelgridindex)
{
  //double ionstagepop(int modelgridindex, int element, int ion);
  //int nions;
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


/*****************************************************************/
double get_oneoverw(int element, int ion, int modelgridindex)
{
  /* Routine to compute the work per ion pair for doing the NT ionization calculation. Makes use of EXTREMELY SIMPLE approximations - high energy limits only */

  double Zbar;

  /* Work in terms of 1/W since this is actually what we want. It is given by sigma/(Latom + Lelec).
     We are going to start by taking all the high energy limits and ignoring Lelec, so that the
     denominator is extremely simplified. Need to get the mean Z value. */

  Zbar = 0.0;
  for (int ielement = 0; ielement < nelements; ielement++)
  {
    Zbar += modelgrid[modelgridindex].composition[ielement].abundance * elements[ielement].anumber;
  }
  //printout("cell %d has Zbar of %g\n", modelgridindex, Zbar);

  double binding = get_mean_binding_energy(element, ion);
  double Aconst = 1.33e-14 * EV * EV;
  double oneoverW = Aconst * binding / Zbar / (2*3.14159*pow(QE,4.0));
  //printout("For element %d ion %d I got W of %g (eV)\n", element, ion, 1./oneoverW/EV);

  return oneoverW;
}


/*****************************************************************/
double get_mean_binding_energy(int element, int ion)
{
  int q[M_NT_SHELLS];
  double total, use1, use2, use3;

  int ioncharge = get_ionstage(element,ion) - 1;
  int nbound = elements[element].anumber - ioncharge; //number of bound electrons

  if (nbound > 0)
  {
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
            printout("Huh? I'm trying to use a binding energy when I have no data.\n");
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


/*****************************************************************/
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


/*****************************************************************/
double nt_ionization_rate(int modelgridindex, int element, int ion)
{
  double Y_nt = rpkt_emiss[modelgridindex] * 1.e20 * 4. * PI;
  // Above is the gamma-ray bit. Below is *supposed* to be the kinetic energy of positrons created by 56Co and 48V. These formulae should be checked, however.
  Y_nt += (0.610*0.19*MEV)*(exp(-1.*time_step[nts_global].mid/TCOBALT) - exp(-1.*time_step[nts_global].mid/TNICKEL))/(TCOBALT-TNICKEL)*modelgrid[modelgridindex].fni*get_rho(modelgridindex)/MNI56;
  Y_nt += (0.290*0.499*MEV)*(exp(-1.*time_step[nts_global].mid/T48V) - exp(-1.*time_step[nts_global].mid/T48CR))/(T48V-T48CR)*modelgrid[modelgridindex].f48cr*get_rho(modelgridindex)/MCR48;

  //this is the energy deposited per unit volume per unit time in the grid cell
  //to get the non-thermal ionization rate we need to divide this by the total ion number density and the "work per ion pair"
  Y_nt = Y_nt / get_tot_nion(modelgridindex) * get_oneoverw(element, ion, modelgridindex);

  return Y_nt;
}


/***************************************************************/
double superlevel_boltzmann(int modelgridindex, int element, int ion, int level)
{
  double T_exc = get_TJ(modelgridindex);
  double E_level = epsilon(element,ion,level);
  double E_ground = epsilon(element,ion,get_nlevels_nlte(element, ion));

  return (stat_weight(element,ion,level) * exp(-(E_level-E_ground)/KB/T_exc));
}
