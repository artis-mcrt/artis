#include <stdio.h>
//#include <math.h>
//#include <stdlib.h>
#include "assert.h"
#include "sn3d.h"
#include "atomic.h"
#include "gamma.h"
#include "grid_init.h"
#include "input.h"
#include "linelist.h"
#include "nltepop.h"
#include "rpkt.h"
#ifdef DO_EXSPEC
  #include "exspec.h"
#endif

const int nlevels_requiretransitions = 0; // first n levels will be collisionally coupled to all other levels

const int groundstate_index_in = 1; // starting level index in the input files

typedef struct
{
  int lower;
  int upper;
  double A;
  double coll_str;
  bool forbidden;
} transitiontable_entry;  /// only used temporarily during input


static void read_phixs_data(void)
{
  printout("readin phixs data\n");

  FILE *restrict phixsdata = fopen("phixsdata_v2.txt", "r");
  if (phixsdata == NULL)
  {
    printout("Cannot open phixsdata_v2.txt.\n");
    abort();
  }

  fscanf(phixsdata,"%d\n",&NPHIXSPOINTS);
  assert(NPHIXSPOINTS > 0);
  fscanf(phixsdata,"%lg\n",&NPHIXSNUINCREMENT);
  assert(NPHIXSNUINCREMENT > 0.);
  int Z;
  int upperionstage;
  int upperlevel_in;
  int lowerionstage;
  int lowerlevel_in;
  double phixs_threshold_ev;
  while (fscanf(phixsdata,"%d %d %d %d %d %lg\n",&Z,&upperionstage,&upperlevel_in,&lowerionstage,&lowerlevel_in,&phixs_threshold_ev) != EOF)
  {
    assert(Z > 0);
    assert(upperionstage >= 2);
    assert(lowerionstage >= 1);
    bool skip_this_phixs_table = false;
    //printout("[debug] Z %d, upperion %d, upperlevel %d, lowerion %d, lowerlevel, %d\n",Z,upperion,upperlevel,lowerion,lowerlevel);
    /// translate readin anumber to element index
    int element = get_elementindex(Z);

    /// store only photoionization crosssections for elements that are part of the current model atom
    if (element >= 0)
    {
      /// translate readin ionstages to ion indices

      const int upperion = upperionstage - get_ionstage(element, 0);
      const int lowerion = lowerionstage - get_ionstage(element, 0);
      const int lowerlevel = lowerlevel_in - groundstate_index_in;
      assert(upperion >= 0);
      assert(lowerion >= 0);
      assert(lowerlevel >= 0);
      /// store only photoionization crosssections for ions that are part of the current model atom
      if (lowerion >= 0 && lowerlevel < get_nlevels(element, lowerion) && upperion < get_nions(element))
      {
          //double phixs_threshold_ev = (epsilon(element, upperion, upperlevel) - epsilon(element, lowerion, lowerlevel)) / EV;
          elements[element].ions[lowerion].levels[lowerlevel].phixs_threshold = phixs_threshold_ev * EV;
          if (upperlevel_in >= 0) // file gives photoionisation to a single target state only
          {
            int upperlevel = upperlevel_in - groundstate_index_in;
            assert(upperlevel >= 0);
            elements[element].ions[lowerion].levels[lowerlevel].nphixstargets = 1;
            if ((elements[element].ions[lowerion].levels[lowerlevel].phixstargets = calloc(1, sizeof(phixstarget_entry))) == NULL)
            {
              printout("[fatal] input: not enough memory to initialize phixstargets... abort\n");
              abort();
            }
            if (upperion == get_nions(element) - 1) // top ion has only one level, so send it to that level
              upperlevel = 0;
            elements[element].ions[lowerion].levels[lowerlevel].phixstargets[0].levelindex = upperlevel;
            elements[element].ions[lowerion].levels[lowerlevel].phixstargets[0].probability = 1.0;
          }
          else // upperlevel < 0, indicating that a table of upper levels and their probabilities will follow
          {
            int in_nphixstargets;
            fscanf(phixsdata,"%d\n", &in_nphixstargets);
            assert(in_nphixstargets >= 0);
            // read in a table of target states and probabilities and store them
            if (upperion < get_nions(element) - 1) // not the top ion, which has nlevelsmax = 1
            {
              if ((elements[element].ions[lowerion].levels[lowerlevel].phixstargets = (phixstarget_entry *) calloc(in_nphixstargets, sizeof(phixstarget_entry))) == NULL)
              {
                printout("[fatal] input: not enough memory to initialize phixstargets list... abort\n");
                abort();
              }
              double probability_sum = 0.;
              for (int i = 0; i < in_nphixstargets; i++)
              {
                double phixstargetprobability;
                fscanf(phixsdata, "%d %lg\n", &upperlevel_in, &phixstargetprobability);
                const int upperlevel = upperlevel_in - groundstate_index_in;
                assert(upperlevel >= 0);
                assert(phixstargetprobability > 0);
                elements[element].ions[lowerion].levels[lowerlevel].phixstargets[i].levelindex = upperlevel;
                elements[element].ions[lowerion].levels[lowerlevel].phixstargets[i].probability = phixstargetprobability;
                probability_sum += phixstargetprobability;
              }
              if (fabs(probability_sum - 1.0) > 0.01)
              {
                printout("WARNING: photoionisation table for Z=%d ionstage %d has probabilities that sum to %g", Z, get_ionstage(element, lowerion), probability_sum);
              }
              elements[element].ions[lowerion].levels[lowerlevel].nphixstargets = in_nphixstargets;
            }
            else // file has table of target states and probabilities but our top ion is limited to one level
            {
              for (int i = 0; i < in_nphixstargets; i++)
              {
                double phixstargetprobability;
                fscanf(phixsdata, "%d %lg\n", &upperlevel_in, &phixstargetprobability);
              }
              if ((elements[element].ions[lowerion].levels[lowerlevel].phixstargets = calloc(1, sizeof(phixstarget_entry))) == NULL)
              {
                printout("[fatal] input: not enough memory to initialize phixstargets... abort\n");
                abort();
              }
              // send it to the ground state of the top ion
              elements[element].ions[lowerion].levels[lowerlevel].nphixstargets = 1;
              elements[element].ions[lowerion].levels[lowerlevel].phixstargets[0].levelindex = 0;
              elements[element].ions[lowerion].levels[lowerlevel].phixstargets[0].probability = 1.0;
            }
          }

          /// The level contributes to the ionisinglevels if its energy
          /// is below the ionisiation potential and the level doesn't
          /// belong to the topmost ion included.
          /// Rate coefficients are only available for ionising levels.
          //  also need (levelenergy < ionpot && ...)?
          if (lowerion < get_nions(element)-1) ///thats only an option for pure LTE && level < TAKE_N_BFCONTINUA)
          {
            for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,lowerion,lowerlevel); phixstargetindex++)
            {
              if ((elements[element].ions[lowerion].levels[lowerlevel].phixstargets[phixstargetindex].spontrecombcoeff = calloc(TABLESIZE, sizeof(double))) == NULL)
              {
                printout("[fatal] input: not enough memory to initialize spontrecombcoeff table for element %d, ion %d, level %d\n",element,lowerion,lowerlevel);
                abort();
              }
              #if (!NO_LUT_PHOTOION)
              if ((elements[element].ions[lowerion].levels[lowerlevel].phixstargets[phixstargetindex].corrphotoioncoeff = calloc(TABLESIZE, sizeof(double))) == NULL)
              {
                printout("[fatal] input: not enough memory to initialize photoioncoeff table for element %d, ion %d, level %d\n",element,lowerion,lowerlevel);
                abort();
              }
              #endif
              #if (!NO_LUT_BFHEATING)
              if ((elements[element].ions[lowerion].levels[lowerlevel].phixstargets[phixstargetindex].bfheating_coeff = calloc(TABLESIZE, sizeof(double))) == NULL)
              {
                printout("[fatal] input: not enough memory to initialize modified_photoioncoeff table for element %d, ion %d, level %d\n",element,lowerion,lowerlevel);
                abort();
              }
              #endif
              if ((elements[element].ions[lowerion].levels[lowerlevel].phixstargets[phixstargetindex].bfcooling_coeff = calloc(TABLESIZE, sizeof(double))) == NULL)
              {
                printout("[fatal] input: not enough memory to initialize modified_photoioncoeff table for element %d, ion %d, level %d\n",element,lowerion,lowerlevel);
                abort();
              }
            }
          }

          if ((elements[element].ions[lowerion].levels[lowerlevel].photoion_xs = (float *) calloc(NPHIXSPOINTS,sizeof(float))) == NULL)
          {
            printout("[fatal] input: not enough memory to initialize photoion_xslist... abort\n");
            abort();
          }
          for (int i = 0; i < NPHIXSPOINTS; i++)
          {
            float phixs;
            fscanf(phixsdata,"%g\n", &phixs);
            assert(phixs >= 0);

            ///the photoionisation cross-sections in the database are given in Mbarn = 1e6 * 1e-28m^2
            ///to convert to cgs units multiply by 1e-18
            elements[element].ions[lowerion].levels[lowerlevel].photoion_xs[i] = phixs * 1e-18;
            //fprintf(database_file,"%g %g\n", nutable[i], phixstable[i]);
          }

          //nbfcontinua++;
          //printout("[debug] element %d, ion %d, level %d: phixs exists %g\n",element,lowerion,lowerlevel,phixs*1e-18);
          skip_this_phixs_table = false;
      }
      else
      {
        skip_this_phixs_table = true;
      }
    }
    else
    {
      skip_this_phixs_table = true;
    }
    if (skip_this_phixs_table) // for ions or elements that are not part of the current model atom, proceed through the lines and throw away the data
    {
      if (upperlevel_in < 0) // a table of target states and probabilities will follow, so read past those lines
      {
        int nphixstargets;
        fscanf(phixsdata,"%d\n",&nphixstargets);
        for (int i = 0; i < nphixstargets; i++)
        {
          double phixstargetprobability;
          fscanf(phixsdata,"%d %lg\n",&upperlevel_in,&phixstargetprobability);
        }
      }
      for (int i = 0; i < NPHIXSPOINTS; i++) //skip through cross section list
      {
        float phixs;
        fscanf(phixsdata,"%g\n",&phixs);
      }
    }
  }
  fclose(phixsdata);
}


static void read_ion_levels(
  FILE* adata, const int element, const int ion, const int nions, const int nlevels, const int nlevelsmax,
  const double energyoffset, const double ionpot,
  int *cont_index, int *nuparr, int *ndownarr)
{
  for (int level = 0; level < nlevels; level++)
  {
    int levelindex_in;
    double levelenergy;
    double statweight;
    int ntransitions;
    fscanf(adata, "%d %lg %lg %d%*[^\n]\n", &levelindex_in, &levelenergy, &statweight, &ntransitions);
    assert(levelindex_in == level + groundstate_index_in);
    assert((ion < nions - 1) || (ntransitions > 0) || (nlevels == 1));
    //if (element == 1 && ion == 0) printf("%d %16.10f %g %d\n",levelindex,levelenergy,statweight,ntransitions);
    if (level < nlevelsmax)
    {
      ndownarr[level] = 1;
      nuparr[level] = 1;
      //elements[element].ions[ion].levels[level].epsilon = (energyoffset + levelenergy) * EV;
      const double currentlevelenergy = (energyoffset + levelenergy) * EV;
      //if (element == 1 && ion == 0) printf("%d %16.10e\n",levelindex,currentlevelenergy);
      //printout("energy for level %d of ionstage %d of element %d is %g\n",level,ionstage,element,currentlevelenergy/EV);
      elements[element].ions[ion].levels[level].epsilon = currentlevelenergy;
      //printout("epsilon(%d,%d,%d)=%g",element,ion,level,elements[element].ions[ion].levels[level].epsilon);

      //if (level == 0 && ion == 0) energyoffset = levelenergy;
      elements[element].ions[ion].levels[level].stat_weight = statweight;
      ///Moved to the section with ionising levels below
      //elements[element].ions[ion].levels[level].cont_index = cont_index;
      //cont_index--;
      /// Initialise the metastable flag to true. Set it to false if a downward transition exists.
      elements[element].ions[ion].levels[level].metastable = true;
      //elements[element].ions[ion].levels[level].main_qn = mainqn;

      /// The level contributes to the ionisinglevels if its energy
      /// is below the ionization potential and the level doesn't
      /// belong to the topmost ion included.
      /// Rate coefficients are only available for ionising levels.
      if (levelenergy < ionpot && ion < nions - 1) ///thats only an option for pure LTE && level < TAKE_N_BFCONTINUA)
      {
        elements[element].ions[ion].ionisinglevels++;

        elements[element].ions[ion].levels[level].cont_index = *cont_index;
        (*cont_index)--;
      }


      /// store the possible downward transitions from the current level in following order to memory
      ///     A_level,level-1; A_level,level-2; ... A_level,1
      /// entries which are not explicitly set are zero (the zero is set/initialized by calloc!)
      if ((transitions[level].to = calloc(level, sizeof(int))) == NULL)
      {
        printout("[fatal] input: not enough memory to initialize transitionlist ... abort\n");
        abort();
      }
      for (int i = 0; i < level; i++)
      {
        transitions[level].to[i] = -99.;
      }
      if ((elements[element].ions[ion].levels[level].downtrans = (transitionlist_entry *) malloc(sizeof(transitionlist_entry))) == NULL)
      {
        printout("[fatal] input: not enough memory to initialize downtranslist ... abort\n");
        abort();
      }
      /// initialize number of downward transitions to zero
      elements[element].ions[ion].levels[level].downtrans[0].targetlevel = 0;
      if ((elements[element].ions[ion].levels[level].uptrans = (transitionlist_entry *) malloc(sizeof(transitionlist_entry))) == NULL)
      {
        printout("[fatal] input: not enough memory to initialize uptranslist ... abort\n");
        abort();
      }
      /// initialize number of upward transitions to zero
      elements[element].ions[ion].levels[level].uptrans[0].targetlevel = 0;
    }
  }
}


static void read_ion_transitions(
  FILE *transitiondata, const int nlevelsmax, const int tottransitions_in,
  int *tottransitions, transitiontable_entry *transitiontable)
{
  if (tottransitions == 0)
  {
    for (int i = 0; i < tottransitions_in; i++)
    {
      double A, coll_str;
      int lower,upper,intforbidden;
      fscanf(transitiondata, "%d %d %lg %lg %d\n", &lower, &upper, &A, &coll_str, &intforbidden);
      //printout("index %d, lower %d, upper %d, A %g\n",transitionindex,lower,upper,A);
    }
  }
  else
  {
    int prev_upper = -1;
    int prev_lower = 0;
    for (int i = 0; i < *tottransitions; i++)
    {
      int lower_in;
      int upper_in;
      double A;
      double coll_str;
      int intforbidden;
      fscanf(transitiondata, "%d %d %lg %lg %d\n", &lower_in, &upper_in, &A, &coll_str, &intforbidden);
      const int lower = lower_in - groundstate_index_in;
      const int upper = upper_in - groundstate_index_in;

      // this entire block can be removed if we don't want to add in extra collisonal
      // transitions between the first n levels
      if (prev_lower < nlevels_requiretransitions && prev_upper < nlevelsmax - 1)
      {
        if (lower == prev_lower && upper > prev_upper + 1) // same lower level, but skipped some upper levels
        {
          printout("+adding transitions from lower level %d starting at lineindex %d\n", prev_lower, i);
          for (int tmplevel = prev_upper + 1; tmplevel < upper; tmplevel++)
          {
            if (tmplevel == prev_lower)
              continue;
            (*tottransitions)++;
            transitiontable = realloc(transitiontable, *tottransitions * sizeof(transitiontable_entry));
            if (transitiontable == NULL)
            {
              printout("Could not reallocate transitiontable\n");
              abort();
            }
            transitiontable[i].lower = prev_lower;
            transitiontable[i].upper = tmplevel;
            transitiontable[i].A = 0.;
            transitiontable[i].coll_str = -2.;
            transitiontable[i].forbidden = true;
            // printout("+adding transition index %d lower %d upper %d\n", i, prev_lower, tmplevel);
            i++;
          }
        }
        else if (lower > prev_lower) // we've moved onto another lower level, but the previous one was missing some required transitions
        {
          for (int tmplevel = prev_upper + 1; tmplevel < nlevelsmax; tmplevel++)
          {
            (*tottransitions)++;
            transitiontable = realloc(transitiontable, *tottransitions * sizeof(transitiontable_entry));
            if (transitiontable == NULL)
            {
              printout("Could not reallocate transitiontable\n");
              abort();
            }
            transitiontable[i].lower = prev_lower;
            transitiontable[i].upper = tmplevel;
            transitiontable[i].A = 0.;
            transitiontable[i].coll_str = -2.;
            transitiontable[i].forbidden = true;
            printout("+adding transition index %d lower %d upper %d\n", i, prev_lower, tmplevel);
            i++;
          }
        }
      }

      transitiontable[i].lower = lower;
      transitiontable[i].upper = upper;
      transitiontable[i].A = A;
      transitiontable[i].coll_str = coll_str;
      transitiontable[i].forbidden = (intforbidden == 1);
      //printout("index %d, lower %d, upper %d, A %g\n",transitionindex,lower,upper,A);
      // printout("reading transition index %d lower %d upper %d\n", i, transitiontable[i].lower, transitiontable[i].upper);
      prev_lower = lower;
      prev_upper = upper;
    }
  }
}


static int compare_linelistentry(const void *restrict p1, const void *restrict p2)
/// Helper function to sort the linelist by frequency.
{
  linelist_entry *restrict a1 = (linelist_entry *)(p1);
  linelist_entry *restrict a2 = (linelist_entry *)(p2);
  //printf("%d %d %d %d %g\n",a1->elementindex,a1->ionindex,a1->lowerlevelindex,a1->upperlevelindex,a1->nu);
  //printf("%d %d %d %d %g\n",a2->elementindex,a2->ionindex,a2->lowerlevelindex,a2->upperlevelindex,a2->nu);
  //printf("%g\n",a2->nu - a1->nu);
  if (fabs(a2->nu - a1->nu) < (1.e-10 * a1->nu))
  {
    a2->nu = a1->nu;
    if (a1->lowerlevelindex > a2->lowerlevelindex)
    {
      return -1;
    }
    else if (a1->lowerlevelindex < a2->lowerlevelindex)
    {
      return 1;
    }
    else if (a1->upperlevelindex > a2->upperlevelindex)
    {
      return -1;
    }
    else if (a1->upperlevelindex < a2->upperlevelindex)
    {
      return 1;
    }
    else
    {
      printout("Duplicate atomic line?\n");
      printout(
        "Z=%d ionstage %d lower %d upper %d nu %g\n",
        get_element(a1->elementindex),
        get_ionstage(a1->elementindex, a1->ionindex),
        a1->lowerlevelindex,
        a1->upperlevelindex,
        a1->nu);
      printout(
        "Z=%d ionstage %d lower %d upper %d nu %g\n",
        get_element(a2->elementindex),
        get_ionstage(a2->elementindex, a2->ionindex),
        a2->lowerlevelindex,
        a2->upperlevelindex,
        a2->nu);
      return 0;
    }
  }
  else
  {
    if ((a1->nu < a2->nu) || (a1->nu == a2->nu))
      return 1;
    else if (a1->nu > a2->nu)
      return -1;
    else
      return 0;
  }
}


static int transitioncheck(int upper, int lower)
{
  const int index = (upper - lower) - 1;
  const int flag = transitions[upper].to[index];

  return flag;
}


static void add_transitions_to_linelist(
  const int element, const int ion, const int nlevelsmax, const int tottransitions,
  transitiontable_entry *transitiontable, int *nuparr, int *ndownarr, int *lineindex)
{
  for (int ii = 0; ii < tottransitions; ii++)
  {
    const int level = transitiontable[ii].upper;
    const double epsilon_upper = epsilon(element, ion, level);
    if (level < nlevelsmax)
    {
      const int targetlevel = transitiontable[ii].lower;
      const double epsilon_lower = epsilon(element, ion, targetlevel);

      //if (level == transitiontable[ii].upper && level-i-1 == transitiontable[ii].lower)
      //{
      //printout("ii %d\n",ii);
      //printout("transtable upper %d, lower %d, A %g, iii %d\n",transitiontable[ii].upper,transitiontable[ii].lower, transitiontable[ii].A,iii);
      /// Make sure that we don't allow duplicate. In that case take only the lines
      /// first occurrence
      if (transitioncheck(level,targetlevel) == -99)
      {
        transitions[level].to[level-targetlevel-1] = *lineindex;
        const double A_ul = transitiontable[ii].A;
        const double coll_str = transitiontable[ii].coll_str;
        //elements[element].ions[ion].levels[level].transitions[level-targetlevel-1].einstein_A = A_ul;

        const double nu_trans = (epsilon(element,ion,level) - epsilon(element,ion,targetlevel)) / H;
        const double g = stat_weight(element,ion,level) / stat_weight(element,ion,targetlevel);
        const double f_ul = g * ME * pow(CLIGHT,3) / (8 * pow(QE * nu_trans * PI, 2)) * A_ul;
        //f_ul = g * OSCSTRENGTHCONVERSION / pow(nu_trans,2) * A_ul;
        //elements[element].ions[ion].levels[level].transitions[level-targetlevel-1].oscillator_strength = g * ME*pow(CLIGHT,3)/(8*pow(QE*nu_trans*PI,2)) * A_ul;

        //printout("lineindex %d, element %d, ion %d, lower %d, upper %d, nu %g\n",*lineindex,element,ion,level-i-1,level,nu_trans);
        linelist[*lineindex].elementindex = element;
        linelist[*lineindex].ionindex = ion;
        linelist[*lineindex].lowerlevelindex = targetlevel;
        linelist[*lineindex].upperlevelindex = level;
        linelist[*lineindex].nu = nu_trans;
        linelist[*lineindex].einstein_A = A_ul;
        linelist[*lineindex].osc_strength = f_ul;
        linelist[*lineindex].coll_str = coll_str;
        linelist[*lineindex].forbidden = transitiontable[ii].forbidden;
        (*lineindex)++;
        if (*lineindex % MLINES == 0)
        {
          printout("[info] read_atomicdata: increase linelistsize from %d to %d\n", *lineindex, *lineindex + MLINES);
          if ((linelist = realloc(linelist, (*lineindex + MLINES) * sizeof(linelist_entry))) == NULL)
          {
            printout("[fatal] input: not enough memory to reallocate linelist ... abort\n");
            abort();
          }
        }

        /// This is not a metastable level.
        elements[element].ions[ion].levels[level].metastable = false;

        elements[element].ions[ion].levels[level].downtrans[0].targetlevel = ndownarr[level];
        if ((elements[element].ions[ion].levels[level].downtrans
            = realloc(elements[element].ions[ion].levels[level].downtrans, (ndownarr[level] + 1) * sizeof(transitionlist_entry))) == NULL)
        {
          printout("[fatal] input: not enough memory to reallocate downtranslist ... abort\n");
          abort();
        }
        elements[element].ions[ion].levels[level].downtrans[ndownarr[level]].targetlevel = targetlevel;
        elements[element].ions[ion].levels[level].downtrans[ndownarr[level]].epsilon_trans = epsilon_upper - epsilon_lower;
        //elements[element].ions[ion].levels[level].downtrans[ndownarr[level]].einstein_A = A_ul;
        //elements[element].ions[ion].levels[level].downtrans[ndownarr[level]].oscillator_strength = f_ul;
        ndownarr[level]++;

        elements[element].ions[ion].levels[targetlevel].uptrans[0].targetlevel = nuparr[targetlevel];
        if ((elements[element].ions[ion].levels[targetlevel].uptrans
            = realloc(elements[element].ions[ion].levels[targetlevel].uptrans, (nuparr[targetlevel] + 1) * sizeof(transitionlist_entry))) == NULL)
        {
          printout("[fatal] input: not enough memory to reallocate uptranslist ... abort\n");
          abort();
        }
        elements[element].ions[ion].levels[targetlevel].uptrans[nuparr[targetlevel]].targetlevel = level;
        elements[element].ions[ion].levels[targetlevel].uptrans[nuparr[targetlevel]].epsilon_trans = epsilon_upper - epsilon_lower;
        nuparr[targetlevel]++;
      }
      else
      {
        // This is a new branch to deal with lines that have different types of transition. It should trip after a transition is already known.
        const int linelistindex = transitions[level].to[level - targetlevel - 1];
        const double A_ul = transitiontable[ii].A;
        const double coll_str = transitiontable[ii].coll_str;
        //elements[element].ions[ion].levels[level].transitions[level-targetlevel-1].einstein_A = A_ul;

        const double nu_trans = (epsilon(element, ion, level) - epsilon(element, ion, targetlevel)) / H;
        const double g = stat_weight(element, ion, level) / stat_weight(element, ion, targetlevel);
        const double f_ul = g * ME * pow(CLIGHT,3) / (8 * pow(QE * nu_trans * PI, 2)) * A_ul;
        //f_ul = g * OSCSTRENGTHCONVERSION / pow(nu_trans,2) * A_ul;
        //elements[element].ions[ion].levels[level].transitions[level-targetlevel-1].oscillator_strength = g * ME*pow(CLIGHT,3)/(8*pow(QE*nu_trans*PI,2)) * A_ul;

        if ((linelist[linelistindex].elementindex != element) || (linelist[linelistindex].ionindex != ion) || (linelist[linelistindex].upperlevelindex != level) || (linelist[linelistindex].lowerlevelindex != targetlevel))
        {
          printout("[input.c] Failure to identify level pair for duplicate bb-transition ... going to abort now\n");
          printout("[input.c]   element %d ion %d targetlevel %d level %d\n", element, ion, targetlevel, level);
          printout("[input.c]   transitions[level].to[level-targetlevel-1]=linelistindex %d\n", transitions[level].to[level - targetlevel - 1]);
          printout("[input.c]   A_ul %g, coll_str %g\n", A_ul, coll_str);
          printout("[input.c]   linelist[linelistindex].elementindex %d, linelist[linelistindex].ionindex %d, linelist[linelistindex].upperlevelindex %d, linelist[linelistindex].lowerlevelindex %d\n", linelist[linelistindex].elementindex, linelist[linelistindex].ionindex, linelist[linelistindex].upperlevelindex,linelist[linelistindex].lowerlevelindex);
          abort();
        }
        linelist[linelistindex].einstein_A += A_ul;
        linelist[linelistindex].osc_strength += f_ul;
        if (coll_str > linelist[linelistindex].coll_str)
        {
          linelist[linelistindex].coll_str = coll_str;
        }
      }
    }
  }
}


static void read_atomicdata_files(void)
{
  int totaluptrans = 0;
  int totaldowntrans = 0;
  int cont_index = -1;

  ///open atomic data file
  FILE *restrict compositiondata = fopen("compositiondata.txt", "r");
  if (compositiondata == NULL)
  {
    printout("Cannot open compositiondata.txt.\n");
    abort();
  }

  FILE *restrict adata = fopen("adata.txt", "r");
  if (adata == NULL)
  {
    printout("Cannot open adata.txt.\n");
    abort();
  }

  /// initialize atomic data structure to number of elements
  fscanf(compositiondata,"%d",&nelements);
  if ((elements = calloc(nelements, sizeof(elementlist_entry))) == NULL)
  {
    printout("[fatal] input: not enough memory to initialize elementlist ... abort\n");
    abort();
  }
  //printout("elements initialized\n");

  /// Initialize the linelist
  if ((linelist = calloc(MLINES, sizeof(linelist_entry))) == NULL)
  {
    printout("[fatal] input: not enough memory to initialize linelist ... abort\n");
    abort();
  }

  /// temperature to determine relevant ionstages
  int T_preset;
  fscanf(compositiondata,"%d",&T_preset);
  int homogeneous_abundances_in;
  fscanf(compositiondata,"%d",&homogeneous_abundances_in);
  homogeneous_abundances = (homogeneous_abundances_in != 0);
  if (homogeneous_abundances)
    printout("[info] read_atomicdata: homogeneous abundances as defined in compositiondata.txt are active\n");

  /// open transition data file
  FILE *restrict transitiondata;
  if ((transitiondata = fopen("transitiondata.txt", "r")) == NULL)
  {
    printout("Cannot open transitiondata.txt.\n");
    abort();
  }
  int lineindex = 0;  ///counter to determine the total number of lines, initialisation

  /// readin
  int nbfcheck = 0;
  int heatingcheck = 0;
  int coolingcheck = 0;
  for (int element = 0; element < nelements; element++)
  {
    /// read information about the next element which should be stored to memory
    int Z, nions, lowermost_ionstage, uppermost_ionstage, nlevelsmax_readin;
    double abundance, mass_amu;
    fscanf(compositiondata,"%d %d %d %d %d %lg %lg", &Z, &nions, &lowermost_ionstage, &uppermost_ionstage, &nlevelsmax_readin, &abundance, &mass_amu);
    printout("readin compositiondata: next element Z %d, nions %d, lowermost %d, uppermost %d, nlevelsmax %d\n",Z,nions,lowermost_ionstage,uppermost_ionstage,nlevelsmax_readin);
    assert(Z > 0);
    assert(nions > 0);
    assert(nions == uppermost_ionstage - lowermost_ionstage + 1);
    assert(abundance >= 0);
    assert(mass_amu >= 0);

    /// write this element's data to memory
    elements[element].anumber = Z;
    elements[element].nions = nions;
    elements[element].abundance = abundance;       /// abundances are expected to be given by mass
    elements[element].mass = mass_amu * MH;
    includedions += nions;

    /// Initialize the elements ionlist
    if ((elements[element].ions = (ionlist_entry *) calloc(nions, sizeof(ionlist_entry))) == NULL)
    {
        printout("[fatal] input: not enough memory to initialize ionlist ... abort\n");
        abort();
    }

    /// now read in data for all ions of the current element. before doing so initialize
    /// energy scale for the current element (all level energies are stored relative to
    /// the ground level of the neutral ion)
    double energyoffset = 0.;
    double ionpot = 0.;
    for (int ion = 0; ion < nions; ion++)
    {
      int nlevelsmax = nlevelsmax_readin;
      printout("element %d ion %d\n", element, ion);
      /// calculate the current levels ground level energy
      energyoffset += ionpot;

      /// read information for the elements next ionstage
      int adata_Z_in = -1;
      int ionstage = -1;
      int nlevels = 0;
      // fscanf(adata,"%d %d %d %lg\n",&adata_Z_in,&ionstage,&nlevels,&ionpot);
      while (adata_Z_in != Z || ionstage != lowermost_ionstage + ion) // skip over this ion block
      {
        if (adata_Z_in == Z)
        {
          printout("increasing energyoffset by ionpot %g\n", ionpot);
          energyoffset += ionpot;
        }
        for (int i = 0; i < nlevels; i++)
        {
          double levelenergy;
          double statweight;
          int levelindex;
          int ntransitions;
          fscanf(adata, "%d %lg %lg %d%*[^\n]\n", &levelindex, &levelenergy, &statweight, &ntransitions);
        }
        const int fscanfadata = fscanf(adata, "%d %d %d %lg\n", &adata_Z_in, &ionstage, &nlevels, &ionpot);

        if (fscanfadata == EOF)
        {
          printout("End of file in adata not expected");
          abort();
        }
      }

      printout("adata header matched: Z %d, ionstage %d, nlevels %d\n", adata_Z_in, ionstage, nlevels);
      if (nlevelsmax < 0)
      {
        nlevelsmax = nlevels;
      }
      else if (nlevels >= nlevelsmax)
      {
        printout("[info] read_atomicdata: reduce number of levels from %d to %d for ion %d of element %d\n", nlevels, nlevelsmax, ion, element);
      }
      else
      {
        printout("[warning] read_atomicdata: requested nlevelsmax=%d > nlevels=%d for ion %d of element %d ... reduced nlevelsmax to nlevels\n",
                 nlevelsmax, nlevels, ion, element);
        nlevelsmax = nlevels;
      }


      /// and proceed through the transitionlist till we match this ionstage (if it was not the neutral one)
      int transdata_Z_in = -1;
      int transdata_ionstage_in = -1;
      int tottransitions_in = 0;
      while (transdata_Z_in != Z || transdata_ionstage_in != ionstage)
      {
        for (int i = 0; i < tottransitions_in; i++)
        {
          int lower;
          int upper;
          double A;
          double coll_str;
          int intforbidden;
          fscanf(transitiondata, "%d %d %lg %lg %d\n", &lower, &upper, &A, &coll_str, &intforbidden);
        }
        const int readtransdata = fscanf(transitiondata,"%d %d %d", &transdata_Z_in, &transdata_ionstage_in, &tottransitions_in);
        if (readtransdata == EOF)
        {
          printout("End of file in transition data");
          abort();
        }
      }

      printout("transdata header matched: transdata_Z_in %d, transdata_ionstage_in %d, tottransitions %d\n",
               transdata_Z_in, transdata_ionstage_in, tottransitions_in);
      assert(tottransitions_in >= 0);

      int tottransitions = tottransitions_in;

      if (ion == nions - 1) // limit the top ion to one level and no transitions
      {
        nlevelsmax = 1;
        tottransitions = 0;
      }

      /// then read in its level and transition data
      if (transdata_Z_in == Z && transdata_ionstage_in == ionstage)
      {
        transitiontable_entry *transitiontable = calloc(tottransitions, sizeof(transitiontable_entry));

        /// load transition table for the CURRENT ion to temporary memory
        if (transitiontable == NULL)
        {
          if (tottransitions > 0)
          {
            printout("[fatal] input: not enough memory to initialize transitiontable ... abort\n");
            abort();
          }
        }

        read_ion_transitions(transitiondata, nlevelsmax, tottransitions_in, &tottransitions, transitiontable);

        /// store the ions data to memory and set up the ions zeta and levellist
        elements[element].ions[ion].ionstage = ionstage;
        elements[element].ions[ion].nlevels = nlevelsmax;
        elements[element].ions[ion].ionisinglevels = 0;
        elements[element].ions[ion].ionpot = ionpot * EV;
//           if ((elements[element].ions[ion].zeta = calloc(TABLESIZE, sizeof(float))) == NULL)
//           {
//             printout("[fatal] input: not enough memory to initialize zetalist for element %d, ion %d ... abort\n",element,ion);
//             abort();
//           }
        if ((elements[element].ions[ion].levels = calloc(nlevelsmax, sizeof(levellist_entry))) == NULL)
        {
          printout("[fatal] input: not enough memory to initialize levelist of element %d, ion %d ... abort\n",element,ion);
          abort();
        }


        /// now we need to readout the data for all those levels, write them to memory
        /// and set up the list of possible transitions for each level
        int nuparr[nlevelsmax];
        int ndownarr[nlevelsmax];
        if ((transitions = calloc(nlevelsmax, sizeof(transitions_t))) == NULL)
        {
          printout("[fatal] input: not enough memory to allocate transitions ... abort\n");
          abort();
        }

        read_ion_levels(adata, element, ion, nions, nlevels, nlevelsmax, energyoffset, ionpot, &cont_index, nuparr, ndownarr);

        add_transitions_to_linelist(element, ion, nlevelsmax, tottransitions, transitiontable, nuparr, ndownarr, &lineindex);

        //printf("A %g\n",elements[element].ions[ion].levels[level].transitions[i].einstein_A );
        //printout("%d -> %d has A %g\n",level,level-i-1,elements[element].ions[ion].levels[level].transitions[i].einstein_A );

        for (int level = 0; level < nlevelsmax; level++)
        {
          totaldowntrans += elements[element].ions[ion].levels[level].downtrans[0].targetlevel;
          totaluptrans += elements[element].ions[ion].levels[level].uptrans[0].targetlevel;
          free(transitions[level].to);
        }
        free(transitiontable);
        free(transitions);

        /// Also the phixslist
        if (ion < nions - 1)
        {
//            elements[element].ions[ion].nbfcontinua = elements[element].ions[ion].ionisinglevels;//nlevelsmax;
          nbfcheck += elements[element].ions[ion].ionisinglevels; //nlevelsmax;
/*            if ((elements[element].ions[ion].phixslist = calloc(nlevelsmax, sizeof(ionsphixslist_t))) == NULL)
          {
            printout("[fatal] input: not enough memory to initialize phixslist for element %d, ion %d ... abort\n",element,ion);
            abort();
          }*/
          nbfcontinua += get_bfcontinua(element,ion);//elements[element].ions[ion].ionisinglevels;//nlevelsmax;
        }
      }
      else
      {
        printout("corrupted atomic data\n");
        abort();
      }
    }
  }
  fclose(adata);
  fclose(transitiondata);
  fclose(compositiondata);
  printout("nbfcheck %d\n",nbfcheck);
  printout("heatingcheck %d\n",heatingcheck);

  /// Save the linecounters value to the global variable containing the number of lines
  nlines = lineindex;
  printout("nlines %d\n",nlines);
  if (nlines > 0)
  {
    /// and release empty memory from the linelist
    if ((linelist = realloc(linelist, nlines * sizeof(linelist_entry))) == NULL)
    {
      printout("[fatal] input: not enough memory to reallocate linelist ... abort\n");
      abort();
    }
  }

  if (T_preset > 0)
    abort();


  /// Set up the list of allowed upward transitions for each level
  printout("total uptrans %d\n", totaluptrans);
  printout("total downtrans %d\n", totaldowntrans);
  printout("coolingcheck %d\n", coolingcheck);


  ///debug output
  /*
  FILE *linelist_file;
  if ((linelist_file = fopen("linelist_unsorted.out", "w")) == NULL){
    printout("Cannot open linelist.out.\n");
    abort();
  }
  for (i = 0; i < nlines; i++)
    fprintf(linelist_file,"element %d, ion %d, ionstage %d, upperindex %d, lowerindex %d, nu %g\n",linelist[i].elementindex, linelist[i].ionindex, elements[linelist[i].elementindex].ions[linelist[i].ionindex].ionstage, linelist[i].upperlevelindex, linelist[i].lowerlevelindex, linelist[i].nu);
  fclose(linelist_file);
  //abort();
  */

  /// then sort the linelist by decreasing frequency
  qsort(linelist, nlines, sizeof(linelist_entry), compare_linelistentry);

  /// Save sorted linelist into a file
  // if (rank_global == 0)
  // {
  //   FILE *restrict linelist_file;
  //   if ((linelist_file = fopen("linelist.dat", "w")) == NULL)
  //   {
  //     printout("Cannot open linelist.out.\n");
  //     abort();
  //   }
  //   fprintf(linelist_file,"%d\n",nlines);
  //   for (int i = 0; i < nlines; i++)
  //   {
  //     fprintf(linelist_file,"%d %d %d %d %d %lg %lg %lg %lg %d\n",
  //             i, linelist[i].elementindex, linelist[i].ionindex,
  //             linelist[i].upperlevelindex, linelist[i].lowerlevelindex,
  //             linelist[i].nu, linelist[i].einstein_A, linelist[i].osc_strength,
  //             linelist[i].coll_str, linelist[i].forbidden);
  //   }
  //   fclose(linelist_file);
  // }


  ///Establish connection between transitions and sorted linelist
  //printout("[debug] init line counter list\n");
  printout("establish connection between transitions and sorted linelist\n");
  for (int i = 0; i < nlines; i++)
  {
    const int element = linelist[i].elementindex;
    const int ion = linelist[i].ionindex;
    const int lowerlevel = linelist[i].lowerlevelindex;
    const int upperlevel = linelist[i].upperlevelindex;
    for (int ii = 1; ii <= elements[element].ions[ion].levels[upperlevel].downtrans[0].targetlevel; ii++)
    {
      if (elements[element].ions[ion].levels[upperlevel].downtrans[ii].targetlevel == lowerlevel)
      {
        elements[element].ions[ion].levels[upperlevel].downtrans[ii].lineindex = i;
      }
    }
    for (int ii = 1; ii <= elements[element].ions[ion].levels[lowerlevel].uptrans[0].targetlevel; ii++)
    {
      if (elements[element].ions[ion].levels[lowerlevel].uptrans[ii].targetlevel == upperlevel)
      {
        elements[element].ions[ion].levels[lowerlevel].uptrans[ii].lineindex = i;
      }
    }
  }


  /// Photoionisation cross-sections
  ///======================================================
  ///finally read in photoionisation cross sections and store them to the atomic data structure
  read_phixs_data();

  printout("cont_index %d\n",cont_index);
}


static int search_groundphixslist(double nu_edge, int *index_in_groundlevelcontestimator, int el, int in, int ll)
/// Return the closest ground level continuum index to the given edge
/// frequency. If the given edge frequency is redder than the reddest
/// continuum return -1.
/// NB: groundphixslist must be in ascending order.
{
  int index;

  if (nu_edge < phixslist[tid].groundcont[0].nu_edge)
  {
    index = -1;
    *index_in_groundlevelcontestimator = -1;
  }
  else
  {
    int i,element,ion;
    for (i = 1; i < nbfcontinua_ground; i++)
    {
      if (nu_edge < phixslist[tid].groundcont[i].nu_edge) break;
    }
/*    if (i == nbfcontinua_ground)
    {
      printout("[fatal] search_groundphixslist: i %d, nu_edge %g, phixslist[tid].groundcont[i-1].nu_edge %g ... abort\n",i,nu_edge,phixslist[tid].groundcont[i-1].nu_edge);
      printout("[fatal] search_groundphixslist: this is element %d, ion %d, level %d in groundphixslist at i-1\n",el,in,ll);
      //printout("[fatal] search_groundphixslist: this is element %d, ion %d, level %d in groundphixslist at i-1\n",phixslist[tid].groundcont[i-1].element,phixslist[tid].groundcont[i-1].ion,groundphixslist[i-1].level);
      abort();
    }*/
    if (i == nbfcontinua_ground)
    {
      element = phixslist[tid].groundcont[i-1].element;
      ion = phixslist[tid].groundcont[i-1].ion;
      int level = phixslist[tid].groundcont[i-1].level;
      if (element == el && ion == in && level == ll)
      {
        index = i - 1;
      }
      else
      {
        printout("[fatal] search_groundphixslist: element %d, ion %d, level %d has edge_frequency %g equal to the bluest ground-level continuum\n",el,in,ll,nu_edge);
        printout("[fatal] search_groundphixslist: bluest ground level continuum is element %d, ion %d, level %d at nu_edge %g\n",element,ion,level,phixslist[tid].groundcont[i-1].nu_edge);
        printout("[fatal] search_groundphixslist: i %d, nbfcontinua_ground %d\n",i,nbfcontinua_ground);
        printout("[fatal] This shouldn't happen, is hoewever possible if there are multiple levels in the adata file at energy=0\n");
        for (int looplevels = 0; looplevels < get_nlevels(el,in); looplevels++)
        {
          printout("[fatal]   element %d, ion %d, level %d, energy %g\n",el,in,looplevels,epsilon(el,in,looplevels));
        }
        printout("[fatal] Abort omitted ... MAKE SURE ATOMIC DATA ARE CONSISTENT\n");
        index = i-1;
        //abort();
      }
    }
    else
    {
      double left = nu_edge - phixslist[tid].groundcont[i-1].nu_edge;
      double right = phixslist[tid].groundcont[i].nu_edge - nu_edge;
      if (left <= right)
      {
        element = phixslist[tid].groundcont[i-1].element;
        ion = phixslist[tid].groundcont[i-1].ion;
        index = i-1;
      }
      else
      {
        element = phixslist[tid].groundcont[i].element;
        ion = phixslist[tid].groundcont[i].ion;
        index = i;
      }
    }
    *index_in_groundlevelcontestimator = element*maxion + ion;
  }

  return index;
}


static void setup_cellhistory(void)
{
  /// SET UP THE CELL HISTORY
  ///======================================================
  /// Stack which holds information about population and other cell specific data
  /// ===> move to update_packets
  if ((cellhistory = (cellhistory_struct *) malloc(nthreads*sizeof(cellhistory_struct))) == NULL)
  {
    printout("[fatal] input: not enough memory to initialize cellhistory of size %d... abort\n",nthreads);
    abort();
  }
  #ifdef _OPENMP
    #pragma omp parallel
    {
  #endif
      printout("[info] input: initializing cellhistory for thread %d ...\n",tid);

      cellhistory[tid].cellnumber = -99;

      cellhistory[tid].coolinglist = (cellhistorycoolinglist_t *) malloc(ncoolingterms*sizeof(cellhistorycoolinglist_t));
      if (cellhistory[tid].coolinglist == NULL)
      {
        printout("[fatal] input: not enough memory to initialize cellhistory's coolinglist ... abort\n");
        abort();
      }

      if ((cellhistory[tid].chelements = (chelements_struct *) malloc(nelements*sizeof(chelements_struct))) == NULL)
      {
        printout("[fatal] input: not enough memory to initialize cellhistory's elementlist ... abort\n");
        abort();
      }
      for (int element = 0; element < nelements; element++)
      {
        const int nions = get_nions(element);
        if ((cellhistory[tid].chelements[element].chions = (chions_struct *) malloc(nions*sizeof(chions_struct))) == NULL)
        {
          printout("[fatal] input: not enough memory to initialize cellhistory's ionlist ... abort\n");
          abort();
        }
        for (int ion = 0; ion < nions; ion++)
        {
          const int nlevels = get_nlevels(element,ion);
          if ((cellhistory[tid].chelements[element].chions[ion].chlevels = (chlevels_struct *) malloc(nlevels*sizeof(chlevels_struct))) == NULL)
          {
            printout("[fatal] input: not enough memory to initialize cellhistory's levellist ... abort\n");
            abort();
          }
          for (int level = 0; level < nlevels; level++)
          {
            const int nphixstargets = get_nphixstargets(element,ion,level);
            if ((cellhistory[tid].chelements[element].chions[ion].chlevels[level].chphixstargets = (chphixstargets_struct *) malloc(nphixstargets*sizeof(chphixstargets_struct))) == NULL)
            {
              printout("[fatal] input: not enough memory to initialize cellhistory's chphixstargets ... abort\n");
              abort();
            }

            int ndowntrans = elements[element].ions[ion].levels[level].downtrans[0].targetlevel;
            int nuptrans = elements[element].ions[ion].levels[level].uptrans[0].targetlevel;

            if ((cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_rad_deexc = (double *) malloc((ndowntrans+1)*sizeof(double))) == NULL)
            {
              printout("[fatal] input: not enough memory to initialize cellhistory's individ_rad_deexc ... abort\n");
              abort();
            }
            cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_rad_deexc[0] = ndowntrans;

            if ((cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_internal_down_same = (double *) malloc((ndowntrans+1)*sizeof(double))) == NULL)
            {
              printout("[fatal] input: not enough memory to initialize cellhistory's individ_internal_down_same ... abort\n");
              abort();
            }
            cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_internal_down_same[0] = ndowntrans;

            if ((cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_internal_up_same = (double *) malloc((nuptrans+1)*sizeof(double))) == NULL)
            {
              printout("[fatal] input: not enough memory to initialize cellhistory's individ_internal_up_same ... abort\n");
              abort();
            }
            cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_internal_up_same[0] = nuptrans;
          }
        }
      }
  #ifdef _OPENMP
    }
  #endif
}


static void write_bflist_file(int includedionisinglevels)
{
  FILE *bflist_file;
  if (rank_global == 0)
  {
    if ((bflist_file = fopen("bflist.dat", "w")) == NULL){
      printout("Cannot open bflist.dat.\n");
      abort();
    }
    fprintf(bflist_file,"%d\n",includedionisinglevels);
  }
  int i = 0;
  for (int element = 0; element < nelements; element++)
  {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++)
    {
      const int nlevels = get_ionisinglevels(element,ion);
      for (int level = 0; level < nlevels; level++)
      {
        bflist[i].elementindex = element;
        bflist[i].ionindex = ion;
        bflist[i].levelindex = level;
        if (rank_global == 0)
          fprintf(bflist_file,"%d %d %d %d\n",i,element,ion,level);
        i++;
      }
    }
  }
  if (rank_global == 0)
    fclose(bflist_file);
}


static void setup_coolinglist(void)
{
  /// SET UP THE COOLING LIST
  ///======================================================
  /// Determine number of processes which allow kpkts to convert to something else.
  /// This number is given by the collisional excitations (so far determined from the oscillator strengths
  /// by the van Regemorter formula, therefore totaluptrans), the number of free-bound emissions and collisional ionisations
  /// (as long as we only deal with ionisation to the ground level this means for both of these
  /// \sum_{elements,ions}get_nlevels(element,ion) and free-free which is \sum_{elements} get_nions(element)-1
  /*ncoolingterms = totaluptrans;
  for (element = 0; element < nelements; element++)
  {
    nions = get_nions(element);
    for (ion=0; ion < nions; ion++)
    {
      if (get_ionstage(element,ion) > 1) ncoolingterms++;
      if (ion < nions - 1) ncoolingterms += 2 * get_ionisinglevels(element,ion);
    }
  }
  printout("[info] read_atomicdata: number of coolingterms %d\n",ncoolingterms);*/

  ncoolingterms = 0;
  for (int element = 0; element < nelements; element++)
  {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++)
    {
      int add = 0; /// Helper variable to count coolingterms per ion
      elements[element].ions[ion].coolingoffset = ncoolingterms;
      /// Ionised ions add one ff-cooling term
      if (get_ionstage(element,ion) > 1)
        add++;
      /// Ionisinglevels below the closure ion add to bf and col ionisation
      if (ion < nions - 1)
      //  add += 2 * get_ionisinglevels(element,ion);
        add += 2 * get_bfcontinua(element,ion);
      /// All the levels add number of col excitations
      const int nlevels = get_nlevels(element,ion);
      for (int level = 0; level < nlevels; level++)
      {
        add += elements[element].ions[ion].levels[level].uptrans[0].targetlevel;
        //if (ion < nions - 1) and (level < get_ionisinglevels(element,ion))
        //  add += get_nphixstargets(element,ion,level)
      }
      elements[element].ions[ion].ncoolingterms = add;
      ncoolingterms += add;
    }
  }
  printout("[info] read_atomicdata: number of coolingterms %d\n",ncoolingterms);

  /// And setup the global coolinglist. A shorter coolinglist with only the important contributors
  /// is part of the cellhistory.
  /*
  if ((globalcoolinglist = malloc(ncoolingterms*sizeof(coolinglist_entry))) == NULL)
  {
    printout("[fatal] read_atomicdata: not enough memory to initialize globalcoolinglist ... abort\n");
    abort();
  }
  */
}


static int compare_phixslistentry_bynuedge(const void *restrict p1, const void *restrict p2)
/// Helper function to sort the phixslist by ascending threshold frequency.
{
  const fullphixslist_t *restrict a1 = (fullphixslist_t *)(p1);
  const fullphixslist_t *restrict a2 = (fullphixslist_t *)(p2);

  double edge_diff = a1->nu_edge - a2->nu_edge;
  if (edge_diff < 0)
    return -1;
  else if (edge_diff > 0)
    return 1;
  else
    return 0;
}


static int compare_groundphixslistentry_bynuedge(const void *restrict p1, const void *restrict p2)
/// Helper function to sort the groundphixslist by ascending threshold frequency.
{
  const groundphixslist_t *restrict a1 = (groundphixslist_t *)(p1);
  const groundphixslist_t *restrict a2 = (groundphixslist_t *)(p2);

  double edge_diff = a1->nu_edge - a2->nu_edge;
  if (edge_diff < 0)
    return -1;
  else if (edge_diff > 0)
    return 1;
  else
    return 0;
}


static void setup_phixs_list(void)
{
  /// SET UP THE PHIXSLIST
  ///======================================================
  printout("[info] read_atomicdata: number of bfcontinua %d\n",nbfcontinua);
  nbfcontinua_ground = includedions - nelements;
  printout("[info] read_atomicdata: number of ground-level bfcontinua %d\n",nbfcontinua_ground);

  phixslist = (phixslist_t *) malloc(nthreads * sizeof(phixslist_t));
  if (phixslist == NULL)
  {
    printout("[fatal] read_atomicdata: not enough memory to initialize phixslist... abort\n");
    abort();
  }

  /// MK: 2012-01-19
  /// To fix the OpenMP problem on BlueGene machines this parallel section was removed and replaced by
  /// a serial loop which intializes the phixslist data structure for all threads in a loop. I'm still
  /// not sure why this causes a problem at all and on BlueGene architectures in particular. However,
  /// it seems to fix the problem.
  //#ifdef _OPENMP
  //  #pragma omp parallel private(i,element,ion,level,nions,nlevels,epsilon_upper,E_threshold,nu_edge)
  //  {
  //#endif
  for (int itid = 0; itid < nthreads; itid++)
  {
    /// Number of ground level bf-continua equals the total number of included ions minus the number
    /// of included elements, because the uppermost ionisation stages can't ionise.
    //printout("groundphixslist nbfcontinua_ground %d\n",nbfcontinua_ground);
    printout("initialising groundphixslist for itid %d\n",itid);
    phixslist[itid].groundcont = (groundphixslist_t *) malloc(nbfcontinua_ground*sizeof(groundphixslist_t));
    if (phixslist[itid].groundcont == NULL)
    {
      printout("[fatal] read_atomicdata: not enough memory to initialize phixslist[%d].groundcont... abort\n",itid);
      abort();
    }

    int i = 0;
    for (int element = 0; element < nelements; element++)
    {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions-1; ion++)
      {
        const int level = 0;
        // const int upperlevel = get_phixsupperlevel(element, ion, level, 0);
        // const double E_threshold = epsilon(element, ion + 1, upperlevel) - epsilon(element, ion, level);
        const double E_threshold = get_phixs_threshold(element, ion, level, 0);
        const double nu_edge = E_threshold / H;
        phixslist[itid].groundcont[i].element = element;
        phixslist[itid].groundcont[i].ion = ion;
        phixslist[itid].groundcont[i].level = level;
        phixslist[itid].groundcont[i].nu_edge = nu_edge;
        //printout("phixslist.groundcont nbfcontinua_ground %d, i %d, element %d, ion %d, level %d, nu_edge %g\n",nbfcontinua_ground,i,element,ion,level,nu_edge);
        i++;
      }
    }
    qsort(phixslist[itid].groundcont,nbfcontinua_ground,sizeof(groundphixslist_t),compare_groundphixslistentry_bynuedge);


    //if (TAKE_N_BFCONTINUA >= 0) phixslist = malloc(includedions*TAKE_N_BFCONTINUA*sizeof(phixslist_t));
    //else
    phixslist[itid].allcont = (fullphixslist_t *) malloc(nbfcontinua*sizeof(fullphixslist_t));
    if (phixslist[itid].allcont == NULL)
    {
      printout("[fatal] read_atomicdata: not enough memory to initialize phixslist... abort\n");
      abort();
    }

    i = 0;
    for (int element = 0; element < nelements; element++)
    {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions - 1; ion++)
      {
        const int nlevels = get_bfcontinua(element,ion);
        //nlevels = get_ionisinglevels(element,ion);
        ///// The following line reduces the number of bf-continua per ion
        //if (nlevels > TAKE_N_BFCONTINUA) nlevels = TAKE_N_BFCONTINUA;
        for (int level = 0; level < nlevels; level++)
        {
          // const int upperlevel = get_phixsupperlevel(element, ion,level, 0);
          // const double E_threshold = epsilon(element, ion + 1, upperlevel) - epsilon(element, ion, level);
          const double E_threshold = get_phixs_threshold(element, ion, level, 0);
          const double nu_edge = E_threshold / H;

          int index_in_groundlevelcontestimator;

          phixslist[itid].allcont[i].element = element;
          phixslist[itid].allcont[i].ion = ion;
          phixslist[itid].allcont[i].level = level;
          phixslist[itid].allcont[i].nu_edge = nu_edge;
          phixslist[itid].allcont[i].index_in_groundphixslist = search_groundphixslist(nu_edge,&index_in_groundlevelcontestimator,element,ion,level);
          #if (!NO_LUT_PHOTOION || !NO_LUT_BFHEATING)
            if (itid == 0)
              elements[element].ions[ion].levels[level].closestgroundlevelcont = index_in_groundlevelcontestimator;
          #endif
          i++;
        }
      }
    }
    //nbfcontinua = i;
    //printout("number of bf-continua reduced to %d\n",nbfcontinua);
    qsort(phixslist[itid].allcont,nbfcontinua,sizeof(fullphixslist_t),compare_phixslistentry_bynuedge);
  }
  //#ifdef _OPENMP
  //  }
  //#endif
}


static void read_atomicdata(void)
/// Subroutine to read in input parameters.
{
  ///new atomic data scheme by readin of adata////////////////////////////////////////////////////////////////////////
  nbfcontinua = 0;
  includedions = 0;

  read_atomicdata_files();
  last_phixs_nuovernuedge = (1.0 + NPHIXSNUINCREMENT * (NPHIXSPOINTS - 1));

  printout("included ions %d\n",includedions);

  /// INITIALISE THE ABSORPTION/EMISSION COUNTERS ARRAYS
  ///======================================================
  #ifdef RECORD_LINESTAT
    if ((ecounter  = (int *) malloc(nlines*sizeof(int))) == NULL)
    {
      printout("[fatal] input: not enough memory to initialise ecounter array ... abort\n");
      abort();
    }
    if ((acounter  = (int *) malloc(nlines*sizeof(int))) == NULL)
    {
      printout("[fatal] input: not enough memory to initialise ecounter array ... abort\n");
      abort();
    }
    if ((linestat_reduced  = (int *) malloc(nlines*sizeof(int))) == NULL)
    {
      printout("[fatal] input: not enough memory to initialise ecounter array ... abort\n");
      abort();
    }
  #endif

  setup_coolinglist();

  setup_cellhistory();

  /// Printout some information about the read-in model atom
  ///======================================================
  //includedions = 0;
  int includedlevels = 0;
  int includedionisinglevels = 0;
  int includedphotoiontransitions = 0;
  printout("[input.c] this simulation contains\n");
  printout("----------------------------------\n");
  for (int element = 0; element < nelements; element++)
  {
    printout("[input.c]   element %d (Z=%2d)\n", element, get_element(element));
    const int nions = get_nions(element);
    //includedions += nions;
    for (int ion = 0; ion < nions; ion++)
    {
      int photoiontransitions = 0;
      for (int level = 0; level < get_nlevels(element,ion); level++)
        photoiontransitions += get_nphixstargets(element,ion,level);
      printout("[input.c]      ion_stage %d with %d levels (%d ionising) and %d photoionisation transitions\n",
               get_ionstage(element,ion),get_nlevels(element,ion),get_ionisinglevels(element,ion),photoiontransitions);
      includedlevels += get_nlevels(element,ion);
      includedionisinglevels += get_ionisinglevels(element,ion);
      includedphotoiontransitions += photoiontransitions;
    }
  }
  printout("[input.c]   in total %d ions, %d levels (%d ionising), %d lines, %d photoionisation transitions\n",
           includedions,includedlevels,includedionisinglevels,nlines,includedphotoiontransitions);

  if ((bflist = (bflist_t *) malloc(includedionisinglevels*sizeof(bflist_t))) == NULL)
  {
    printout("[fatal] input: not enough memory to initialize bflist ... abort\n");
    abort();
  }

  write_bflist_file(includedionisinglevels);

  setup_phixs_list();

  ///set-up/gather information for nlte stuff

  total_nlte_levels = 0;
  n_super_levels = 0;

  if (NLTE_POPS_ON)
  {
    for (int element = 0; element < nelements; element++)
    {
      const int nions = get_nions(element);
      //includedions += nions;
      for (int ion = 0; ion < nions; ion++)
      {
        elements[element].ions[ion].first_nlte = total_nlte_levels;
        const int nlevels = get_nlevels(element,ion);
        int count = 0;
        if (nlevels > 1)
        {
          for (int level = 1; level < nlevels; level++)
          {
            if (is_nlte(element,ion,level))
            {
              count++;
              total_nlte_levels++;
            }
          }
        }

        if (count < (nlevels - 1))
        {
          // If there are more levels that the ground state + the number of NLTE levels then we need an extra
          // slot to store data for the "superlevel", which is a representation of all the other levels that
          // are not treated in detail.
          total_nlte_levels++;
          n_super_levels++;
        }

        elements[element].ions[ion].nlevels_nlte = count;

        printout("[input.c]  element %2d Z=%2d ion_stage %2d has %5d NLTE levels. Starting at %d\n",
                 element, get_element(element), get_ionstage(element,ion),
                 get_nlevels_nlte(element,ion), elements[element].ions[ion].first_nlte);
      }
    }
  }

  printout("[input.c] Total NLTE levels: %d, of which %d are superlevels\n", total_nlte_levels, n_super_levels);
}


static void read_1d_model(void)
/// Subroutine to read in a 1-D model.
{
  FILE *model_input;
  if ((model_input = fopen("model.txt", "r")) == NULL)
  {
    printout("Cannot open model.txt.\n");
    abort();
  }

  /* 1st read the number of data points in the table of input model. */
  int dum1;
  fscanf(model_input, "%d", &dum1);
  npts_model = dum1;
  if (npts_model > MMODELGRID)
  {
    printout("Too many points in input model. Abort.\n");
    abort();
  }
  /* Now read the time (in days) at which the model is specified. */
  float dum2;
  fscanf(model_input, "%g", &dum2);
  t_model = dum2 * DAY;

  /* Now read in the lines of the model. Each line has 5 entries: the
     cell number (integer) the velocity at outer boundary of cell (float),
     the mass density in the cell (float), the abundance of Ni56 by mass
     in the cell (float) and the total abundance of all Fe-grp elements
     in the cell (float). For now, the last number is recorded but never
     used. */

  for (int n = 0; n < npts_model; n++)
  {
    int dum1;
    float dum2,dum3,dum4,dum5,dum6,dum7,dum8;
    fscanf(model_input, "%d %g %g %g %g %g %g %g", &dum1, &dum2, &dum3, &dum4, &dum5, &dum6, &dum7, &dum8);
    vout_model[n] = dum2 * 1.e5;
    rho_model[n] = pow(10., dum3);
    ffegrp_model[n] = dum4;
    fni_model[n] = dum5;
    fco_model[n] = dum6;
    f52fe_model[n] = dum7;
    f48cr_model[n] = dum8;
  }

  fclose(model_input);

  /* Now compute mtot and mni56 for the model. */

  mtot = 0.0;
  mni56 = 0.0;
  mfe52 = 0.0;
  mcr48 = 0.0;
  mfeg = 0.0;

  /* do inner most cell 1st */
  double mass_in_shell = rho_model[0] * (pow(vout_model[0],3.)) * 4. * PI * pow(t_model,3.) / 3.;
  mtot += mass_in_shell;
  mni56 += mass_in_shell * fni_model[0];
  mfe52 += mass_in_shell * f52fe_model[0];
  mcr48 += mass_in_shell * f48cr_model[0];
  mfeg += mass_in_shell * ffegrp_model[0];

  /* Now do the rest. */
  for (int n = 1; n < npts_model; n++)
  {
    mass_in_shell = rho_model[n] * (pow(vout_model[n], 3) - pow(vout_model[n - 1], 3)) * 4 * PI * pow(t_model, 3) / 3;
    mtot += mass_in_shell;
    mni56 += mass_in_shell * fni_model[n];
    mfe52 += mass_in_shell * f52fe_model[n];
    mcr48 += mass_in_shell * f48cr_model[n];
    mfeg += mass_in_shell * ffegrp_model[n];
  }

  printout("Total mass: %g. 56Ni mass: %g. 52Fe mass: %g. 48Cr mass: %g. Fe-grp mass: %g.\n", mtot/MSUN, mni56/MSUN, mfe52/MSUN, mcr48/MSUN, mfeg/MSUN);

  vmax = vout_model[npts_model - 1];
  rmax = vmax * tmin;
  xmax = ymax = zmax = rmax;
  printout("rmax %g\n", rmax);
}


static void read_2d_model(void)
/// Subroutine to read in a 2-D model.
{
  FILE *model_input;
  if ((model_input = fopen("model.txt", "r")) == NULL)
  {
    printout("Cannot open model.txt.\n");
    abort();
  }

  /* 1st read the number of data points in the table of input model. */
  int dum1, dum12;
  fscanf(model_input, "%d %d", &dum1, &dum12);
  ncoord1_model = dum1; //r (cylindrical polar)
  ncoord2_model = dum12;//z (cylindrical polar)

  npts_model = dum1 * dum12;
  if (npts_model > MMODELGRID)
  {
    printout("Too many points in input model. Abort.\n");
    abort();
  }
  /* Now read the time (in days) at which the model is specified. */
  float dum2, dum3, dum4, dum5, dum6, dum7, dum8, dum9;
  fscanf(model_input, "%g", &dum2);
  t_model = dum2 * DAY;

  /* Now read in vmax (in cm/s) */
  fscanf(model_input, "%g", &dum2);
  vmax = dum2;
  dcoord1 = dum2 * t_model / ncoord1_model; //dr for input model
  dcoord2 = 2. * dum2 * t_model / ncoord2_model;//dz for input model


  /* Now read in the model. Each point in the model has two lines of input.
     First is an index for the cell then its r-mid point then its z-mid point
     then its total mass density.
     Second is the total FeG mass, initial 56Ni mass, initial 56Co mass */

  for (int n = 0; n < npts_model; n++)
  {
    fscanf(model_input, "%d %g %g %g", &dum1, &dum2,  &dum3, &dum4);
    fscanf(model_input, "%g %g %g %g %g",&dum5, &dum6, &dum7, &dum8, &dum9);
    rho_model[n] = dum4;
    ffegrp_model[n] = dum5;
    fni_model[n] = dum6;
    fco_model[n] = dum7;
    f52fe_model[n] = dum8;
    f48cr_model[n] = dum9;
  }

  fclose(model_input);

  /* Now compute mtot and mni56 for the model. */

  mtot = 0.0;
  mni56 = 0.0;
  mfe52 = 0.0;
  mcr48 = 0.0;
  mfeg = 0.0;

  int n1 = 0;
  /* Now do the rest. */
  for (int n = 0; n < npts_model; n++)
  {
    double mass_in_shell = rho_model[n] * ((2*n1) + 1) * PI * dcoord2 * pow(dcoord1,2.);
    mtot += mass_in_shell;
    mni56 += mass_in_shell * fni_model[n];
    mfe52 += mass_in_shell * f52fe_model[n];
    mcr48 += mass_in_shell * f48cr_model[n];
    mfeg += mass_in_shell * ffegrp_model[n];
    n1++;
    if (n1 == ncoord1_model)
    {
      n1 = 0;
    }
  }

  printout("Total mass: %g. Ni mass: %g.  52Fe mass: %g.  48Cr mass: %g. Fe-grp mass: %g.\n", mtot/MSUN, mni56/MSUN, mfe52/MSUN, mcr48/MSUN, mfeg/MSUN);

  rmax = vmax * tmin;
  xmax = ymax = zmax = rmax;
  printout("rmax %g\n", rmax);
}


static void read_3d_model(void)
/// Subroutine to read in a 3-D model.
{
  float dum2, dum3, dum4, dum5, dum6;
  float rho_model;
  double helper;

  FILE *model_input;
  if ((model_input = fopen("model.txt", "r")) == NULL)
  {
    printout("Cannot open model.txt.\n");
    abort();
  }

  /// 1st read the number of data points in the table of input model.
  /// This MUST be the same number as the maximum number of points used in the grid - if not, abort.
  int dum1;
  fscanf(model_input, "%d", &dum1);
  npts_model = dum1;
  if (npts_model > MMODELGRID)
  {
    printout("Too many points in input model. Abort. (%d > %d)\n",npts_model,MMODELGRID);
    abort();
  }
  if (npts_model != nxgrid * nygrid * nzgrid)
  {
    printout("3D model/grid mismatch. Abort. %d != %d\n",npts_model, nxgrid * nygrid * nzgrid);
    abort();
  }

  /// Now read the time (in days) at which the model is specified.
  fscanf(model_input, "%g", &dum2);
  t_model = dum2 * DAY;

  /// Now read in vmax for the model (in cm s^-1).
  fscanf(model_input, "%g", &dum2);
  vmax = dum2;

  /// Now read in the lines of the model.
  min_den=1.e99;

  /*mgi is the index to the model grid - empty cells are sent to MMODELGRID, otherwise each input cell is one modelgrid cell */

  int mgi = 0;
  for (int n = 0; n < npts_model; n++)
  {
    fscanf(model_input, "%d %g %g %g %g", &dum1, &dum3,  &dum4, &dum5, &rho_model);
    //printout("cell %d, posz %g, posy %g, posx %g, rho %g, rho_init %g\n",dum1,dum3,dum4,dum5,rho_model,rho_model* pow( (t_model/tmin), 3.));
    if (rho_model < 0)
    {
      printout("negative input density %g %d\n", rho_model, n);
      abort();
    }

    if (mgi > MMODELGRID-1)
    {
      printout("3D model wants more modelgrid cells than MMODELGRID. Abort.\n");
      abort();
    }

    if (rho_model > 0)
    {
      mg_associated_cells[mgi] = 1;
      cell[n].modelgridindex = mgi;
      helper = rho_model * pow( (t_model/tmin), 3.);
      //printout("mgi %d, helper %g\n",mgi,helper);
      set_rhoinit(mgi,helper);
      //printout("mgi %d, rho_init %g\n",mgi,get_rhoinit(mgi));
      set_rho(mgi,helper);

      if (min_den > rho_model)
      {
        min_den = rho_model;
      }
    }
    else
    {
      cell[n].modelgridindex = MMODELGRID;
    }

    fscanf(model_input, "%g %g %g %g %g", &dum2, &dum3,  &dum4, &dum5, &dum6);
    //printout("ffe %g, fni %g, fco %g, ffe52 %g, fcr48 %g\n",dum2,dum3,dum4,dum5,dum6);
    if (rho_model > 0)
    {
      set_ffe(mgi,dum2);
      set_f56ni(mgi,dum3);
      set_fco(mgi,dum4);
      set_f52fe(mgi,dum5);
      set_f48cr(mgi,dum6);
      allocate_compositiondata(mgi);
      allocate_cooling(mgi);
      //printout("mgi %d, control rho_init %g\n",mgi,get_rhoinit(mgi));
      mgi++;
    }
  }

  printout("min_den %g\n", min_den);
  printout("Effectively used model grid cells %d\n", mgi);

  /// Now, reduce the actual size of the modelgrid to the number of non-empty cells.
  /// Actually this doesn't reduce the memory consumption since all the related
  /// arrays are allocated statically at compile time with size MMODELGRID+1.
  /// However, it ensures that update_grid causes no idle tasks due to empty cells!
  npts_model = mgi;

  fclose(model_input);

  set_rhoinit(MMODELGRID,0.);
  set_rho(MMODELGRID,0.);
  set_nne(MMODELGRID,0.);
  set_ffe(MMODELGRID,0.);
  set_fco(MMODELGRID,0.);
  set_f56ni(MMODELGRID,0.);
  set_f52fe(MMODELGRID,0.);
  set_f48cr(MMODELGRID,0.);
  set_Te(MMODELGRID,MINTEMP);
  set_TJ(MMODELGRID,MINTEMP);
  set_TR(MMODELGRID,MINTEMP);
  allocate_compositiondata(MMODELGRID);

  /// Now compute mtot and mni56 for the model.
  /// Assumes cells are cubes here - all same volume.

  mtot = 0.0;
  mni56 = 0.0;
  mfe52 = 0.0;
  mcr48 = 0.0;
  mfeg = 0.0;

  //for (n = 0; n < npts_model; n++)
  for (int mgi = 0; mgi < npts_model; mgi++)
  {
    //mgi = cell[n].modelgridindex;
    double mass_in_shell = get_rhoinit(mgi);
    //printout("n %d, mgi %d, rho_init %g\n",n,mgi,mass_in_shell);
    mtot += mass_in_shell;
    mni56 += mass_in_shell * get_f56ni(mgi);
    mfe52 += mass_in_shell * get_f52fe(mgi);
    mcr48 += mass_in_shell * get_f48cr(mgi);
    mfeg += mass_in_shell * get_ffe(mgi);
  }

  double cellvolume = pow((2 * vmax * tmin),3.) / (nxgrid*nygrid*nzgrid);
  mtot = mtot * cellvolume;
  mni56 = mni56 * cellvolume;
  mfe52 = mfe52 * cellvolume;
  mcr48 = mcr48 * cellvolume;
  mfeg = mfeg  * cellvolume;

  printout("Total mass: %g. Ni mass: %g.  52Fe mass: %g.  48Cr mass: %g. Fe-group mass: %g.\n", mtot/MSUN, mni56/MSUN, mfe52/MSUN, mcr48/MSUN, mfeg/MSUN);

  rmax = vmax * tmin;
  xmax = ymax = zmax = rmax;
  printout("rmax %g\n", rmax);


  ///Opacity treatment moved to grid_init & update_grid
//    if (opacity_case == 0)
//      {
//        for (n = 0; n < npts_model; n++)
// 	 {
// 	   cell[n].kappa_grey = GREY_OP;
// 	 }
//      }
//    else if (opacity_case == 1)
//      {
//        for (n = 0; n < npts_model; n++)
// 	 {
// 	   cell[n].kappa_grey = ((0.9 * cell[n].f_fe) + 0.1) * GREY_OP / ((0.9 *  mfeg / mtot) + 0.1);
// 	 }
//      }
//    else if (opacity_case == 2)
//      {
//        for (n = 0; n < npts_model; n++)
// 	 {
// //           cell[n].kappa_grey = ((0.9 * cell[n].f_fe) + 0.1)/cell[n].rho_init * GREY_OP*rhotot / ((0.9 *  mfeg / mtot) + 0.1);
// //           cell[n].kappa_grey = ((0.9 * cell[n].f_fe) + 0.1)/cell[n].rho_init * GREY_OP*rhom / ((0.9 *  mfeg) + (0.1 * mtot));
//            cell[n].kappa_grey = ((0.9 * cell[n].f_fe) + 0.1)/cell[n].rho_init * GREY_OP*rho_sum / ((0.9 *  ni56_sum) + (0.1 * npts_model));
// 	 }
//      }
//    else
//      {
//        printout("Unknown opacity case. Abort.\n");
//        abort();
//      }
}


static void read_binding_energies(void)
{
  FILE *binding;
  if ((binding = fopen("binding_energies.txt", "r")) == NULL)
  {
    printout("Cannot open binding_energies.txt.\n");
    abort();
  }

  int dum1, dum2;
  fscanf(binding, "%d %d", &dum1, &dum2); //dimensions of the table
  if ((dum1 != M_NT_SHELLS) || (dum2 != MAX_Z_BINDING))
  {
    printout("Wrong size for the binding energy tables!\n");
    abort();
  }

  for (int index1 = 0; index1 < dum2; index1++)
  {
    float dum[10];
    fscanf(binding, "%g %g %g %g %g %g %g %g %g %g",
           &dum[0], &dum[1], &dum[2], &dum[3], &dum[4], &dum[5], &dum[6], &dum[7], &dum[8], &dum[9]);

    for (int index2 = 0; index2 < 10; index2++)
    {
      electron_binding[index1][index2] = dum[index2] * EV;
    }
  }

  fclose(binding);
}


void input(int rank)
/// To govern the input. For now hardwire everything.
{
  homogeneous_abundances = false;
  t_model = 0.0;

  /// Select grid type
  grid_type = GRID_UNIFORM;
  model_type = RHO_UNIFORM;

  maxion = MIONS;
  /// Set grid size
  //nxgrid = 4; //pow(MGRID,1./3.); //10;
  //nygrid = 4; //pow(MGRID,1./3.); //10;
  //nzgrid = 4; //pow(MGRID,1./3.); //10;
  nxgrid = 50;
  nygrid = 50;
  nzgrid = 50;
  printout("nxgrid %d\n",nxgrid);
  ngrid = nxgrid * nygrid * nzgrid; ///Moved to input.c
  if (ngrid > MGRID)
  {
    printout("[fatal] input: Error: too many grid cells. Abort. %d>%d",ngrid,MGRID);
    abort();
  }

  /// Set number of packets, outer and middle iterations
  npkts = MPKTS;
  n_out_it = 10;
  n_middle_it = 1;
/*  #ifdef FORCE_LTE
    n_titer = 1;
  #else
    n_titer = 6;
  #endif*/
  n_titer = 1;
  initial_iteration = false;

  printout("[info] input: do n_titer %d iterations per timestep\n",n_titer);
  if (n_titer > 1)
  {
    #ifndef DO_TITER
      printout("[fatal] input: n_titer > 1, but DO_TITER not defined ... abort\n");
      abort();
    #endif
  }
  else if (n_titer == 1)
  {
    #ifdef DO_TITER
      printout("[warning] input: n_titer = 1 but DO_TITER defined, remove DO_TITER to save memory\n");
    #endif
  }
  else
  {
    printout("[fatal] input: no valid value for n_titer selected\n");
    abort();
  }

  nu_min_r = 1e14;   /// lower frequency boundary for UVOIR spectra and BB sampling
  nu_max_r = 5e15;   /// upper frequency boundary for UVOIR spectra and BB sampling
  #ifdef DO_EXSPEC
    /// Spectra settings
    nnubins = MNUBINS; //1000;  /// frequency bins for spectrum
    /// lower and upper frequency boundaries for make_spectrum_res: gamma spectra?
    nu_min = 0.05 * MEV / H;
    nu_max = 4 * MEV / H;
  #endif

  /// Lightcurve setting
  do_r_lc = false;    /// default to no lc = gamma-ray spectrum
  do_rlc_est = 0; /// ^^

  nfake_gam = 1; ///# of fake gamma ray lines for syn

  /// Read in parameters from input.txt
  ///======================================================
  read_parameterfile(rank);
  ntbins = ntstep;   ///time bins for spectrum equal #(timesteps)
  ntlcbins = ntstep; ///time bins for light curve #(timesteps)

  /// Read in atomic data
  ///======================================================
  read_atomicdata();

  //#ifndef DO_EXSPEC
    /// Read in input model
    ///======================================================
    if (model_type == RHO_UNIFORM)
    {
      mtot = 1.39 * MSUN;
      mni56 = 0.625 * MSUN;
      vmax = 1.e9;
      rmax = vmax * tmin;
  //  rhotot = 3 * mtot / 4 / PI / rmax /rmax /rmax; //MK
      xmax = ymax = zmax = rmax;
    }
    else if (model_type == RHO_1D_READ)
    {
      printout("Read 1D model!\n");
      read_1d_model();
    }
    else if (model_type == RHO_2D_READ)
    {
      printout("Read 2D model!\n");
      read_2d_model();
    }
    else if (model_type == RHO_3D_READ)
    {
      printout("Read 3D model!\n");
      read_3d_model();
    }
    else
    {
      printout("Unknown model. Abort.\n");
      abort();
    }
  //#endif


  /// Read in data for gamma ray lines.
  ///======================================================
  FILE *co_lines;
  if ((co_lines = fopen("co_lines.txt", "r")) == NULL)
  {
    printout("Cannot open co_lines.txt.\n");
    abort();
  }
  int dum1 = 0;
  fscanf(co_lines, "%d", &dum1);
  cobalt_spec.nlines = dum1;

  ECOBALT = 0.0;
  for (int n = 0; n < dum1; n++)
  {
    float dum2, dum3;
    fscanf(co_lines, "%g %g", &dum2, &dum3);
    cobalt_spec.energy[n] = dum2 * MEV;
    cobalt_spec.probability[n] = dum3;
    ECOBALT += dum2 * MEV * dum3;
  }
  /// Average energy per gamma line of Co56 decay and positron annihilation
  ECOBALT_GAMMA = ECOBALT;
  /// For total deposited energy we need to add the kinetic energy per emitted positron
  ECOBALT += 0.63*MEV * 0.19;
  //ECOBALT = ECOBALT/5; /// DEBUGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG

  fclose(co_lines);

  FILE *ni_lines;
  if ((ni_lines = fopen("ni_lines.txt", "r")) == NULL)
  {
    printout("Cannot open ni_lines.txt.\n");
    abort();
  }
  fscanf(ni_lines, "%d", &dum1);
  nickel_spec.nlines = dum1;

  ENICKEL = 0.0;
  for (int n = 0; n < dum1; n++)
  {
    float dum2, dum3;
    fscanf(ni_lines, "%g %g", &dum2, &dum3);
    nickel_spec.energy[n] = dum2 * MEV;
    nickel_spec.probability[n] = dum3;
    ENICKEL += dum2 * MEV * dum3;
  }
  //ENICKEL = ENICKEL/5; /// DEBUGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG

  fclose(ni_lines);

  FILE *v48_lines;
  if ((v48_lines = fopen("v48_lines.txt", "r")) == NULL)
  {
    printout("Cannot open v48_lines.txt.\n");
    abort();
  }
  fscanf(v48_lines, "%d", &dum1);
  v48_spec.nlines = dum1;

  E48V = 0.0;
  for (int n = 0; n < dum1; n++)
  {
    float dum2, dum3;
    fscanf(v48_lines, "%g %g", &dum2, &dum3);
    v48_spec.energy[n] = dum2 * MEV;
    v48_spec.probability[n] = dum3;
    E48V += dum2 * MEV * dum3;
  }

  fclose(v48_lines);

  FILE *cr48_lines;
  if ((cr48_lines = fopen("cr48_lines.txt", "r")) == NULL)
  {
    printout("Cannot open cr48_lines.txt.\n");
    abort();
  }
  fscanf(cr48_lines, "%d", &dum1);
  cr48_spec.nlines = dum1;

  E48CR = 0.0;
  for (int n = 0; n < dum1; n++)
  {
    float dum2, dum3;
    fscanf(cr48_lines, "%g %g", &dum2, &dum3);
    cr48_spec.energy[n] = dum2 * MEV;
    cr48_spec.probability[n] = dum3;
    E48CR += dum2 * MEV * dum3;
  }

  fclose(cr48_lines);

  /// With the lines all ready in, now make a list of them in energy order.
  get_gam_ll();

  /// Now that the list exists use it to find values for spectral synthesis
  /// stuff.
  const int lindex_max = get_nul(nusyn_max);
  const int lindex_min = get_nul(nusyn_min);
  printout("lindex_max %d, lindex_min %d\n", lindex_max, lindex_min);

  emiss_offset = lindex_min;
  emiss_max = lindex_max - lindex_min + 1;
  printout("emiss_max using %d of a possible %d\n", emiss_max, EMISS_MAX);

  if (emiss_max > EMISS_MAX)
  {
    printout("Too many points needed for emissivities. Use smaller frequency range or increase EMISS_MAX. Abort.\n");
    abort();
  }

  if (NT_ON)
    read_binding_energies();

  #ifdef DO_EXSPEC
    /// Check if enough  memory for spectra has been assigned
    /// and allocate memory for the emission statistics
    if (nnubins > MNUBINS)
    {
      printout("Too many frequency bins in spectrum - reducing.\n");
      nnubins = MNUBINS;
    }
    if (ntbins > MTBINS)
    {
      printout("Too many time bins in spectrum - reducing.\n");
      ntbins = MTBINS;
    }
    for (int n = 0; n < ntbins; n++)
    {
      for (int m = 0; m < nnubins; m++)
      {
        if ((spectra[n].stat[m].absorption = (double *) malloc((nelements*maxion)*sizeof(double))) == NULL)
        {
          printout("[fatal] input: not enough memory to spectrum structure ... abort\n");
          abort();
        }
        if ((spectra[n].stat[m].emission = (double *) malloc((2*nelements*maxion+1)*sizeof(double))) == NULL)
        {
          printout("[fatal] input: not enough memory for spectrum structure ... abort\n");
          abort();
        }
        if ((spectra[n].stat[m].trueemission = (double *) malloc((2*nelements*maxion+1)*sizeof(double))) == NULL)
        {
          printout("[fatal] input: not enough memory for spectrum structure ... abort\n");
          abort();
        }
        /*
        if (do_emission_res == 1)
        {
          for (nn = 0; nn < MABINS; nn++)
          {
            if ((spectra_res[n][nn].emission[m].count = malloc((2*nelements*maxion+1)*sizeof(int))) == NULL)
            {
              printout("[fatal] input: not enough memory to initialise spectra_res structure ... abort\n");
              abort();
            }
          }
        }
        */
      }
    }
  #endif
}


void read_parameterfile(int rank)
/// Subroutine to read in input parameters from input.txt.
{
  unsigned long int pre_zseed;

  FILE *input_file = fopen("input.txt", "r");
  if (input_file == NULL)
  {
    printout("Cannot open input.txt.\n");
    abort();
  }

  int dum1;
  fscanf(input_file, "%d", &dum1);
  if (dum1 > 0)
  {
    pre_zseed = dum1; // random number seed
  }
  else
  {
    pre_zseed = time(NULL);
    printout("[debug] randomly-generated random number seed is %d\n", pre_zseed);
  }

  #ifdef _OPENMP
    #pragma omp parallel
    {
  #endif
      /// For MPI parallelisation, the random seed is changed based on the rank of the process
      /// For OpenMP parallelisation rng is a threadprivate variable and the seed changed according
      /// to the thread-ID tid.
      unsigned long int zseed = pre_zseed + (13 * rank) + (17 * tid); /* rnum generator seed */
      printout("rank %d: thread %d has zseed %d\n", rank, tid, zseed);
      /// start by setting up the randon number generator
      rng = gsl_rng_alloc(gsl_rng_ran3);
      gsl_rng_set(rng, zseed);
      /// call it a few times to get it in motion.
      for (int n = 0; n < 100; n++)
      {
        //double x = gsl_rng_uniform(rng);
        //printout("zrand %g\n", x);
        gsl_rng_uniform(rng);
      }
      printout("rng is a '%s' generator\n", gsl_rng_name(rng));
  #ifdef _OPENMP
    }
  #endif

  fscanf(input_file, "%d", &ntstep); // number of time steps

  fscanf(input_file, "%d %d", &itstep, &ftstep); // number of start and end time step

  float dum2, dum3;
  fscanf(input_file, "%g %g", &dum2, &dum3); // start and end times
  tmin = dum2 * DAY;
  tmax = dum3 * DAY;
  //tlimit = dum4 * DAY;

  fscanf(input_file, "%g %g", &dum2, &dum3);
  nusyn_min = dum2 * MEV / H; // lowest frequency to synthesise
  nusyn_max = dum3 * MEV / H; // highest frequency to synthesise

  fscanf(input_file, "%d", &nsyn_time); // number of times for synthesis

  fscanf(input_file, "%g %g", &dum2, &dum3); // start and end times for synthesis
  for (int i = 0; i < nsyn_time; i++)
  {
    time_syn[i] = exp(log(dum2) + (dum3 * i)) * DAY;
  }

  fscanf(input_file, "%d", &dum1); // model type
  if (dum1 == 1)
  {
    model_type = RHO_1D_READ;
  }
  else if (dum1 == 2)
  {
    model_type = RHO_2D_READ;
  }
  else if (dum1 == 3)
  {
    model_type = RHO_3D_READ;
  }

  fscanf(input_file, "%d", &dum1); ///compute the r-light curve?
  if (dum1 == 1) /// lc no estimators
  {
    do_r_lc = true;
    do_rlc_est = 0;
  }
  else if (dum1 == 2)/// lc case with thin cells
  {
    do_r_lc = true;
    do_rlc_est = 1;
  }
  else if (dum1 == 3)/// lc case with thick cells
  {
    do_r_lc = true;
    do_rlc_est = 2;
  }
  else if (dum1 == 4) /// gamma-ray heating case
  {
    do_r_lc = true;
    do_rlc_est = 3;
  }
  else if (dum1 != 0)
  {
    printout("Unknown rlc mode. Abort.\n");
    abort();
  }

  fscanf(input_file, "%d", &n_out_it); // number of iterations

  fscanf(input_file, "%g", &dum2); // change speed of light?
  CLIGHT_PROP = dum2 * CLIGHT;

  fscanf(input_file, "%lg", &gamma_grey); ///use grey opacity for gammas?

  float syn_dir_in[3];
  fscanf(input_file, "%g %g %g", &syn_dir_in[0], &syn_dir_in[1], &syn_dir_in[2]); // components of syn_dir

  const double rr = (syn_dir_in[0] * syn_dir_in[0]) + (syn_dir_in[1] * syn_dir_in[1]) + (syn_dir_in[2] * syn_dir_in[2]);
  /// ensure that this vector is normalised.
  if (rr > 1.e-6)
  {
    syn_dir[0] = syn_dir_in[0] / sqrt(rr);
    syn_dir[1] = syn_dir_in[1] / sqrt(rr);
    syn_dir[2] = syn_dir_in[2] / sqrt(rr);
  }
  else
  {
    const double z1 = 1. - (2 * gsl_rng_uniform(rng));
    const double z2 = gsl_rng_uniform(rng) * 2.0 * PI;
    syn_dir[2] = z1;
    syn_dir[0] = sqrt( (1. - (z1 * z1))) * cos(z2);
    syn_dir[1] = sqrt( (1. - (z1 * z1))) * sin(z2);
  }

  fscanf(input_file, "%d", &opacity_case); // opacity choice

  ///MK: begin
  fscanf(input_file, "%lg", &rho_crit_para); ///free parameter for calculation of rho_crit
  printout("input: rho_crit_para %g\n", rho_crit_para);
  ///the calculation of rho_crit itself depends on the time, therfore it happens in grid_init and update_grid

  fscanf(input_file, "%d", &debug_packet); /// activate debug output for packet
                                           /// select a negative value to deactivate
  ///MK: end

  /// Do we start a new simulation or, continue another one?
  int continue_flag;
  fscanf(input_file, "%d", &continue_flag);
  simulation_continued_from_saved = (continue_flag == 1);
  if (simulation_continued_from_saved)
    printout("input: resuming simulation from saved point\n");
  else
    printout("input: starting a new simulation\n");

  /// Wavelength (in Angstroms) at which the parameterisation of the radiation field
  /// switches from the nebular approximation to LTE.
  fscanf(input_file, "%g", &dum2);  // free parameter for calculation of rho_crit
  nu_rfcut = CLIGHT / (dum2 * 1e-8);
  printout("input: nu_rfcut %g\n",nu_rfcut);

  /// Sets the number of initial LTE timesteps for NLTE runs
  fscanf(input_file, "%d", &n_lte_timesteps);
  #ifdef FORCE_LTE
    printout("input: this is a pure LTE run\n");
  #else
    printout("input: this is a NLTE run\n");
    printout("input: do the first %d timesteps in LTE\n", n_lte_timesteps);
  #endif

  if (NT_ON)
  {
    if (NT_SOLVE_SPENCERFANO)
      printout("input: Non-thermal ionisation with a Spencer-Fano solution is switched on for this run.\n");
    else
      printout("input: Non-thermal ionisation with the work function approximation is switched on for this run.\n");
    #ifdef FORCE_LTE
      printout("input: Non-thermal ionisation requires the code to run in non-LTE mode. Remove macro FORCE_LTE and recompile!\n");
      abort();
    #endif
  }
  else
    printout("input: No non-thermal ionisation is used in this run.\n");

  #if (NO_LUT_PHOTOION)
    printout("Corrphotoioncoeff is calculated from the radiation field at each timestep in each modelgrid cell (no LUT).\n");
  #else
    printout("Corrphotoioncoeff is calculated from LTE lookup tables (ratecoeff.dat) and corrphotoionrenorm estimator.\n");
  #endif

  #if (NO_LUT_BFHEATING)
    printout("bfheating coefficients are calculated from the radiation field at each timestep in each modelgrid cell (no LUT).\n");
  #else
    printout("bfheating coefficients are calculated from LTE lookup tables (ratecoeff.dat) and bfheatingestimator.\n");
  #endif

  if (MULTIBIN_RADFIELD_MODEL_ON)
    printout("The multibin radiation field estimators are being used instead of the whole-spectrum fit from timestep %d onwards.\n", FIRST_NLTE_RADFIELD_TIMESTEP);
  else
    printout("The radiation field model is a whole-spectrum fit to a single diluted blackbody.\n");

  /// Set up initial grey approximation?
  fscanf(input_file, "%lg %d", &cell_is_optically_thick, &n_grey_timesteps);
  printout("input: cells with Thomson optical depth > %g are treated in grey approximation for the first %d timesteps\n",cell_is_optically_thick,n_grey_timesteps);

  /// Limit the number of bf-continua
  fscanf(input_file, "%d", &max_bf_continua);
  if (max_bf_continua == -1)
  {
    printout("input: use all bf-continua\n");
    max_bf_continua = 1e6;
  }
  else
  {
    printout("input: use only %d bf-continua per ion\n", max_bf_continua);
  }

  /// The following parameters affect the DO_EXSPEC mode only /////////////////
  /// Read number of MPI tasks for exspec
  fscanf(input_file, "%d", &dum1);
  #ifdef DO_EXSPEC
    nprocs_exspec = dum1;
    printout("input: DO_EXSPEC ... extract spectra for %d MPI tasks\n",nprocs_exspec);
    printout("input: DO_EXSPEC ... and %d packets per task\n",npkts);
  #endif

  /// Extract line-of-sight dependent information of last emission for spectrum_res
  ///   if 1, yes
  ///   if 0, no
  /// Only affects runs with DO_EXSPEC. But then it needs huge amounts of memory!!!
  /// Thus be aware to reduce MNUBINS for this mode of operation!
  fscanf(input_file, "%d", &dum1);
  #ifdef DO_EXSPEC
    do_emission_res = dum1;
    if (do_emission_res == 1) printout("input: DO_EXSPEC ... extract LOS dependent emission information\n");
  #endif

  /// To reduce the work imbalance between different MPI tasks I introduced a diffusion
  /// for kpkts, since it turned out that this work imbalance was largely dominated
  /// by continuous collisional interactions. By introducing a diffusion time for kpkts
  /// this loop is broken. The following two parameters control this approximation.
  /// Parameter one (a float) gives the relative fraction of a time step which individual
  /// kpkts live. Parameter two (an int) gives the number of time steps for which we
  /// want to use this approximation
  fscanf(input_file, "%g %d", &kpktdiffusion_timescale, &n_kpktdiffusion_timesteps);
  printout("input: kpkts diffuse %g of a time step's length for the first %d time steps\n", kpktdiffusion_timescale, n_kpktdiffusion_timesteps);

  fclose(input_file);
}



/*
int compare_linelistentry(const void *p1, const void *p2)
{
  linelist_entry *a1, *a2;
  a1 = (linelist_entry *)(p1);
  a2 = (linelist_entry *)(p2);
  //printf("%d %d %d %d %g\n",a1->elementindex,a1->ionindex,a1->lowerlevelindex,a1->upperlevelindex,a1->nu);
  //printf("%d %d %d %d %g\n",a2->elementindex,a2->ionindex,a2->lowerlevelindex,a2->upperlevelindex,a2->nu);
  //printf("%g\n",a2->nu - a1->nu);
  if (a1->nu - a2->nu < 0)
    return 1;
  else if (a1->nu - a2->nu > 0)
    return -1;
  else
    return 0;
}
*/


void update_parameterfile(int nts)
/// Subroutine to read in input parameters from input.txt.
{
  FILE *input_file;
  if ((input_file = fopen("input.txt", "r+")) == NULL)
  {
    printout("Cannot open input.txt.\n");
    abort();
  }
  //setvbuf(input_file, NULL, _IOLBF, 0);

  int dum1;
  fscanf(input_file, "%d\n", &dum1); /// Random number seed
  fscanf(input_file, "%d\n", &dum1); /// Number of time steps

  /// Number of start and end time step, update start time step
  fseek(input_file,0,SEEK_CUR);
  fprintf(input_file, "%3.3d %3.3d\n", nts, ftstep);
  fseek(input_file,0,SEEK_CUR);

  float dum2, dum3, dum4;
  fscanf(input_file, "%g %g\n", &dum2, &dum3);  ///start and end times
  fscanf(input_file, "%g %g\n", &dum2, &dum3);  ///lowest and highest frequency to synthesise in MeV
  fscanf(input_file, "%d\n", &dum1);            ///number of times for synthesis
  fscanf(input_file, "%g %g\n", &dum2, &dum3);  ///start and end times for synthesis
  fscanf(input_file, "%d\n", &dum1);            ///model type
  fscanf(input_file, "%d\n", &dum1);            ///compute the r-light curve?
  fscanf(input_file, "%d\n", &dum1);            ///number of iterations
  fscanf(input_file, "%g\n", &dum2);            ///change speed of light?
  fscanf(input_file, "%g\n", &dum2);            ///use grey opacity for gammas?
  fscanf(input_file, "%g %g %g\n", &dum2, &dum3, &dum4); ///components of syn_dir
  fscanf(input_file, "%d\n", &dum1);            ///opacity choice
  fscanf(input_file, "%g\n", &dum2);            ///free parameter for calculation of rho_crit
  fscanf(input_file, "%d\n", &dum1);            /// activate debug output for packet

  /// Flag continuation parameter as active
  fseek(input_file,0,SEEK_CUR);
  fprintf(input_file, "%d\n", 1); /// Force continuation
  //  fseek(input_file,0,SEEK_CUR);
  //fscanf(input_file, "%d", &dum1);  /// Do we start a new simulation or, continue another one?

  fclose(input_file);
}



