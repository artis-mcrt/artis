#include "sn3d.h"
#include "atomic.h"
#include "grid_init.h"
#include "rpkt.h"
#include "vectors.h"


///****************************************************************************
/// Routine for doing a uniform density grid.
/*int uniform_density_setup ()
{
  int n;
  double vec_len();
  double dcen[3];
  double fni(CELL *grid_ptr);
  void allocate_compositiondata(int cellnumber);
  double opcase2_normal, opcase3_sum, check1, check2;
  int empty_cells;

  //Calculate the critical opacity at which opacity_case 3 switches from a
  //regime proportional to the density to a regime independent of the density
  //This is done by solving for tau_sobolev == 1
  //tau_sobolev = PI*QE*QE/(ME*C) * rho_crit_para * rho/MNI56 * 3000e-8 * time_step[m].mid;
  rho_crit = ME*CLIGHT*MNI56 / (PI*QE*QE * rho_crit_para * 3000e-8 * tmin);
  printout("grid_init: rho_crit = %g\n", rho_crit);

  rho_sum = 0.0;
  fe_sum = 0.0;
  empty_cells = 0;
  check1 = check2 = 0.0;
  opcase3_sum = 0;

  //MK: first loop to initialize basic variables
  for (n=0;n < ngrid; n++)
  {
    dcen[0] = cell[n].pos_init[0] + (0.5*wid_init);
    dcen[1] = cell[n].pos_init[1] + (0.5*wid_init);
    dcen[2] = cell[n].pos_init[2] + (0.5*wid_init);


    if (vec_len(dcen) < rmax)
    {
      cell[n].rho_init = cell[n].rho = 3 * mtot / 4 / PI / rmax /rmax /rmax;
      cell[n].f_ni = cell[n].f_fe = fni(&cell[n]);
      cell[n].f_co = 0.;
      allocate_compositiondata(n);
    //printout("[debug] uniform_density_setup: cell[%d].rho_init: %g\n",n,cell[n].rho_init);
    }
     //MK: added missing else branch for correct initialization
    else
    {
      cell[n].rho_init = cell[n].rho = 0.0;
      cell[n].kappa_grey = 0.0;
      empty_cells += 1;
    }
    rho_sum += cell[n].rho_init;
    fe_sum += fni(&cell[n]);  //MK: use for this model fni as values for iron group mass fractions
  }


  if (opacity_case == 3)
  {
    for (n = 0; n < ngrid; n++)
    {
      if (cell[n].rho_init > 0)
      {
        if (cell[n].rho_init > rho_crit)
        {
          cell[n].kappa_grey = (0.9 * fni(&cell[n]) + 0.1) * rho_crit/cell[n].rho_init;
        }
        else
        {
          cell[n].kappa_grey = (0.9 * fni(&cell[n]) + 0.1);
        }
      }
      else if (cell[n].rho_init == 0)
      {
        cell[n].kappa_grey = 0;
      }
      else if (cell[n].rho_init < 0)
      {
        printout("Error: negative density. Abort.\n");
        exit(0);
      }
      opcase3_sum += cell[n].kappa_grey*cell[n].rho_init;
    }
  }


  //MK: second loop to set up opacities which need integrated basic quantities for normalisation
  for (n=0;n < ngrid; n++)
  {
    if (cell[n].rho_init > 0)
    {
      if (opacity_case == 0)
      {
        cell[n].kappa_grey = GREY_OP;
      }
      else if (opacity_case == 1)
      {
        cell[n].kappa_grey = ((0.9 * fni(&cell[n])) + 0.1) * GREY_OP / ((0.9 *  mni56 / mtot) + 0.1);
      }
      else if (opacity_case == 2)
      {
        opcase2_normal = GREY_OP*rho_sum / ((0.9 *  fe_sum) + (0.1 * (ngrid - empty_cells)));
        cell[n].kappa_grey = opcase2_normal/cell[n].rho_init * ((0.9 * fni(&cell[n])) + 0.1);
      }
      else if (opacity_case == 3)
      {
        opcase3_normal = GREY_OP*rho_sum / opcase3_sum;
        cell[n].kappa_grey = cell[n].kappa_grey * opcase3_normal;
      }
      else if (opacity_case == 4)
      {
        cell[n].kappa_grey = SIGMA_T;
      }
      else
      {
        printout("Unknown opacity case. Abort.\n");
        exit(0);
      }
    }
    else if (cell[n].rho_init == 0)
    {
      cell[n].kappa_grey = 0.0;
    }
    else if (cell[n].rho_init < 0)
    {
      printout("Error: negative density. Abort.\n");
      exit(0);
    }

    check1 = check1 + (cell[n].kappa_grey  * cell[n].rho_init);
    check2 = check2 + cell[n].rho_init;
  }


  printout("Initial densities assigned uniformly.\n");
  printout("Grey normalisation check: %g\n", check1/check2);

  return 0;
}*/



///****************************************************************************
/// Routine for doing a density grid read from a 1-D model.
static int density_1d_read(void)
{
  //int renorm[MMODELGRID];
  //double den_norm[MMODELGRID];

  //Calculate the critical opacity at which opacity_case 3 switches from a
  //regime proportional to the density to a regime independent of the density
  //This is done by solving for tau_sobolev == 1
  //tau_sobolev = PI*QE*QE/(ME*C) * rho_crit_para * rho/MNI56 * 3000e-8 * time_step[m].mid;
  rho_crit = ME*CLIGHT*MNI56 / (PI*QE*QE * rho_crit_para * 3000e-8 * tmin);
  printout("grid_init: rho_crit = %g\n", rho_crit);

  /*for (n=0;n<MMODELGRID;n++)
  {
    renorm[n]=0;
    den_norm[n]=0.0;
  }*/

  rho_sum = 0.0;
  fe_sum = 0.0;
  double opcase3_sum = 0.0;
  int empty_cells = 0;

  /*
    for (n=0;n < ngrid; n++)
  {
    dcen[0] = cell[n].pos_init[0] + (0.5*wid_init);
    dcen[1] = cell[n].pos_init[1] + (0.5*wid_init);
    dcen[2] = cell[n].pos_init[2] + (0.5*wid_init);

    radial_pos = vec_len(dcen);
    if (radial_pos < rmax)
    {
      mkeep =0;
      cell[n].modelgridindex = 0;

      helper = rho_model[0] * pow( (t_model/tmin), 3.);
      set_rhoinit(0,helper);
      set_rho(0,helper);
      set_ffe(0,ffegrp_model[0]);
      set_fni(0,fni_model[0]);
      set_fco(0,fco_model[0]);
      //allocate_compositiondata(0);
      for (element = 0; element < nelements; element++)
      {
        ///now set the abundances (by mass) of included elements, i.e.
        ///read out the abundances specified in the atomic data file
        anumber = get_element(element);
        abundance = abund_model[0][anumber-1];
        //cell[n].composition[element].abundance = abundance;
        modelgrid[0].composition[element].abundance = abundance;
        //if (anumber == 8)
        //  cell[n].composition[element].abundance = 1-cell[n].f_ni;
        //if (anumber == 20)
        //  cell[n].composition[element].abundance = 0.1039; ///force Ca to have higher abundance (equal to O abundance) in a pure Ca run

        //printout("element %d has abundance %g in cell %d\n",element,cell[n].composition[element].abundance,n);

        if (anumber == 28)
          set_fnistable(0,abundance - get_fni(n));
        if (anumber == 27)
          set_fcostable(0,abundance - get_fco(n));
        if (anumber == 26)
          set_ffeinit(0,abundance);
      }

      for (m = 0; m < (npts_model-1); m++)
      {
        if (radial_pos > (vout_model[m] * tmin))
        {
          mkeep = m+1;
          cell[n].modelgridindex = m+1;

          helper = rho_model[m+1] * pow( (t_model/tmin), 3.);
          set_rhoinit(m+1,helper);
          set_rho(m+1,helper);
          set_ffe(m+1,ffegrp_model[m+1]);
          set_fni(m+1,fni_model[m+1]);
          set_fco(m+1,fco_model[m+1]);
          //allocate_compositiondata(m+1);
    for (element = 0; element < nelements; element++)
          {
            ///now set the abundances (by mass) of included elements, i.e.
            ///read out the abundances specified in the atomic data file
            anumber = get_element(element);
            abundance = abund_model[m+1][anumber-1];
            modelgrid[m+1].composition[element].abundance = abundance;
            //if (anumber == 8)
            //  cell[n].composition[element].abundance = 1-cell[n].f_ni;
            //printout("element %d has abundance %g in cell %d\n",element,cell[n].composition[element].abundance,n);
            //if (anumber == 20)
            //  cell[n].composition[element].abundance = 0.1039; ///force Ca to have higher abundance (equal to O abundance) in a pure Ca run

            if (anumber == 28)
              set_fnistable(m+1,abundance - get_fni(n));
            if (anumber == 27)
              set_fcostable(m+1,abundance - get_fco(n));
            if (anumber == 26)
              set_ffeinit(m+1,abundance);
          }
        }
      }
      renorm[mkeep] += 1;
    }
    else
    {
      cell[n].modelgridindex = MMODELGRID;
      /// This initialisation is not needed if the cell MMODELGRID (which is a representation
      /// of empty cells) is once initialised to zero values
      set_rhoinit(MMODELGRID,0.);
      set_rho(MMODELGRID,0.);
      set_nne(MMODELGRID,0.);
      set_ffe(MMODELGRID,0.);
      set_fco(MMODELGRID,0.);
      set_fni(MMODELGRID,0.);
      set_Te(MMODELGRID,MINTEMP);
      set_TJ(MMODELGRID,MINTEMP);
      set_TR(MMODELGRID,MINTEMP);
      allocate_compositiondata(MMODELGRID);
      empty_cells += 1;
    }

    mgi = cell[n].modelgridindex;
    if (get_rhoinit(mgi) < 0)
    {
      printout("Error: negative density. Abort.\n");
      exit(0);
    }
    rho_sum += get_rhoinit(mgi);
    fe_sum += get_ffe(mgi);
  }
  */


  for (int n = 0; n < ngrid; n++)
  {
    double dcen[3];
    dcen[0] = cell[n].pos_init[0] + (0.5*wid_init);
    dcen[1] = cell[n].pos_init[1] + (0.5*wid_init);
    dcen[2] = cell[n].pos_init[2] + (0.5*wid_init);

    double radial_pos = vec_len(dcen);
    if (radial_pos < rmax)
    {
      //int mkeep = 0;
      cell[n].modelgridindex = 0;

      for (int m = 0; m < (npts_model-1); m++)
      {
        if (radial_pos > (vout_model[m] * tmin))
        {
          //mkeep = m + 1;
          cell[n].modelgridindex = m + 1;
        }
      }
      modelgrid[cell[n].modelgridindex].initial_radial_pos += radial_pos;
      //renorm[mkeep] += 1;
    }
    else
    {
      cell[n].modelgridindex = MMODELGRID;
    }
  }


  /// Determine the number of simulation cells associated with the model cells
  for (int i = 0; i < npts_model; i++)
  {
    int count = 0;
    for (int ii = 0; ii < ngrid; ii++)
    {
      if (cell[ii].modelgridindex == i) count++;
    }
    modelgrid[i].associated_cells = count;
  }


//  for (n=0;n < ngrid; n++)
//  {
//    mgi = cell[n].modelgridindex;
  for (int mgi = 0; mgi < npts_model; mgi++)
  {
    if (modelgrid[mgi].associated_cells > 0)
    {
      double helper = rho_model[mgi] * pow( (t_model/tmin), 3.);
      set_rhoinit(mgi,helper);
      set_rho(mgi,helper);
      set_ffe(mgi,ffegrp_model[mgi]);
      set_fni(mgi,fni_model[mgi]);
      set_fco(mgi,fco_model[mgi]);
      set_f52fe(mgi,f52fe_model[mgi]);
      set_f48cr(mgi,f48cr_model[mgi]);
      allocate_compositiondata(mgi);
      allocate_cooling(mgi);
      for (int element = 0; element < nelements; element++)
      {
        ///now set the abundances (by mass) of included elements, i.e.
        ///read out the abundances specified in the atomic data file
        int anumber = get_element(element);
        float abundance = abund_model[mgi][anumber-1];
        //cell[n].composition[element].abundance = abundance;
        modelgrid[mgi].composition[element].abundance = abundance;
        //if (anumber == 8)
        //  cell[n].composition[element].abundance = 1-cell[n].f_ni;
        //if (anumber == 20)
        //  cell[n].composition[element].abundance = 0.1039; ///force Ca to have higher abundance (equal to O abundance) in a pure Ca run

        //printout("element %d has abundance %g in cell %d\n",element,cell[n].composition[element].abundance,n);

        if (anumber == 28)
        {
          set_fnistable(mgi,abundance - get_fni(mgi));
          //printout("mgi %d, ni_abund %g, fni %g, fnistable %g\n",mgi,abundance,get_fni(mgi),abundance - get_fni(mgi));
        }
        else if (anumber == 27)
          set_fcostable(mgi,abundance - get_fco(mgi));
        else if (anumber == 26)
          set_ffestable(mgi,abundance - get_f52fe(mgi));
        else if (anumber == 25)
          set_fmnstable(mgi,abundance);
        else if (anumber == 24)
          set_fcrstable(mgi,abundance - get_f48cr(mgi));
        else if (anumber == 23)
          set_fvstable(mgi,abundance);
        else if (anumber == 22)
          set_ftistable(mgi,abundance);
      }

      if (get_rhoinit(mgi) < 0)
      {
        printout("Error: negative density. Abort.\n");
        exit(0);
      }
    }
  }

  /// This is the placeholder for empty cells. Temperatures must be positive
  /// as long as ff opacities are calculated.
  set_rhoinit(MMODELGRID,0.);
  set_rho(MMODELGRID,0.);
  set_nne(MMODELGRID,0.);
  set_ffe(MMODELGRID,0.);
  set_fco(MMODELGRID,0.);
  set_fni(MMODELGRID,0.);
  set_f48cr(MMODELGRID,0.);
  set_f52fe(MMODELGRID,0.);
  set_Te(MMODELGRID,MINTEMP);
  set_TJ(MMODELGRID,MINTEMP);
  set_TR(MMODELGRID,MINTEMP);
  allocate_compositiondata(MMODELGRID);

  /// First pass through to get normalisation coefficients
  for (int n = 0; n < ngrid; n++)
  {
    int mgi = cell[n].modelgridindex;
    rho_sum += get_rhoinit(mgi);
    fe_sum += get_ffe(mgi);

    if (opacity_case == 3)
    {
      if (get_rhoinit(mgi) > 0.)
      {
        if (get_rhoinit(mgi) > rho_crit)
        {
          set_kappagrey(mgi, (0.9 * get_ffe(mgi) + 0.1) * rho_crit / get_rhoinit(mgi));
        }
        else
        {
          set_kappagrey(mgi, (0.9 * get_ffe(mgi) + 0.1));
        }
      }
      else if (get_rhoinit(mgi) == 0.)
      {
        set_kappagrey(mgi,0.);
      }
      else if (get_rhoinit(mgi) < 0.)
      {
        printout("Error: negative density. Abort.\n");
        exit(0);
      }
      opcase3_sum += get_kappagrey(mgi)*get_rhoinit(mgi);
    }
  }

  FILE *grid_file;
  if (rank_global == 0)
  {
    if ((grid_file = fopen("grid.out", "w")) == NULL)
    {
      printf("Cannot open grid file.\n");
      exit(0);
    }
  }

  double check1 = 0.0;
  double check2 = 0.0;
  /// Second pass through allows calculation of normalized kappa_grey
  for (int n = 0; n < ngrid; n++)
  {
    int mgi = cell[n].modelgridindex;
    if (rank_global == 0 && mgi != MMODELGRID) fprintf(grid_file,"%d %d\n",n,mgi); ///write only non emtpy cells to grid file
    if (get_rhoinit(mgi) > 0)
    {
      if (opacity_case == 0)
      {
        set_kappagrey(mgi, GREY_OP);
      }
      else if (opacity_case == 1)
      {
        set_kappagrey(mgi, ((0.9 * get_ffe(mgi)) + 0.1) * GREY_OP / ((0.9 *  mfeg / mtot) + 0.1));
      }
      else if (opacity_case == 2)
      {
        double opcase2_normal = GREY_OP * rho_sum / ((0.9 *  fe_sum) + (0.1 * (ngrid - empty_cells)));
        set_kappagrey(mgi, opcase2_normal/get_rhoinit(mgi) * ((0.9 * get_ffe(mgi)) + 0.1));
      }
      else if (opacity_case == 3)
      {
        double opcase3_normal = GREY_OP * rho_sum / opcase3_sum;
        set_kappagrey(mgi, get_kappagrey(mgi) * opcase3_normal);
      }
      else if (opacity_case == 4)
      {
        ///kappagrey used for initial grey approximation in this case
        set_kappagrey(mgi, ((0.9 * get_ffe(mgi)) + 0.1) * GREY_OP / ((0.9 *  mfeg / mtot) + 0.1));
        //set_kappagrey(mgi, SIGMA_T);
      }
      else
      {
        printout("Unknown opacity case. Abort.\n");
        exit(0);
      }
    }
    else if (get_rhoinit(mgi) == 0.)
    {
      set_kappagrey(mgi, 0.);
    }
    else if (get_rhoinit(mgi) < 0.)
    {
      printout("Error: negative density. Abort.\n");
      exit(0);
    }

    check1 = check1 + (get_kappagrey(mgi)  * get_rhoinit(mgi));
    check2 = check2 + get_rhoinit(mgi);
  }
  if (rank_global == 0) fclose(grid_file);




/*backup copy of old version
  for (n=0;n < ngrid; n++)
  {
  dcen[0] = cell[n].pos_init[0] + (0.5*wid_init);
  dcen[1] = cell[n].pos_init[1] + (0.5*wid_init);
  dcen[2] = cell[n].pos_init[2] + (0.5*wid_init);

  radial_pos = vec_len(dcen);
  if (radial_pos < rmax)
  {
  mkeep =0;
  cell[n].rho_init = cell[n].rho = rho_model[0] * pow( (t_model/tmin), 3.);
  if (opacity_case == 0)
  {
  cell[n].kappa_grey = GREY_OP;
}
  else if (opacity_case == 1)
  {
  cell[n].kappa_grey = ((0.9 * ffegrp_model[0]) + 0.1) * GREY_OP / ((0.9 *  mfeg / mtot) + 0.1);
}
  else if (opacity_case == 2)
  {
  Q = GREY_OP*rho_sum / ((0.9 *  fe_sum) + (0.1 * (npts_model)));
  cell[n].kappa_grey = Q/rho_model[0] * ((0.9 * ffegrp_model[0]) + 0.1);
}
  else
  {
  printout("Unknown opacity case. Abort.\n");
  exit(0);
}
  for (m = 0; m < (npts_model-1); m++)
  {
  if (radial_pos > (vout_model[m] * tmin))
  {
  mkeep = m+1;
  cell[n].rho_init = cell[n].rho = rho_model[m+1] * pow( (t_model/tmin), 3.);
  if (opacity_case == 0)
  {
  cell[n].kappa_grey = GREY_OP;
}
  else if (opacity_case == 1)
  {
  cell[n].kappa_grey = ((0.9 * ffegrp_model[m+1]) + 0.1) * GREY_OP / ((0.9 *  mfeg / mtot) + 0.1);
}
  else if (opacity_case == 2)
  {
  Q = GREY_OP*rho_sum / ((0.9 *  fe_sum) + (0.1 * (npts_model)));
  cell[n].kappa_grey = Q/rho_model[m+1] * ((0.9 * ffegrp_model[m+1]) + 0.1);
}
  else
  {
  printout("Unknown opacity case. Abort.\n");
  exit(0);
}
}
}
  check1 = check1 + (cell[n].kappa_grey  * cell[n].rho_init);
  check2 = check2 + cell[n].rho_init;
  renorm[mkeep] += 1;

}
  else
  {
  cell[n].rho_init = cell[n].rho = 0.0;
  cell[n].kappa_grey = 0.0;
}

  if (cell[n].rho_init < 0)
  {
  printout("Error: negative density. Abort.\n");
  exit(0);
}

}
*/

  /* The following block of code was for a check in the 1d models - currently commented out */

  /* *******************************************************************************************

  check1 = 0.0;
  check2 = 0.0;

  if (renorm[0] > 0 )
  {
  den_norm[0] = rho_model[0] * (pow(vout_model[0],3.)) * 4. * PI * pow(tmin,3.) / 3. / (renorm[0] * wid_init * wid_init * wid_init);
}
  else
  {
  printout("No cells assigned in model shell 0\n");
}

  for (n=1;n<MMODELGRID;n++)
  {
  if (renorm[n] > 0 )
  {
  den_norm[n] = rho_model[n] * (pow(vout_model[n],3.) - pow(vout_model[n-1],3.)) * 4. * PI * pow(tmin,3.) / 3. / (renorm[n] * wid_init * wid_init * wid_init);
}
  else if (n < npts_model)
  {
  printout("No cells assigned in model shell %d\n",n);
}
}

  for (n=0;n < ngrid; n++)
  {
  dcen[0] = cell[n].pos_init[0] + (0.5*wid_init);
  dcen[1] = cell[n].pos_init[1] + (0.5*wid_init);
  dcen[2] = cell[n].pos_init[2] + (0.5*wid_init);

  radial_pos = vec_len(dcen);
  if (radial_pos < rmax)
  {
  mkeep =0;
  cell[n].rho_init = cell[n].rho = den_norm[0] * pow( (t_model/tmin), 3.);
  if (opacity_case == 1)
  {
  cell[n].kappa_grey = ((0.9 *  ffegrp_model[0]) + 0.1) * GREY_OP / ((0.9 *  mfeg / mtot) + 0.1);
}
  else
  {
  cell[n].kappa_grey = GREY_OP;
}
  for (m = 0; m < (npts_model-1); m++)
  {
  if (radial_pos > (vout_model[m] * tmin))
  {
  mkeep = m+1;
  cell[n].rho_init = cell[n].rho = den_norm[m+1] * pow( (t_model/tmin), 3.);
  if (opacity_case == 1)
  {
  cell[n].kappa_grey = ((0.9 *  ffegrp_model[m+1]) + 0.1) * GREY_OP / ((0.9 *  mfeg / mtot) + 0.1);
}
  else
  {
  cell[n].kappa_grey = GREY_OP;
}
}
}
  check1 = check1 + (cell[n].kappa_grey  * cell[n].rho_init);
  check2 = check2 + cell[n].rho_init;
  renorm[mkeep] += 1;

}
  else
  {
  cell[n].rho_init = cell[n].rho = 0.0;
  cell[n].kappa_grey = 0.0;
}

  if (cell[n].rho_init < 0)
  {
  printout("Error: negative density. Abort.\n");
  exit(0);
}

}


  ******************************************************************************************************** */

  printout("Grey normalisation check: %g\n", check1/check2);
  printout("Total mass check: %g\n", check2 * wid_init * wid_init * wid_init / MSUN);

  return 0;
}

///****************************************************************************
/// Routine for doing a density grid read from a 2-D model.
static int density_2d_read(void)
{
  double radial_pos;
  int mkeep1, mkeep2;
  //int renorm[MMODELGRID];
  //double den_norm[MMODELGRID];
  double opcase2_normal;
  int anumber;
  float abundance;
  double zcylindrical, rcylindrical;

  //Calculate the critical opacity at which opacity_case 3 switches from a
  //regime proportional to the density to a regime independent of the density
  //This is done by solving for tau_sobolev == 1
  //tau_sobolev = PI*QE*QE/(ME*C) * rho_crit_para * rho/MNI56 * 3000e-8 * time_step[m].mid;
  rho_crit = ME*CLIGHT*MNI56 / (PI*QE*QE * rho_crit_para * 3000e-8 * tmin);
  printout("grid_init: rho_crit = %g\n", rho_crit);

  double check1 = 0.0;
  double check2 = 0.0;

  rho_sum = 0.0;
  fe_sum = 0.0;
  double opcase3_sum = 0.0;
  int empty_cells = 0;

  for (int n = 0; n < ngrid; n++)
  {
    double dcen[3];
    dcen[0] = cell[n].pos_init[0] + (0.5*wid_init);
    dcen[1] = cell[n].pos_init[1] + (0.5*wid_init);
    dcen[2] = cell[n].pos_init[2] + (0.5*wid_init);

    radial_pos = vec_len(dcen);

    if (radial_pos < rmax)
    {
      mkeep1 = mkeep2 = 0;
      cell[n].modelgridindex = 0;
      zcylindrical = dcen[2];
      dcen[2] = 0.0;
      rcylindrical = vec_len(dcen);

      /*Grid is uniform so only need to search in 1d to get r and z positions */

      for (int m = 0; m < ncoord1_model; m++)
      {
        if (rcylindrical > (m * dcoord1 * tmin/t_model))
        {
          mkeep1 = m;
          //cell[n].modelgridindex = m+1;
        }
      }
      for (int m = 0; m < ncoord2_model; m++)
      {
        if (zcylindrical > (((m * dcoord2) * tmin/t_model) - rmax))
        {
          mkeep2 = m;
          //cell[n].modelgridindex = m+1;
        }
      }
      cell[n].modelgridindex = (mkeep2 * ncoord1_model) + mkeep1;
      modelgrid[cell[n].modelgridindex].initial_radial_pos += radial_pos;

      //renorm[mkeep] += 1;
    }
    else
    {
      cell[n].modelgridindex = MMODELGRID;
    }
  }

  /// Determine the number of simulation cells associated with the model cells
  for (int i = 0; i < npts_model; i++)
  {
    int count = 0;
    for (int ii = 0; ii < ngrid; ii++)
    {
      if (cell[ii].modelgridindex == i) count++;
    }
    modelgrid[i].associated_cells = count;
  }

//  for (n=0;n < ngrid; n++)
//  {
//    mgi = cell[n].modelgridindex;
  for (int mgi = 0; mgi < npts_model; mgi++)
  {
    if (modelgrid[mgi].associated_cells > 0)
    {
      double helper = rho_model[mgi] * pow( (t_model/tmin), 3.);
      set_rhoinit(mgi,helper);
      set_rho(mgi,helper);
      set_ffe(mgi,ffegrp_model[mgi]);
      set_fni(mgi,fni_model[mgi]);
      set_fco(mgi,fco_model[mgi]);
      set_f52fe(mgi,f52fe_model[mgi]);
      set_f48cr(mgi,f48cr_model[mgi]);
      allocate_compositiondata(mgi);
      allocate_cooling(mgi);
      for (int element = 0; element < nelements; element++)
      {
        ///now set the abundances (by mass) of included elements, i.e.
        ///read out the abundances specified in the atomic data file
        anumber = get_element(element);
        abundance = abund_model[mgi][anumber-1];
        //cell[n].composition[element].abundance = abundance;
        modelgrid[mgi].composition[element].abundance = abundance;
        //if (anumber == 8)
        //  cell[n].composition[element].abundance = 1-cell[n].f_ni;
        //if (anumber == 20)
        //  cell[n].composition[element].abundance = 0.1039; ///force Ca to have higher abundance (equal to O abundance) in a pure Ca run

        //printout("element %d has abundance %g in cell %d\n",element,cell[n].composition[element].abundance,n);

        if (anumber == 28)
          set_fnistable(mgi,abundance - get_fni(mgi));

        if (anumber == 27)
          set_fcostable(mgi,abundance - get_fco(mgi));

        if (anumber == 26)
          set_ffestable(mgi,abundance - get_f52fe(mgi));

        if (anumber == 25)
          set_fmnstable(mgi,abundance);

        if (anumber == 24)
          set_fcrstable(mgi,abundance - get_f48cr(mgi));

        if (anumber == 23)
          set_fvstable(mgi,abundance);

        if (anumber == 22)
          set_ftistable(mgi,abundance);
      }

      if (get_rhoinit(mgi) < 0)
      {
        printout("Error: negative density. Abort.\n");
        exit(0);
      }
    }
  }


  /// This is the placeholder for empty cells. Temperatures must be positive
  /// as long as ff opacities are calculated.
  set_rhoinit(MMODELGRID,0.);
  set_rho(MMODELGRID,0.);
  set_nne(MMODELGRID,0.);
  set_ffe(MMODELGRID,0.);
  set_fco(MMODELGRID,0.);
  set_fni(MMODELGRID,0.);
  set_f48cr(MMODELGRID,0.);
  set_f52fe(MMODELGRID,0.);
  set_Te(MMODELGRID,MINTEMP);
  set_TJ(MMODELGRID,MINTEMP);
  set_TR(MMODELGRID,MINTEMP);
  allocate_compositiondata(MMODELGRID);


  /// First pass through to get normalisation coefficients
  for (int n = 0; n < ngrid; n++)
  {
    int mgi = cell[n].modelgridindex;
    rho_sum += get_rhoinit(mgi);
    fe_sum += get_ffe(mgi);

    if (opacity_case == 3)
    {
      if (get_rhoinit(mgi) > 0.)
      {
        if (get_rhoinit(mgi) > rho_crit)
        {
          set_kappagrey(mgi, (0.9 * get_ffe(mgi) + 0.1) * rho_crit/get_rhoinit(mgi));
        }
        else
        {
          set_kappagrey(mgi, (0.9 * get_ffe(mgi) + 0.1));;
        }
      }
      else if (get_rhoinit(mgi) == 0.)
      {
        set_kappagrey(mgi,0.);
      }
      else if (get_rhoinit(mgi) < 0.)
      {
        printout("Error: negative density. Abort.\n");
        exit(0);
      }
      opcase3_sum += get_kappagrey(mgi)*get_rhoinit(mgi);
    }
  }


  FILE *grid_file;
  if (rank_global == 0)
  {
    if ((grid_file = fopen("grid.out", "w")) == NULL)
    {
      printf("Cannot open grid file.\n");
      exit(0);
    }
  }

  /// Second pass through allows calculation of normalized kappa_grey
  for (int n = 0; n < ngrid; n++)
  {
    int mgi = cell[n].modelgridindex;
    if (rank_global == 0 && mgi != MMODELGRID) fprintf(grid_file,"%d %d\n",n,mgi); ///write only non emtpy cells to grid file
    if (get_rhoinit(mgi) > 0)
    {
      if (opacity_case == 0)
      {
        set_kappagrey(mgi, GREY_OP);
      }
      else if (opacity_case == 1)
      {
        set_kappagrey(mgi, ((0.9 * get_ffe(mgi)) + 0.1) * GREY_OP / ((0.9 *  mfeg / mtot) + 0.1));
      }
      else if (opacity_case == 2)
      {
        opcase2_normal = GREY_OP*rho_sum / ((0.9 *  fe_sum) + (0.1 * (ngrid - empty_cells)));
        set_kappagrey(mgi, opcase2_normal/get_rhoinit(mgi) * ((0.9 * get_ffe(mgi)) + 0.1));
      }
      else if (opacity_case == 3)
      {
        opcase3_normal = GREY_OP*rho_sum / opcase3_sum;
        set_kappagrey(mgi, get_kappagrey(mgi) * opcase3_normal);
      }
      else if (opacity_case == 4)
      {
        ///kappagrey used for initial grey approximation in this case
        set_kappagrey(mgi, ((0.9 * get_ffe(mgi)) + 0.1) * GREY_OP / ((0.9 *  mfeg / mtot) + 0.1));
        //set_kappagrey(mgi, SIGMA_T);
      }
      else
      {
        printout("Unknown opacity case. Abort.\n");
        exit(0);
      }
    }
    else if (get_rhoinit(mgi) == 0.)
    {
      set_kappagrey(mgi, 0.);
    }
    else if (get_rhoinit(mgi) < 0.)
    {
      printout("Error: negative density. Abort.\n");
      exit(0);
    }

    check1 = check1 + (get_kappagrey(mgi)  * get_rhoinit(mgi));
    check2 = check2 + get_rhoinit(mgi);
  }
  if (rank_global == 0) fclose(grid_file);


  printout("Grey normalisation check: %g\n", check1/check2);
  printout("Total mass check: %g\n", check2 * wid_init * wid_init * wid_init / MSUN);

  return 0;
}



///****************************************************************************
/// Routine for doing a density grid read from a 3-D model.
static int density_3d_read ()
{
  rho_sum = 0.0;
  fe_sum = 0.0;
  int empty_cells = 0;
  double opcase3_sum = 0;

  //Calculate the critical opacity at which opacity_case 3 switches from a
  //regime proportional to the density to a regime independent of the density
  //This is done by solving for tau_sobolev == 1
  //tau_sobolev = PI*QE*QE/(ME*C) * rho_crit_para * rho/MNI56 * 3000e-8 * time_step[m].mid;
  rho_crit = ME*CLIGHT*MNI56 / (PI*QE*QE * rho_crit_para * 3000e-8 * tmin);
  printout("grid_init: rho_crit = %g\n", rho_crit);

  for (int n = 0; n < ngrid; n++)
  {
    //printout("grid_init: n = %d, grid_type %d, mgi %d, ngrid %d\n", n,grid_type,cell[n].modelgridindex,ngrid);
    double dcen[3];
    dcen[0] = cell[n].pos_init[0] + (0.5*wid_init);
    dcen[1] = cell[n].pos_init[1] + (0.5*wid_init);
    dcen[2] = cell[n].pos_init[2] + (0.5*wid_init);
    //printout("grid_init2: n = %d\n", n);

    double radial_pos = vec_len(dcen);
    modelgrid[cell[n].modelgridindex].initial_radial_pos = radial_pos;
//     printout("grid_init: n = %d, grid_type %d, mgi %d, ngrid %d, r %g, sum %g\n", n,grid_type,cell[n].modelgridindex,ngrid,radial_pos,modelgrid[cell[n].modelgridindex].initial_radial_pos);
//     modelgrid[cell[n].modelgridindex].initial_radial_pos += radial_pos;
//     printout("grid_init: n = %d, grid_type %d, mgi %d, ngrid %d, r %g, sum %g\n", n,grid_type,cell[n].modelgridindex,ngrid,radial_pos,modelgrid[cell[n].modelgridindex].initial_radial_pos);

    /*
    radial_pos = vec_len(dcen);
    if (radial_pos > rmax)
    {
      cell[n].rho_init = cell[n].rho = 0.0;
      cell[n].f_ni = 0.0;
      cell[n].f_fe = 0.0;
      empty_cells += 1;
    }
    else if (cell[n].rho_init < min_den)
    {
      printout("test\n");
      cell[n].rho_init = cell[n].rho = min_den/100.;
      allocate_compositiondata(n);
    }
    else
    {
      allocate_compositiondata(n);
    }

    */
  }

  /*Remainder of routine copied from the 1D case. Previous version removed (but still given below */

  for (int n = 0; n < ngrid; n++)
  {
    int mgi = cell[n].modelgridindex;
    rho_sum += get_rhoinit(mgi);
    fe_sum += get_ffe(mgi);

    if (opacity_case == 3)
    {
      if (get_rhoinit(mgi) > 0.)
      {
        if (get_rhoinit(mgi) > rho_crit)
        {
          set_kappagrey(mgi, (0.9 * get_ffe(mgi) + 0.1) * rho_crit/get_rhoinit(mgi));
        }
        else
        {
          set_kappagrey(mgi, (0.9 * get_ffe(mgi) + 0.1));;
        }
      }
      else if (get_rhoinit(mgi) == 0.)
      {
        set_kappagrey(mgi,0.);
      }
      else if (get_rhoinit(mgi) < 0.)
      {
        printout("Error: negative density. Abort.\n");
        exit(0);
      }
      opcase3_sum += get_kappagrey(mgi)*get_rhoinit(mgi);
    }
  }

  FILE *grid_file;
  if (rank_global == 0)
  {
    if ((grid_file = fopen("grid.out", "w")) == NULL)
    {
      printf("Cannot open grid file.\n");
      exit(0);
    }
  }

  /// Second pass through allows calculation of normalized kappa_grey
  double check1 = 0.0;
  double check2 = 0.0;
  for (int n = 0; n < ngrid; n++)
  {
    int mgi = cell[n].modelgridindex;
    if (rank_global == 0 && mgi != MMODELGRID) fprintf(grid_file,"%d %d\n",n,mgi); ///write only non emtpy cells to grid file
    if (get_rhoinit(mgi) > 0)
    {
      if (opacity_case == 0)
      {
        set_kappagrey(mgi, GREY_OP);
      }
      else if (opacity_case == 1)
      {
        set_kappagrey(mgi, ((0.9 * get_ffe(mgi)) + 0.1) * GREY_OP / ((0.9 *  mfeg / mtot) + 0.1));
      }
      else if (opacity_case == 2)
      {
        double opcase2_normal = GREY_OP*rho_sum / ((0.9 *  fe_sum) + (0.1 * (ngrid - empty_cells)));
        set_kappagrey(mgi, opcase2_normal/get_rhoinit(mgi) * ((0.9 * get_ffe(mgi)) + 0.1));
      }
      else if (opacity_case == 3)
      {
        opcase3_normal = GREY_OP*rho_sum / opcase3_sum;
        set_kappagrey(mgi, get_kappagrey(mgi) * opcase3_normal);
      }
      else if (opacity_case == 4)
      {
        ///kappagrey used for initial grey approximation in this case
        set_kappagrey(mgi, ((0.9 * get_ffe(mgi)) + 0.1) * GREY_OP / ((0.9 *  mfeg / mtot) + 0.1));
        //set_kappagrey(mgi, SIGMA_T);
      }
      else
      {
        printout("Unknown opacity case. Abort.\n");
        exit(0);
      }
    }
    else if (get_rhoinit(mgi) == 0.)
    {
      set_kappagrey(mgi, 0.);
    }
    else if (get_rhoinit(mgi) < 0.)
    {
      printout("Error: negative density. Abort.\n");
      exit(0);
    }

    check1 = check1 + (get_kappagrey(mgi)  * get_rhoinit(mgi));
    check2 = check2 + get_rhoinit(mgi);
  }
  if (rank_global == 0) fclose(grid_file);


  /*
  if (opacity_case == 3)
  {
    for (n = 0; n < ngrid; n++)
    {
      if (cell[n].rho_init > 0)
      {
        if (cell[n].rho_init > rho_crit)
        {
          cell[n].kappa_grey = (0.9 * cell[n].f_fe + 0.1) * rho_crit/cell[n].rho_init;
        }
        else
        {
          cell[n].kappa_grey = (0.9 * cell[n].f_fe + 0.1);
        }
      }
      else if (cell[n].rho_init == 0)
      {
        cell[n].kappa_grey = 0;
      }
      else if (cell[n].rho_init < 0)
      {
        printout("Error: negative density. Abort.\n");
        exit(0);
      }
      opcase3_sum += cell[n].kappa_grey*cell[n].rho_init;
    }
  }
  */


//  printout("grid_init: opcase3_sum, rho_sum %g %g\n",opcase3_sum, rho_sum);


  //MK: second loop to set up opacities which need integrated basic quantities for normalisation

  /*
    for (n=0;n < ngrid; n++)
    {
      if (cell[n].rho_init > 0)
      {
        if (opacity_case == 0)
        {
          cell[n].kappa_grey = GREY_OP;
        }
        else if (opacity_case == 1)
        {
          cell[n].kappa_grey = ((0.9 * cell[n].f_fe) + 0.1) * GREY_OP / ((0.9 *  mfeg / mtot) + 0.1);
        }
        else if (opacity_case == 2)
        {
          opcase2_normal = GREY_OP*rho_sum / ((0.9 *  fe_sum) + (0.1 * (ngrid - empty_cells)));
          cell[n].kappa_grey = opcase2_normal/cell[n].rho_init * ((0.9 * cell[n].f_fe) + 0.1);
        }
        else if (opacity_case == 3)
        {
          opcase3_normal = GREY_OP*rho_sum / opcase3_sum;
          cell[n].kappa_grey = cell[n].kappa_grey * opcase3_normal;
        }
        else if (opacity_case == 4)
        {
          cell[n].kappa_grey = SIGMA_T;
        }
        else
        {
          printout("Unknown opacity case. Abort.\n");
          exit(0);
        }
      }
      else if (cell[n].rho_init == 0)
      {
        cell[n].kappa_grey = 0;
      }
      else if (cell[n].rho_init < 0)
      {
        printout("Error: negative density. Abort.\n");
        exit(0);
      }

      check1 = check1 + (cell[n].kappa_grey  * cell[n].rho_init);
      check2 = check2 + cell[n].rho_init;
    }


  printout("grid_init: opcase3_normal %g\n", opcase3_normal);

  */

  printout("Initial densities taken from readin.\n");
  printout("Grey normalisation check: %g\n", check1/check2);

  return 0;
}


///****************************************************************************
/// Initialise composition dependent cell data for the given cell
void allocate_compositiondata(int modelgridindex)
{
  if ((modelgrid[modelgridindex].composition = (compositionlist_entry *) malloc(nelements*sizeof(compositionlist_entry))) == NULL)
  {
    printout("[fatal] input: not enough memory to initialize compositionlist for cell %d... abort\n",modelgridindex);
    exit(0);
  }

  if ((modelgrid[modelgridindex].nlte_pops = (double *) malloc(total_nlte_levels*sizeof(double))) == NULL)
  {
    printout("[fatal] input: not enough memory to initialize nlte memory for cell %d... abort\n",modelgridindex);
    exit(0);
  }

  for (int element = 0; element < total_nlte_levels; element++)
  {
    modelgrid[modelgridindex].nlte_pops[element] = -1.0; ///flag to indicate that there is
                                                         /// currently no information on the nlte populations
  }

  //printout("Managed to allocate memory for %d nlte levels\n", total_nlte_levels);

  for (int element = 0; element < nelements; element++)
  {
    /// Set initial abundances to zero
    modelgrid[modelgridindex].composition[element].abundance = 0.;

    /// and allocate memory to store the ground level populations for each ionisation stage
    if ((modelgrid[modelgridindex].composition[element].groundlevelpop = (float *) malloc(get_nions(element)*sizeof(float))) == NULL)
    {
      printout("[fatal] input: not enough memory to initialize groundlevelpoplist for element %d in cell %d... abort\n",element,modelgridindex);
      exit(0);
    }
    for (int ion_index = 0; ion_index < get_nions(element); ion_index++)
    {
      modelgrid[modelgridindex].composition[element].groundlevelpop[ion_index]=0.0;
    }

    if ((modelgrid[modelgridindex].composition[element].partfunct = (float *) malloc(get_nions(element)*sizeof(float))) == NULL)
    {
      printout("[fatal] input: not enough memory to initialize partfunctlist for element %d in cell %d... abort\n",element,modelgridindex);
      exit(0);
    }
    /*
    if ((modelgrid[n].composition[element].ltepartfunct = malloc(get_nions(element)*sizeof(float))) == NULL)
    {
      printout("[fatal] input: not enough memory to initialize lte partfunctlist for element %d in cell %d... abort\n",element,n);
      exit(0);
    }
    */
  }
}


///****************************************************************************
/// Initialise composition dependent cell data for the given cell
void allocate_cooling(int modelgridindex)
{
  if ((modelgrid[modelgridindex].cooling = (mgicooling_t *) malloc(nelements*sizeof(mgicooling_t))) == NULL)
  {
    printout("[fatal] input: not enough memory to initialize coolinglist for cell %d... abort\n",modelgridindex);
    exit(0);
  }

  for (int element = 0; element < nelements; element++)
  {
    /// and allocate memory to store the ground level populations for each ionisation stage
    if ((modelgrid[modelgridindex].cooling[element].contrib = (double *) malloc(get_nions(element)*sizeof(double))) == NULL)
    {
      printout("[fatal] input: not enough memory to initialize coolinglist for element %d in cell %d... abort\n",element,modelgridindex);
      exit(0);
    }
  }
}


/*void abundances_setup()
{
  double vec_len(double x[3]);
  int get_elementindex(int anumber);
  double dcen[3],m_r;
  int element,n;

  if (homogeneous_abundances == 1)
  {
    for (n = 0; n < ngrid; n++)
    {
      if (cell[n].rho >= MINDENSITY)
      {
        for (element = 0; element < nelements; element++)
        {
          /// Now set the abundances (by mass) of included elements, i.e.
          /// read out the abundances specified in the atomic data file
          cell[n].composition[element].abundance = elements[element].abundance;

          /// No stable Ni and Co isotopes in that case. Their abundances must fit
          /// the given total Ni and Co abundances and we have no separate readin
          /// for this case!
          //if (get_element(element) == 28)
          //  cell[n].f_ni =  elements[element].abundance - cell[n].f_ni;
          //if (get_element(element) == 27)
          //  cell[n].f_co =  elements[element].abundance - cell[n].f_co;
          if (get_element(element) == 26)
            cell[n].f_fe_init =  elements[element].abundance;
        }
      }
    }
  }
  else
  {
    for (n = 0; n < ngrid; n++)
    {
      if (cell[n].rho >= MINDENSITY)
      {
        dcen[0]=cell[n].pos_init[0] + (0.5*wid_init);
        dcen[1]=cell[n].pos_init[1] + (0.5*wid_init);
        dcen[2]=cell[n].pos_init[2] + (0.5*wid_init);

        m_r = vec_len(dcen) / rmax;
        m_r = pow(m_r,3) * mtot / MSUN;

        if (m_r < 0.5)
        {
          /// Inner part consists of pure Nickel
          cell[n].composition[get_elementindex(28)].abundance = 1.0;
          //cell[n].composition[get_elementindex(27)].abundance = 1e-10;
          //cell[n].composition[get_elementindex(26)].abundance = 1e-10;
          //cell[n].composition[get_elementindex(14)].abundance = 1e-10;
          //cell[n].composition[get_elementindex(8)].abundance = 1e-10;
        }
        else if (m_r < 0.75)
        {
          /// In the middle we assume Nickel & intermediate mass elements
          /// which are represented by Si
          cell[n].composition[get_elementindex(28)].abundance = 1.0 - ((m_r - 0.5)*4);
          cell[n].composition[get_elementindex(14)].abundance = (m_r - 0.5)*4;
          //cell[n].composition[get_elementindex(27)].abundance = 1e-10;
          //cell[n].composition[get_elementindex(26)].abundance = 1e-10;
          //cell[n].composition[get_elementindex(8)].abundance = 1e-10;
        }
        else
        {
          /// For the outer shell we assume a mixture of intermediate mass
          /// elements (represented by Si) and unburned material (represented
          /// by O)
          cell[n].composition[get_elementindex(14)].abundance = 0.5;
          cell[n].composition[get_elementindex(8)].abundance = 0.5;
          //cell[n].composition[get_elementindex(28)].abundance = 1e-10;
          //cell[n].composition[get_elementindex(27)].abundance = 1e-10;
          //cell[n].composition[get_elementindex(26)].abundance = 1e-10;
        }
      }
    }
  }
}*/


static void abundances_3d_read(void)
{
  FILE *abundance_file;
  /// Open the abundances file
  if ((abundance_file = fopen("abundances.txt", "r")) == NULL)
  {
    printout("Cannot open abundances.txt.\n");
    exit(0);
  }

  /// and process through the grid to read in the abundances per cell
  /// The abundance file should only contain information for non-empty
  /// cells. Its format must be cellnumber (integer), abundance for
  /// element Z=1 (float) up to abundance for element Z=30 (float)
  /// i.e. in total one integer and 30 floats.
  for (int n = 0; n < ngrid; n++)
  {
    int mgi = cell[n].modelgridindex;

    float dum[30];
    int cellnumber;
    fscanf(abundance_file, "%d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g", &cellnumber, &dum[0], &dum[1], &dum[2], &dum[3], &dum[4], &dum[5], &dum[6], &dum[7], &dum[8], &dum[9], &dum[10], &dum[11], &dum[12], &dum[13], &dum[14], &dum[15], &dum[16], &dum[17], &dum[18], &dum[19], &dum[20], &dum[21], &dum[22], &dum[23], &dum[24], &dum[25], &dum[26], &dum[27], &dum[28], &dum[29]);

    if (n != cellnumber-1)
    {
      printout("[fatal] abundances_3d_read: grid cell mismatch ... abort\n");
      printout("[fatal] n %d, cellnumber %d\n",n,cellnumber);
      abort();
    }

    for (int element = 0; element < nelements; element++)
    {
      ///now set the abundances (by mass) of included elements, i.e.
      ///read out the abundances specified in the atomic data file
      int anumber = get_element(element);
      float abundance = dum[anumber-1];
      modelgrid[mgi].composition[element].abundance = abundance;

      if (anumber == 28)
      {
        set_fnistable(mgi,abundance - get_fni(mgi));
        //printout("mgi %d, ni_abund %g, fni %g, fnistable %g\n",mgi,abundance,get_fni(mgi),abundance - get_fni(mgi));
      }
      if (anumber == 27)
        set_fcostable(mgi,abundance - get_fco(mgi));

      if (anumber == 26)
        set_ffestable(mgi,abundance - get_f52fe(mgi));

      if (anumber == 25)
      set_fmnstable(mgi,abundance);

      if (anumber == 24)
        set_fcrstable(mgi,abundance - get_f48cr(mgi));

      if (anumber == 23)
        set_fvstable(mgi,abundance);

      if (anumber == 22)
        set_ftistable(mgi,abundance);
      }
    }

  fclose(abundance_file);
}


///***************************************************************************/
static void abundances_1d_read(void)
{
  /// Open the abundances file
  FILE *abundance_file;
  if ((abundance_file = fopen("abundances.txt", "r")) == NULL)
  {
    printout("Cannot open abundances.txt.\n");
    exit(0);
  }

  /// and process through the grid to read in the abundances per cell
  /// The abundance file should only contain information for non-empty
  /// cells. Its format must be cellnumber (integer), abundance for
  /// element Z=1 (float) up to abundance for element Z=30 (float)
  /// i.e. in total one integer and 30 floats.
  for (int n = 0; n < npts_model; n++)
  {
    int cellnumber;
    float dum[30];
    fscanf(abundance_file, "%d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g", &cellnumber, &dum[0], &dum[1], &dum[2], &dum[3], &dum[4], &dum[5], &dum[6], &dum[7], &dum[8], &dum[9], &dum[10], &dum[11], &dum[12], &dum[13], &dum[14], &dum[15], &dum[16], &dum[17], &dum[18], &dum[19], &dum[20], &dum[21], &dum[22], &dum[23], &dum[24], &dum[25], &dum[26], &dum[27], &dum[28], &dum[29]);

    float norm = 0.;
    for (int nn = 0; nn < 30; nn++)
    {
      abund_model[n][nn] = dum[nn];
      norm += dum[nn];
    }

    if (norm > 0)
    {
      for (int nn = 0; nn < 30; nn++)
      {
        abund_model[n][nn] *= 1./norm;
      }
    }
  }

  fclose(abundance_file);
}


///***************************************************************************/
static void assign_temperature(void)
/// Routine for assigning temperatures to the grid cells at the start of the
/// simulation.
{
  /// For a simulation started from scratch we estimate the initial temperatures
  if (!continue_simulation)
  {
    /// We assume that for early times the material is so optically thick, that
    /// all the radiation is trapped in the cell it originates from. This
    /// means furthermore LTE, so that both temperatures can be evaluated
    /// according to the local energy density resulting from the 56Ni decay.
    /// The dilution factor is W=1 in LTE.

    double tstart = time_step[0].mid;
    double factor = CLIGHT/4/STEBO * 1./56/MH * pow(tmin/tstart,3);
    factor *= -1./(tstart*(-TCOBALT+TNICKEL));
    factor *= (-ENICKEL*exp(-tstart/TNICKEL)*tstart*TCOBALT - ENICKEL*exp(-tstart/TNICKEL)*TNICKEL*TCOBALT + ENICKEL*exp(-tstart/TNICKEL)*tstart*TNICKEL + pow(TNICKEL,2)*ENICKEL*exp(-tstart/TNICKEL) - TCOBALT*tstart*ECOBALT*exp(-tstart/TCOBALT) - pow(TCOBALT,2)*ECOBALT*exp(-tstart/TCOBALT) + ECOBALT*tstart*TNICKEL*exp(-tstart/TNICKEL) + pow(TNICKEL,2)*ECOBALT*exp(-tstart/TNICKEL) + ENICKEL*TCOBALT*TNICKEL - ENICKEL*pow(TNICKEL,2) - pow(TNICKEL,2)*ECOBALT + ECOBALT*pow(TCOBALT,2));

    double factor52fe = CLIGHT/4/STEBO * 1./52/MH * pow(tmin/tstart,3);
    factor52fe *= -1./(tstart*(-T52MN+T52FE));
    factor52fe *= (-E52FE*exp(-tstart/T52FE)*tstart*T52MN - E52FE*exp(-tstart/T52FE)*T52FE*T52MN + E52FE*exp(-tstart/T52FE)*tstart*T52FE + pow(T52FE,2)*E52FE*exp(-tstart/T52FE) - T52MN*tstart*E52MN*exp(-tstart/T52MN) - pow(T52MN,2)*E52MN*exp(-tstart/T52MN) + E52MN*tstart*T52FE*exp(-tstart/T52FE) + pow(T52FE,2)*E52MN*exp(-tstart/T52FE) + E52FE*T52MN*T52FE - E52FE*pow(T52FE,2) - pow(T52FE,2)*E52MN + E52MN*pow(T52MN,2));

    double factor48cr = CLIGHT/4/STEBO * 1./48/MH * pow(tmin/tstart,3);
    factor48cr *= -1./(tstart*(-T48V+T48CR));
    factor48cr *= (-E48CR*exp(-tstart/T48CR)*tstart*T48V - E48CR*exp(-tstart/T48CR)*T48CR*T48V + E48CR*exp(-tstart/T48CR)*tstart*T48CR + pow(T48CR,2)*E48CR*exp(-tstart/T48CR) - T48V*tstart*E48V*exp(-tstart/T48V) - pow(T48V,2)*E48V*exp(-tstart/T48V) + E48V*tstart*T48CR*exp(-tstart/T48CR) + pow(T48CR,2)*E48V*exp(-tstart/T48CR) + E48CR*T48V*T48CR - E48CR*pow(T48CR,2) - pow(T48CR,2)*E48V + E48V*pow(T48V,2));

    //factor = CLIGHT/4/STEBO * ENICKEL/56/MH;
    /// This works only for the inbuilt Lucy model
    //factor = CLIGHT/4/STEBO * 3*mtot/4/PI * ENICKEL/56/MH  / pow(vmax,3);
    //for (n = 0; n < ngrid; n++)
    for (int n = 0; n < npts_model; n++)
    {
      //mgi = cell[n].modelgridindex;
      double T_initial = pow(((factor * get_fni(n) * get_rhoinit(n))
           + (factor52fe * get_f52fe(n) * get_rhoinit(n))
           + (factor48cr * get_f48cr(n) * get_rhoinit(n))), 1./4.);
      //T_initial = pow(factor * cell[n].f_ni * cell[n].rho_init * (1.-exp(-tmin/TNICKEL)), 1./4.);
      //T_initial = pow(factor * cell[n].f_ni * (1.-exp(-tmin/TNICKEL))/pow(tmin,3), 1./4.);
      //T_initial = 30615.5;
      if (T_initial < MINTEMP)
      {
        set_Te(n, MINTEMP);
        set_TJ(n, MINTEMP);
        set_TR(n, MINTEMP);
      }
      else if (T_initial > MAXTEMP)
      {
        set_Te(n, MAXTEMP);
        set_TJ(n, MAXTEMP);
        set_TR(n, MAXTEMP);
      }
      else
      {
        set_Te(n, T_initial);
        set_TJ(n, T_initial);
        set_TR(n, T_initial);
      }
      set_W(n, 1.);
    }
  }
  /// For continuation of an existing simulation we read the temperatures
  /// at the end of the simulation and write them to the grid.
  else
  {
    FILE *gridsave_file;
    printout("READIN GRID SNAPSHOT\n");
    if ((gridsave_file = fopen("gridsave.dat", "r")) == NULL)
    {
      printf("[fatal] assign_temperature: Cannot open gridsave.dat.\n");
      exit(0);
    }

    for (int n = 0; n < npts_model; n++)
    {
      //fscanf(inputtemperatures_file,"%d %g %g %g %g %g %g %g\n",&cellnumber,&T_R,&T_e,&W,&T_D,&W_D,&dummy,&dummy);
      //fscanf(inputtemperatures_file,"%d %g %g %g %g %g %g %g %g %d\n",&cellnumber,&T_R,&T_e,&W,&T_D,&W_D,&dummy,&dummy,&dummy,&idummy);
      float T_R,T_e,W,T_J;
      int mgi,thick;
      fscanf(gridsave_file,"%d %g %g %g %g %d",&mgi,&T_R,&T_e,&W,&T_J,&thick);
      if (n == mgi)
      {
        set_TR(mgi, T_R);
        set_Te(mgi, T_e);
        set_W(mgi, W);
        set_TJ(mgi, T_J);
        modelgrid[mgi].thick = thick;

        #ifndef FORCE_LTE
          for (int element = 0; element < nelements; element++)
          {
            int nions = get_nions(element);
            for (int ion = 0; ion < nions; ion++)
            {
              double Gamma;
              fscanf(gridsave_file,"%lg ",&Gamma);
              corrphotoionrenorm[n*nelements*maxion+element*maxion+ion] = Gamma;
            }
          }

          for (int element = 0; element < nelements; element++)
          {
            int nions = get_nions(element);
            for (int ion = 0; ion < nions; ion++)
            {
              double Gamma;
              fscanf(gridsave_file,"%lg ",&Gamma);
              gammaestimator[n*nelements*maxion+element*maxion+ion] = Gamma;
            }
          }
        #endif
      }
      else
      {
        printout("[fatal] assign_temperature: cell mismatch in reading input initial_temperatures.dat ... abort\n");
        printout("[fatal] assign_temperature: read cellnumber %d, expected cellnumber %d\n",mgi,n);
        abort();
      }
    }

    fclose(gridsave_file);
  }
}


static int uniform_grid_setup(void)
/// Routine for doing a uniform cuboidal grid.
{
  int nx = 0;
  int ny = 0;
  int nz = 0;
  for (int n = 0; n < ngrid; n++)
  {
    cell[n].pos_init[0] = -xmax + (2 * nx * xmax / nxgrid);
    cell[n].pos_init[1] = -ymax + (2 * ny * ymax / nygrid);
    cell[n].pos_init[2] = -zmax + (2 * nz * zmax / nzgrid);

    wid_init = 2 * xmax / nxgrid;
    wid_init = 2 * ymax / nygrid;
    wid_init = 2 * zmax / nzgrid;

    //cell[n].cen_init[0] = cell[n].pos_init[0] + (0.5 * wid_init);
    //cell[n].cen_init[1] = cell[n].pos_init[1] + (0.5 * wid_init);
    //cell[n].cen_init[2] = cell[n].pos_init[2] + (0.5 * wid_init);

    cell[n].xyz[0] = nx;
    cell[n].xyz[1] = ny;
    cell[n].xyz[2] = nz;

    nx++;
    if (nx == nxgrid)
    {
      nx = 0;
      ny++;
    }
    if (ny == nygrid)
    {
      ny = 0;
      nz++;
    }

    ///Do we need this initialisation anywhere else (after modelgridindex was initialised) ???????????????????????????????????
    //cell[n].f_ni_stable = 0.0;
    //cell[n].f_co_stable = 0.0;
    //cell[n].f_fe_init = 0.0;
  }

  /*
  /// Finally we must also create the composition dependent data structure for
  /// the samplingcell which is located at MGRID (i.e. the MGRID+1th cell)
  n = MGRID;
  if ((cell[n].composition = malloc(nelements*sizeof(compositionlist_entry))) == NULL)
  {
    printout("[fatal] input: not enough memory to initialize compositionlist for cell %d... abort\n",n);
    exit(0);
  }
  for (element = 0; element < nelements; element++)
  {
    ///now set the abundances (by mass) of included elements, i.e.
    ///read out the abundances specified in the atomic data file
    ///and allocate memory to store the ground level populations for each ionisation stage
    if ((cell[n].composition[element].groundlevelpop = malloc(get_nions(element)*sizeof(float))) == NULL)
    {
      printout("[fatal] input: not enough memory to initialize groundlevelpoplist for element %d in cell %d... abort\n",element,n);
      exit(0);
    }
    if ((cell[n].composition[element].partfunct = malloc(get_nions(element)*sizeof(float))) == NULL)
    {
      printout("[fatal] input: not enough memory to initialize partfunctlist for element %d in cell %d... abort\n",element,n);
      exit(0);
    }
//     if ((cell[n].composition[element].ltepartfunct = malloc(get_nions(element)*sizeof(float))) == NULL)
//     {
//       printout("[fatal] input: not enough memory to initialize lte partfunctlist for element %d in cell %d... abort\n",element,n);
//       exit(0);
//     }
  }
  */

  return 0;
}


int grid_init(void)
/// Subroutine that initialises the grid cells. Designed so that grid cells
/// don't need to be uniform but for the moment they are.
{
  /// Start by checking that the number of grid cells is okay */
  //ngrid = nxgrid * nygrid * nzgrid; ///Moved to input.c
  //if (ngrid > MGRID)
  //{
  //  printout("[fatal] grid_init: Error: too many grid cells. Abort.");
  //  exit(0);
  //}
  for (int n = 0; n <= MMODELGRID; n++)
  {
    modelgrid[n].initial_radial_pos = 0;
  }

  /// The cells will be ordered by x then y, then z. Call a routine that
  /// sets up the initial positions and widths of the cells.
  if (grid_type == GRID_UNIFORM)
  {
    uniform_grid_setup();
  }
  else
  {
    printout("[fatal] grid_init: Error: Unknown grid type. Abort.");
    exit(0);
  }

  /// Now set up the density in each cell.
  if (model_type == RHO_UNIFORM)
  {
    //uniform_density_setup ();
    //abundances_setup();
  }
  else if (model_type == RHO_1D_READ)
  {
    abundances_1d_read();
    density_1d_read();
  }
  else if (model_type == RHO_2D_READ)
  {
    abundances_1d_read(); //for 2d can handle abundances exactly as for 1D
    density_2d_read();
  }
  else if (model_type == RHO_3D_READ)
  {
    density_3d_read();
    abundances_3d_read();
  }
  else
  {
    printout("[fatal] grid_init: Error: Unknown density type. Abort.");
    exit(0);
  }

  /// and assign a temperature to the cells
  assign_temperature();

  /// Finally determine which cells are non-empty...
  /*if ((nonemptycells = malloc(ngrid*sizeof(int))) == NULL)
  {
    printout("[fatal] grid_init: not enough memory to initialize the list of non-empty cells\n");
    exit(0);
  }
  i = 0;
  for (n = 0; n < ngrid; n++)
  {
    if (cell[n].rho_init > MINDENSITY)
    {
      nonemptycells[i] = n;
      i += 1;
    }
  }
  nnonemptycells = i;*/
  /*  if ((nonemptycells = realloc(nonemptycells, nnonemptycells*sizeof(int))) == NULL)
  {
    printout("[fatal] grid_init: problem during reallocation of list of non-empty cells\n");
    exit(0);
  }
  */

  return 0;
}
