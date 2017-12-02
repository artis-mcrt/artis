#include "sn3d.h"
#include "atomic.h"
#include "ltepop.h"

extern inline int get_element(int element);
extern inline int get_elementindex(int Z);
extern inline int get_nions(int element);
extern inline int get_ionstage(int element, int ion);
extern inline int get_nlevels(int element, int ion);
extern inline int get_nlevels_nlte(int element, int ion);
extern inline int get_ionisinglevels(int element, int ion);
extern inline double epsilon(int element, int ion, int level);
extern inline double stat_weight(int element, int ion, int level);
extern inline int get_bfcontinua(int element, int ion);
extern inline bool is_nlte(int element, int ion, int level);
extern inline int get_continuumindex(int element, int ion, int level);
extern inline int get_nphixstargets(int element, int ion, int level);
extern inline int get_phixsupperlevel(int element, int ion, int level, int phixstargetindex);
extern inline double get_phixsprobability(int element, int ion, int level, int phixstargetindex);
extern inline double einstein_spontaneous_emission(int lineindex);
extern inline double osc_strength(int lineindex);
extern inline double get_coll_str(int lineindex);
extern inline double statw_upper(int lineindex);
extern inline double statw_lower(int lineindex);
extern inline double photoionization_crosssection_macroatom(double nu_edge, double nu);
extern inline double photoionization_crosssection(int element, int ion, int level, double nu_edge, double nu);
extern inline double get_phixs_threshold(int element, int ion, int level, int phixstargetindex);


double get_tau_sobolev(int modelgridindex, int lineindex, double t_current)
{
  const int element = linelist[lineindex].elementindex;
  const int ion = linelist[lineindex].ionindex;
  const int lower = linelist[lineindex].lowerlevelindex;
  const int upper = linelist[lineindex].upperlevelindex;

  const double statweight_target = statw_upper(lineindex);
  const double statweight_lower = statw_lower(lineindex);

  const double n_l = get_levelpop(modelgridindex,element,ion,lower);
  const double n_u = get_levelpop(modelgridindex,element,ion,upper);

  const double nu_trans = (epsilon(element, ion, upper) - epsilon(element, ion, lower)) / H;
  const double A_ul = einstein_spontaneous_emission(lineindex);
  const double B_ul = CLIGHTSQUAREDOVERTWOH / pow(nu_trans, 3) * A_ul;
  const double B_lu = statweight_target / statweight_lower * B_ul;

  const double tau_sobolev = (B_lu * n_l - B_ul * n_u) * HCLIGHTOVERFOURPI * t_current;
  return tau_sobolev;
}
