#ifndef RPKT_H
#define RPKT_H

#include "artisoptions.h"
#include "constants.h"
#include "grid.h"
#include "sn3d.h"

void do_rpkt(struct packet *pkt_ptr, double t2);
void emit_rpkt(struct packet *pkt_ptr);
auto closest_transition(double nu_cmf, int next_trans) -> int;
auto calculate_chi_bf_gammacontr(int modelgridindex, double nu) -> double;
void calculate_chi_rpkt_cont(double nu_cmf, struct rpkt_continuum_absorptioncoeffs *chi_rpkt_cont_thisthread,
                             int modelgridindex, bool usecellhistupdatephixslist);

constexpr auto get_linedistance(const double prop_time, const double nu_cmf, const double nu_trans,
                                const double d_nu_on_d_l) -> double {
  // distance from packet position to redshifting into line at frequency nu_trans

  if (nu_cmf <= nu_trans) {
    return 0;  /// photon was propagated too far, make sure that we don't miss a line
  }

  if constexpr (USE_RELATIVISTIC_DOPPLER_SHIFT) {
    // With special relativity, the Doppler shift formula has an extra factor of 1/gamma in it,
    // which changes the distance reach a line resonance and creates a dependence
    // on packet position and direction

    // use linear interpolation of frequency along the path
    return (nu_trans - nu_cmf) / d_nu_on_d_l;
  }

  return CLIGHT * prop_time * (nu_cmf / nu_trans - 1);
}

[[nodiscard]] inline auto get_ionestimindex(const int mgi, const int element, const int ion) -> int {
  assert_testmodeonly(ion >= 0);
  assert_testmodeonly(ion < get_nions(element));
  assert_always(grid::get_numassociatedcells(mgi) > 0);
  // const int nonemptymgi = grid::get_modelcell_nonemptymgi(mgi);
  return mgi * get_includedions() + get_uniqueionindex(element, ion);
}
#endif  // RPKT_H
