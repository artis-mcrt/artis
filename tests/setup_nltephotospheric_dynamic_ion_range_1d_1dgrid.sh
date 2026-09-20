#!/usr/bin/env zsh

set -x

source ./setupfuncs.sh

runfolder=nltephotospheric_dynamic_ion_range_1d_1dgrid_testrun

getatomicdata atomicdata_hefeconi_fe_i_to_vii.tar.xz

mkdir -p $runfolder

cd $runfolder

rsync -av --exclude="recombrates.txt" ../nebular_1d_3dgrid_inputfiles/ ./

rsync --ignore-times -av ../nltephotospheric_dynamic_ion_range_1d_1dgrid_inputfiles/ ./

tar -xf ../atomicdata_hefeconi_fe_i_to_vii.tar.xz --directory .

ln -s ../../ artis

cp artis/artisoptions_nltephotospheric_dynamic_ion_range.h artisoptions.h

sedopt "constexpr std::int64_t NUM_PACKETS.*" "constexpr std::int64_t NUM_PACKETS = 1600;"

sedopt 'constexpr std::optional<GridType> GRID_TYPE_OVERRIDE.*' 'constexpr std::optional<GridType> GRID_TYPE_OVERRIDE = GridType::SPHERICAL1D;'

sedopt 'constexpr int NLTE_TE_NNE_MAXITER.*' 'constexpr int NLTE_TE_NNE_MAXITER = 2;'

perl -0777 -i -pe 'my $n = s|^constexpr int ION_NLEVELS_EXCITED_NLTE\(int element_z, int ionstage\) \{.*?^\}$|constexpr int ION_NLEVELS_EXCITED_NLTE(int element_z, int ionstage) {\n  if (element_z == 26 && ionstage == 2) {\n    return 100;\n  }\n  return 50;\n}|ms; die "[error] the pattern for ION_NLEVELS_EXCITED_NLTE did not match once\n" unless $n == 1;' artisoptions.h

sedopt 'constexpr int NLEVELS_REQUIRETRANSITIONS(int element_z, int ionstage) {.*}' 'constexpr int NLEVELS_REQUIRETRANSITIONS(int element_z, int ionstage) { return (element_z < 20) ? 20 : 40; }'

sedopt 'constexpr int RATECOEFF_TABLESIZE.*' 'constexpr int RATECOEFF_TABLESIZE = 40;'

sedopt 'constexpr int FIRST_NLTE_RADFIELD_TIMESTEP.*' 'constexpr int FIRST_NLTE_RADFIELD_TIMESTEP = 4;'

sedopt 'constexpr int RADFIELDBINCOUNT.*' 'constexpr int RADFIELDBINCOUNT = 24;'

sedopt 'constexpr int DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP.*' 'constexpr int DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP = 4;'

rm -f artisoptions.h.bak

cd -

set +x
