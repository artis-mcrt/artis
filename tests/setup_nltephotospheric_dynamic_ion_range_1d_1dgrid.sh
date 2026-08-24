#!/usr/bin/env zsh

set -x

runfolder=nltephotospheric_dynamic_ion_range_1d_1dgrid_testrun

if [ ! -f atomicdata_hefeconi_fe_i_to_vii.tar.xz ]; then curl -O -L https://github.com/artis-mcrt/artis/releases/download/v2026.5.15/atomicdata_hefeconi_fe_i_to_vii.tar.xz; fi

mkdir -p $runfolder

cd $runfolder

rsync -av --exclude="recombrates.txt" ../nebular_1d_3dgrid_inputfiles/ ./

rsync --ignore-times -av ../nltephotospheric_dynamic_ion_range_1d_1dgrid_inputfiles/ ./

tar -xf ../atomicdata_hefeconi_fe_i_to_vii.tar.xz --directory .

ln -s ../../ artis

cp artis/artisoptions_nltephotospheric_dynamic_ion_range.h artisoptions.h

sed -i.bak -e 's/constexpr int MPKTS.*/constexpr int MPKTS = 400;/g' artisoptions.h

sed -i.bak -e 's/constexpr std::optional<GridType> GRID_TYPE_OVERRIDE.*/constexpr std::optional<GridType> GRID_TYPE_OVERRIDE = GridType::SPHERICAL1D;/g' artisoptions.h

sed -i.bak -e 's/constexpr int NLTE_TE_NNE_MAXITER.*/constexpr int NLTE_TE_NNE_MAXITER = 2;/g' artisoptions.h

perl -0777 -i -pe 'my $n = s|^constexpr int ION_NLEVELS_EXCITED_NLTE\(int element_z, int ionstage\) \{.*?^\}$|constexpr int ION_NLEVELS_EXCITED_NLTE(int element_z, int ionstage) {\n  if (element_z == 26 && ionstage == 2) {\n    return 100;\n  }\n  return 50;\n}|ms; die "[error] the pattern for ION_NLEVELS_EXCITED_NLTE did not match once\n" unless $n == 1;' artisoptions.h

sed -i.bak -e 's|constexpr int NLEVELS_REQUIRETRANSITIONS(int element_z, int ionstage) {.*}|constexpr int NLEVELS_REQUIRETRANSITIONS(int element_z, int ionstage) { return (element_z < 20) ? 20 : 40; }|g' artisoptions.h

sed -i.bak -e 's/constexpr int RATECOEFF_TABLESIZE.*/constexpr int RATECOEFF_TABLESIZE = 40;/g' artisoptions.h

sed -i.bak -e 's/constexpr int FIRST_NLTE_RADFIELD_TIMESTEP.*/constexpr int FIRST_NLTE_RADFIELD_TIMESTEP = 4;/g' artisoptions.h

sed -i.bak -e 's/constexpr int RADFIELDBINCOUNT.*/constexpr int RADFIELDBINCOUNT = 24;/g' artisoptions.h

sed -i.bak -e 's/constexpr int DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP.*/constexpr int DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP = 4;/g' artisoptions.h

rm -f artisoptions.h.bak

cd -

set +x
