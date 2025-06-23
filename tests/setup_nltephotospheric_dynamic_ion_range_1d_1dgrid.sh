#!/usr/bin/env zsh

set -x

runfolder=nltephotospheric_dynamic_ion_range_1d_1dgrid_testrun

mkdir -p $runfolder

if [ ! -f atomicdata_hefeconi_fe_I_to_VII.tar.xz ]; then curl -O https://theory.gsi.de/~lshingle/artis_http_public/artis/atomicdata_hefeconi_fe_i_to_vii.tar.xz; fi

tar -xf atomicdata_hefeconi_fe_i_to_vii.tar.xz --directory $runfolder

rsync -av --ignore-times nltephotospheric_dynamic_ion_range_1d_1dgrid_inputfiles/ $runfolder

rsync -av nebular_1d_3dgrid_inputfiles/model.txt $runfolder

cp ../data/* $runfolder

cp ../artisoptions_nltephotospheric_dynamic_ion_range.h $runfolder/artisoptions.h

cd $runfolder

sed -i'' -e 's/constexpr int MPKTS.*/constexpr int MPKTS = 20;/g' artisoptions.h

sed -i'' -e 's/constexpr auto GRID_TYPE.*/constexpr auto GRID_TYPE = GridType::SPHERICAL1D;/g' artisoptions.h

sed -i'' -e 's|constexpr int NLEVELS_REQUIRETRANSITIONS(int Z, int ionstage) {.*}|constexpr int NLEVELS_REQUIRETRANSITIONS(int Z, int ionstage) { return (Z < 20) ? 20 : 40; }|g' artisoptions.h


sed -i'' -e 's/constexpr int TABLESIZE.*/constexpr int TABLESIZE = 20;/g' artisoptions.h
sed -i'' -e 's/constexpr double MINTEMP.*/constexpr double MINTEMP = 3500.;/g' artisoptions.h
sed -i'' -e 's/constexpr double MAXTEMP.*/constexpr double MAXTEMP = 140000.;/g' artisoptions.h

sed -i'' -e 's/constexpr int FIRST_NLTE_RADFIELD_TIMESTEP.*/constexpr int FIRST_NLTE_RADFIELD_TIMESTEP = 4;/g' artisoptions.h

sed -i'' -e 's/constexpr int DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP.*/constexpr int DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP = 4;/g' artisoptions.h

sed -i'' -e 's/constexpr float STRICT_POPULATION_CHECKING_INVERSION_FACTOR_SOLVER_FAIL.*/constexpr float STRICT_POPULATION_CHECKING_INVERSION_FACTOR_SOLVER_FAIL = 50.0;/g' artisoptions.h

sed -i'' -e 's/constexpr float STRICT_POPULATION_CHECKING_INVERSION_FACTOR_PRINTOUT_WARNING.*/constexpr float STRICT_POPULATION_CHECKING_INVERSION_FACTOR_PRINTOUT_WARNING = 2.0;/g' artisoptions.h

sed -i'' -e 's/constexpr double NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION.*/constexpr double NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION = 1e-9;/g' artisoptions.h

sed -i'' -e 's/constexpr bool WRITE_PARTIAL_EMISSIONABSORPTIONSPEC.*/constexpr bool WRITE_PARTIAL_EMISSIONABSORPTIONSPEC = true;/g' artisoptions.h

cd -

set +x
