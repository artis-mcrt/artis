#!/usr/bin/env zsh

set -x

runfolder=kilonova_2d_timedepnlte_testrun

if [ ! -f atomicdata_sryzrlace.tar.zst ]; then curl -O -L https://github.com/artis-mcrt/artis/releases/download/v2026.5.15/atomicdata_sryzrlace.tar.zst; fi

mkdir -p $runfolder

cd $runfolder

# the model files of the kilonova_2d test. The atomic data archive supplies compositiondata.txt
rsync -av --exclude="compositiondata.txt" --exclude="input-*.txt" --exclude="results_md5*" ../kilonova_2d_inputfiles/ ./

tar -xf ../atomicdata_sryzrlace.tar.zst --directory ./

rsync -av ../kilonova_2d_timedepnlte_inputfiles/ ./

ln -s ../../ artis

cp artis/artisoptions_nltenebular.h artisoptions.h

xz -f -d -v -T0 *.xz

sed -i'' -e 's/constexpr int MPKTS.*/constexpr int MPKTS = 200000;/g' artisoptions.h

sed -i'' -e 's/constexpr int TABLESIZE.*/constexpr int TABLESIZE = 20;/g' artisoptions.h
sed -i'' -e 's/constexpr double MINTEMP.*/constexpr double MINTEMP = 1000.;/g' artisoptions.h
sed -i'' -e 's/constexpr double MAXTEMP.*/constexpr double MAXTEMP = 20000.;/g' artisoptions.h

perl -0777 -i -pe 's|constexpr int ION_NLEVELS_EXCITED_NLTE\(int element_z, int ionstage\) \{.*?\n\}|constexpr int ION_NLEVELS_EXCITED_NLTE(int element_z, int ionstage) { return 30; }|s' artisoptions.h

sed -i'' -e 's/constexpr int FIRST_NLTE_RADFIELD_TIMESTEP.*/constexpr int FIRST_NLTE_RADFIELD_TIMESTEP = 2;/g' artisoptions.h
sed -i'' -e 's/constexpr int DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP.*/constexpr int DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP = 2;/g' artisoptions.h

sed -i'' -e 's/constexpr std::optional<int> NLTE_TIME_DEPENDENT_FIRST_TIMESTEP.*/constexpr std::optional<int> NLTE_TIME_DEPENDENT_FIRST_TIMESTEP = 3;/g' artisoptions.h

sed -i'' -e 's/constexpr int NLTE_OUTER_ANDERSON_DEPTH.*/constexpr int NLTE_OUTER_ANDERSON_DEPTH = 2;/g' artisoptions.h
sed -i'' -e 's/constexpr double NLTE_OUTER_RELTOL.*/constexpr double NLTE_OUTER_RELTOL = 0.01;/g' artisoptions.h

sed -i'' -e 's/constexpr bool USE_CALCULATED_MEANATOMICWEIGHT.*/constexpr bool USE_CALCULATED_MEANATOMICWEIGHT = true;/g' artisoptions.h
sed -i'' -e 's/constexpr auto PARTICLE_THERMALISATION_SCHEME.*/constexpr auto PARTICLE_THERMALISATION_SCHEME = ParticleThermalisationScheme::TIMEDEPENDENT;/g' artisoptions.h

sed -i'' -e 's/constexpr bool WRITE_EMISSIONABSORPTION_SPEC_AT_END.*/constexpr bool WRITE_EMISSIONABSORPTION_SPEC_AT_END = true;/g' artisoptions.h

cd -

set +x
