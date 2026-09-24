#!/usr/bin/env zsh

set -x

source ./setupfuncs.sh

runfolder=kilonova_1d_timedepnlte_testrun

getatomicdata atomicdata_sryzrlace.tar.zst

mkdir -p $runfolder

cd $runfolder

# the model files of the kilonova_1d test. The atomic data archive supplies compositiondata.txt, and
# the second rsync below supplies the input and checksum files of this test
rsync -av --exclude="compositiondata.txt" ../kilonova_1d_inputfiles/ ./

tar -xf ../atomicdata_sryzrlace.tar.zst --directory ./

# --ignore-times: the input files of this test have the same size as the files of the kilonova_1d
# test. After a checkout they also have the same modification time. Without the option rsync skips them.
rsync --ignore-times -av ../kilonova_1d_timedepnlte_inputfiles/ ./

ln -s ../../ artis

cp artis/artisoptions_kilonova_nlte.h artisoptions.h

xz -f -d -v -T0 *.xz

# Scale the model to a typical kilonova: ten times the mass and one third of the velocity of the
# kilonova_1d model, which gives 0.044 Msun below 0.16 c. The shape of the density profile does not
# change. The model then holds the density of a 0.05 Msun kilonova below 0.15 c to inside 13 percent
# over the whole time range of the test. The thin fast model of the kilonova_1d test is about 300
# times less dense at the same time, and the thermal balance of such a model has no solution.
awk 'BEGIN{n=0; logscale=log(27*10)/log(10)} /^#/{print; next} {n++; if(n<=2){print; next} if(NF>4){$2=sprintf("%.6f",$2/3); $3=sprintf("%.8f",$3+logscale)}; print}' model.txt > model_scaled.txt
mv model_scaled.txt model.txt

# the 1D model has 25 cells and the 2D model had 128, so fewer packets give the same number of
# packets for each cell
sedopt "constexpr std::int64_t NUM_PACKETS.*" "constexpr std::int64_t NUM_PACKETS = 80'000;"

sedopt 'constexpr int RATECOEFF_TABLESIZE.*' 'constexpr int RATECOEFF_TABLESIZE = 40;'

# a few near-vacuum cells reach NLTE_TE_NNE_MAXITER in every timestep, so this limit keeps the run time of the test low
sedopt 'constexpr int NLTE_TE_NNE_MAXITER.*' 'constexpr int NLTE_TE_NNE_MAXITER = 10;'

# element_z == 58 is cerium. This test gives cerium no NLTE level, so that it tests the hybrid mode:
# the NLTE solver holds Sr, Y, Zr, and La, and calculate_ion_balance_nne() holds cerium with the
# photoionisation balance. The mode covers the ion balance of both kinds of element in one cell, the
# sum of the electron contributions, and the rule of chargetransfer.cc that a reaction needs NLTE
# levels on both sides. Keep at least one element at zero, or the test loses the hybrid mode.
sedopt 'constexpr int ION_NLEVELS_EXCITED_NLTE.*' 'constexpr int ION_NLEVELS_EXCITED_NLTE(int element_z, int ionstage) { return (element_z == 58) ? 0 : 20; }'

perl -0777 -i -pe 'my $n = s|^constexpr int NLEVELS_REQUIRETRANSITIONS\(int element_z, int ionstage\) \{.*?^\}$|constexpr int NLEVELS_REQUIRETRANSITIONS(int element_z, int ionstage) { return 10; }|ms; die "[error] the pattern for NLEVELS_REQUIRETRANSITIONS did not match once\n" unless $n == 1;' artisoptions.h

sedopt 'constexpr int FIRST_NLTE_RADFIELD_TIMESTEP.*' 'constexpr int FIRST_NLTE_RADFIELD_TIMESTEP = 2;'
sedopt 'constexpr int DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP.*' 'constexpr int DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP = 2;'

sedopt 'constexpr bool NLTE_USE_GTH_SOLVER.*' 'constexpr bool NLTE_USE_GTH_SOLVER = false;'

sedopt 'constexpr std::optional<int> NLTE_TIME_DEPENDENT_FIRST_TIMESTEP.*' 'constexpr std::optional<int> NLTE_TIME_DEPENDENT_FIRST_TIMESTEP = 3;'

sedopt 'constexpr bool ENABLE_CHARGE_TRANSFER_REACTIONS.*' 'constexpr bool ENABLE_CHARGE_TRANSFER_REACTIONS = true;'

rm -f artisoptions.h.bak

cd -

set +x
