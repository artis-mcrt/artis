#!/usr/bin/env zsh

set -x

source ./setupfuncs.sh

modelpath=nltephotospheric_dynamic_ion_range_1d_1dgrid_testrun

getatomicdata atomicdata_hefeconi_fe_i_to_vii.tar.xz

mkdir -p $modelpath

cd $modelpath

rsync -av --exclude="recombrates.txt.zst" ../nebular_1d_3dgrid_inputfiles/ ./

rsync --ignore-times -av ../nltephotospheric_dynamic_ion_range_1d_1dgrid_inputfiles/ ./

tar -xf ../atomicdata_hefeconi_fe_i_to_vii.tar.xz --directory .

# The line walk of the r-packets over 3.7 million lines otherwise takes most of the run time. The NLTE
# solver holds 101 levels of Fe II, so the limit also leaves a superlevel.
setnlevelsmax 150

ln -s ../../ artis

cp artis/artisoptions_nltephotospheric.h artisoptions.h

setoption NUM_PACKETS '1.6e3'

setoption GRID_TYPE_OVERRIDE 'GridType::SPHERICAL1D'

setoption NLTE_TE_NNE_MAXITER '2'

setoption ION_NLEVELS_EXCITED_NLTE '{
  if (element_z == 26 && ionstage == 2) {
    return 100;
  }
  return 50;
}'

setoption NLEVELS_ADD_MISSING_COLLISIONAL_TRANSITIONS '{
  return (element_z < 20) ? 20 : 40;
}'

setoption RATECOEFF_TABLESIZE '40'

setoption MULTIBIN_RADFIELD_FIRST_TIMESTEP '4'

setoption RADFIELDBINCOUNT '24'

setoption DETAILED_BF_ESTIMATORS_FIRST_TIMESTEP '4'

setoption KEEP_ESCAPED_GAMMAS 'true'

# This test covers the output files of each rank.
setoption WRITE_COMBINED_ALLRANK_OUT_FILES 'false'

cd -

set +x
