#!/usr/bin/env zsh

set -x

source ./setupfuncs.sh

modelpath=kilonova_2d_barnesthermalisation_testrun

getatomicdata atomicdata_feconi.tar.xz

mkdir -p $modelpath

cd $modelpath

tar -xf ../atomicdata_feconi.tar.xz --directory ./

# same input files as the other test run
rsync -av ../kilonova_2d_inputfiles/ ./

# the phixsdata_v2.txt of this test, and the checksum files
rsync -av --ignore-times ../kilonova_2d_barnesthermalisation_inputfiles/ ./

ln -s ../../ artis

cp artis/artisoptions_kilonova_lte.h artisoptions.h

setoption NUM_PACKETS '3.2e5'

setoption RATECOEFF_TABLESIZE '20'
setoption MINTEMP '1000.'
setoption MAXTEMP '20000.'

setoption PARTICLE_THERMALISATION_SCHEME 'ParticleThermalisationScheme::BARNES'

setoption GAMMA_THERMALISATION_SCHEME 'GammaThermalisationScheme::BARNES'

# the composition-dependent gamma opacity needs the FREQUENCYDEPENDENT scheme
setoption GAMMA_USE_XCOM_PHOTOION 'false'

setoption KEEP_ESCAPED_GAMMAS 'true'

cd -

set +x
