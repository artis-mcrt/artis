#!/bin/bash
## SLURM META DIRECTIVES HERE DON'T WORK UNDER CENTOS VIRTUAL APPLICATION ENVIRONMENT
## So they are located in artis-virgo-submit.sh as cmd-line parameters to sbatch

export APPTAINER_CONTAINER="/cvmfs/vae.gsi.de/vae26/slurm-25-11/container/user_container-develop.sif"
export APPTAINER_NAME="vae26-user_container"
export APPTAINER_SHARENS=true
export APPTAINER_CONFIGDIR=/tmp/$USER

export LUSTRE_HOME="/lustre/theory/$USER"
export PIXI_HOME="$LUSTRE_HOME/.pixi"
export PIXI_BIN_DIR="$PIXI_HOME/bin"
export PIXI_CACHE_DIR="$LUSTRE_HOME/.cache/rattler"
export UV_INSTALL_DIR="$LUSTRE_HOME/.local/bin"
export UV_CACHE_DIR="$LUSTRE_HOME/.cache/uv"
export UV_PYTHON_DIR="$LUSTRE_HOME/.local/share/uv/python"
export UV_PYTHON_BIN_DIR="$LUSTRE_HOME/.local/bin"
export UV_TOOL_DIR="$LUSTRE_HOME/.local/share/uv/tools"
export UV_TOOL_BIN_DIR="$LUSTRE_HOME/.local/bin"
export PATH="$PIXI_BIN_DIR:$UV_INSTALL_DIR:$UV_TOOL_BIN_DIR:$PATH"

if [ ! -x "$PIXI_BIN_DIR/pixi" ]; then
    curl -fsSL https://pixi.sh/install.sh | bash
fi

if [ ! -x "$PIXI_HOME/envs/gxx/bin/g++" ]; then
    "$PIXI_BIN_DIR/pixi" global install "gxx==16.2"
fi

if ! command -v uv >/dev/null 2>&1; then
    "$PIXI_BIN_DIR/pixi" global install uv
fi

eval `spack load --first --sh openmpi%gcc`
# ARTIS no longer uses GSL. An older version needs both of these lines.
#eval `spack load --first --sh gsl%gcc`
#export LD_LIBRARY_PATH=$(gsl-config --prefix)/lib/:$LD_LIBRARY_PATH

export MAKEFLAGS="--check-symlink-times --jobs=$(nproc --all)"
export OMPI_CXX="$PIXI_HOME/envs/gxx/bin/g++"

# The conda linker of pixi ignores the DT_RPATH of libmpi.so. The option -rpath
# finds the libstdc++ of pixi at run time.
mpi_rpath=$(readelf -d "$(mpicxx --showme:libdirs)/libmpi.so" | awk -F'[][]' '/RPATH|RUNPATH/{print $2}')
export LDFLAGS="-Wl,-rpath-link,$mpi_rpath -Wl,-rpath,$PIXI_HOME/envs/gxx/lib"

cd $SLURM_SUBMIT_DIR

cd artis
make sn3d
cd ..

echo "CPU type: $("$OMPI_CXX" -march=native -Q --help=target | grep -- '-march=  ' | cut -f3)"

# decompress any zipped input files
source ./artis/scripts/exspec-before.sh

hoursleft=$(python3 ./artis/scripts/slurmjobhoursleft.py ${SLURM_JOB_ID})
source ./artis/scripts/corehours-before.sh
echo "$(date): before srun sn3d. hours left: $hoursleft"

time srun -- ./artis/sn3d -w $hoursleft -o ${SLURM_JOB_ID}.slurm > out.txt

hoursleftafter=$(python3 ./artis/scripts/slurmjobhoursleft.py ${SLURM_JOB_ID})
echo "$(date): after srun sn3d finished. hours left: $hoursleftafter"
source ./artis/scripts/corehours-after.sh

if grep -q "RESTART_NEEDED" "output_0-0.txt"
then
    source ./artis/scripts/artis-virgo-submit.sh
else
    if [ -f packets00_0000.out ]; then
        source ./artis/scripts/exspec-zip-virgo-submit.sh
    fi
fi
