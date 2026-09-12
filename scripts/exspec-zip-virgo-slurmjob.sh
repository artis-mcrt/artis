#!/bin/bash

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
export UV_PYTHON_INSTALL_DIR="$LUSTRE_HOME/.local/share/uv/python"
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

eval `spack load --first --sh openmpi%gcc`
# ARTIS no longer uses GSL. An older version needs both of these lines.
#eval `spack load --first --sh gsl%gcc`
#export LD_LIBRARY_PATH=$(gsl-config --prefix)/lib/:$LD_LIBRARY_PATH

export MAKEFLAGS="--check-symlink-times --jobs=${SLURM_CPUS_PER_TASK:-$(nproc)}"
export OMPI_CXX="$PIXI_HOME/envs/gxx/bin/g++"

# The conda linker of pixi ignores the DT_RPATH of libmpi.so. The option -rpath
# finds the libstdc++ of pixi at run time.
mpi_rpath=$(readelf -d "$(mpicxx --showme:libdirs)/libmpi.so" | awk -F'[][]' '/RPATH|RUNPATH/{print $2}')
export LDFLAGS="-Wl,-rpath-link,$mpi_rpath -Wl,-rpath,$PIXI_HOME/envs/gxx/lib"

cd "${SLURM_SUBMIT_DIR:?}"

cd artis
make exspec || exit 1
cd ..

echo "CPU type: $("$OMPI_CXX" -march=native -Q --help=target | grep -- '-march=  ' | cut -f3)"


source ./artis/scripts/corehours-before.sh
echo "$(date): before exspec"

source ./artis/scripts/run-exspec-if-needed.sh

source ./artis/scripts/exspec-after.sh

echo "$(date): after exspec finished"
source ./artis/scripts/corehours-after.sh
