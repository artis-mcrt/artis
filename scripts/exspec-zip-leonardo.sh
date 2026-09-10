#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH --exclusive
#SBATCH --partition=dcgp_usr_prod
#SBATCH --qos=normal
#SBATCH --account=EUHPC_R07_209
#SBATCH --mail-type=ALL

projectfolder="/leonardo_work/EUHPC_R07_209"

export PATH="$projectfolder/bin:$PATH"
export PIXI_HOME="$projectfolder/.pixi"
export PIXI_BIN_DIR="$projectfolder/.pixi/bin"
export PIXI_CACHE_DIR="$projectfolder/.cache/rattler"
export UV_CACHE_DIR="$projectfolder/.cache/uv"
export UV_PYTHON_DIR="$projectfolder/.local/share/uv/python"
export UV_PYTHON_BIN_DIR="$projectfolder/.local/bin"
export UV_TOOL_DIR="$projectfolder/.local/share/uv/tools"
export UV_TOOL_BIN_DIR="$projectfolder/.local/bin"
export PATH="$PIXI_BIN_DIR:$UV_TOOL_BIN_DIR:$PATH"

# See artis-leonardo.sh for the reason of the module and of the two exports.
module load openmpi/4.1.6--gcc--12.2.0-cuda-12.2
module list

export OMPI_CXX=g++
export LDFLAGS="-Wl,-rpath-link,/usr/lib64"

cd $SLURM_SUBMIT_DIR

echo "CPU type: $(c++ -march=native -Q --help=target | grep -- '-march=  ' | cut -f3)"

source ./artis/scripts/corehours-before.sh
echo "$(date): before exspec"

source ./artis/scripts/run-exspec-if-needed.sh

source ./artis/scripts/exspec-after.sh

echo "$(date): after exspec finished"
source ./artis/scripts/corehours-after.sh
