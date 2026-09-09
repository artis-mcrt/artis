#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --ntasks=1008
#SBATCH --ntasks-per-node=112
#SBATCH --exclusive
#SBATCH --cpus-per-task=1
#SBATCH --partition=dcgp_usr_prod
#SBATCH --qos=normal
#SBATCH --account=EUHPC_R07_209
#SBATCH --mail-type=ALL
# artis-leonardo-submit.sh gives --mail-user from the EMAIL variable

export PATH="/leonardo_work/EUHPC_R07_209/bin:$PATH"
export PIXI_HOME="/leonardo_work/EUHPC_R07_209/.pixi"
export PIXI_BIN_DIR="/leonardo_work/EUHPC_R07_209/.pixi/bin"
export PIXI_CACHE_DIR="/leonardo_work/EUHPC_R07_209/.cache/rattler"
export UV_CACHE_DIR="/leonardo_work/EUHPC_R07_209/.cache/uv"
export UV_PYTHON_DIR="/leonardo_work/EUHPC_R07_209/.local/share/uv/python"
export UV_PYTHON_BIN_DIR="/leonardo_work/EUHPC_R07_209/.local/bin"
export UV_TOOL_DIR="/leonardo_work/EUHPC_R07_209/.local/share/uv/tools"
export UV_TOOL_BIN_DIR="/leonardo_work/EUHPC_R07_209/.local/bin"
export PATH="$PIXI_BIN_DIR:$UV_TOOL_BIN_DIR:$PATH"

# Load the Open MPI of the system. It has the UCX transport of the InfiniBand
# network. A conda-forge Open MPI has no UCX. It falls back to TCP and the run
# then stops with "received unexpected process identifier". The module also
# gives zstd and xz for the compressed files. Load the module after the exports
# above, because the module must put its own mpicxx first in the PATH.
module load openmpi/4.1.6--gcc--12.2.0-cuda-12.2
module list

# The Open MPI of the system is a build of gcc 12.2.0, which is too old for
# ARTIS. Compile the C++ with the conda-forge gcc instead. The MPI library has
# a C interface, so the two compilers are compatible.
export OMPI_CXX=g++

# libmpi.so of the system needs the libevent and the hwloc of /usr/lib64. The
# linker of conda-forge does not search that folder, so name it here.
export LDFLAGS="-Wl,-rpath-link,/usr/lib64"

cd $SLURM_SUBMIT_DIR

export MAKEFLAGS="--check-symlink-times --jobs=$(nproc --all)"
cd artis
make
cd ..

mpicxx --version
echo "CPU type: $(c++ -march=native -Q --help=target | grep -- '-march=  ' | cut -f3)"

# decompress any zipped input files
source ./artis/scripts/exspec-before.sh
hoursleft=$(python3.14 ./artis/scripts/slurmjobhoursleft.py ${SLURM_JOB_ID})
source ./artis/scripts/corehours-before.sh
echo "$(date): before srun sn3d. hours left: $hoursleft"
time srun --hint=nomultithread -- ./artis/sn3d -w $hoursleft -o ${SLURM_JOB_ID}.slurm > out.txt
hoursleftafter=$(python3.14 ./artis/scripts/slurmjobhoursleft.py ${SLURM_JOB_ID})
echo "$(date): after srun sn3d finished. hours left: $hoursleftafter"
source ./artis/scripts/corehours-after.sh

if grep -q "RESTART_NEEDED" "output_0-0.txt"
then
    # the submit script sets the job name and the address of the next job
    source ./artis/scripts/artis-leonardo-submit.sh
else
    # post-processing can remove restart files, so only queue it when no continuation job was submitted
    if [ -f packets00_0000.out ]; then
        source ./artis/scripts/exspec-zip-leonardo-submit.sh
    fi
fi
