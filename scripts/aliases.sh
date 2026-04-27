# check mtime of symlinks and compile in parallel
export MAKEFLAGS="--check-symlink-times -j$(nproc --all)"

alias squeueu="squeue --user=$USER --format=\"%.9i %.10P %.8T %.12M %22N %5D %4C %j\""
alias qstatu="squeueu"

# enqueue with a job name set to the current directory name
alias qsub='sbatch --job-name="${PWD##*/}"'

alias qdel="scancel"
alias qstat="squeue"
alias jobinfo="scontrol show jobid"
alias rsyncnooutput="rsync -av --exclude='speclc_angle_res' --exclude='exspec*.txt' --exclude='*.pdf' --exclude='out*.txt' --exclude='*.slurm' --exclude='*.out*' --exclude='packets*' --exclude='logfiles*' --exclude='*.tmp'"
