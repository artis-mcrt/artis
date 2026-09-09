#!/bin/bash -x
# Submit the ARTIS job on Leonardo. Use this script in place of a direct sbatch.

# sbatch expands no variable in a #SBATCH line, and Slurm supplies no
# SBATCH_MAIL_USER variable. The address must therefore come from the command
# line. The job script artis-leonardo.sh keeps all of the resource options,
# because a #SBATCH line works correctly on this system.

# Leonardo mails <user>@leonardo.local when the address is absent, and that
# address reaches nobody. So give no --mail-user when EMAIL is empty.
mailargs=()
if [ -n "$EMAIL" ]; then
    mailargs=(--mail-user="$EMAIL")
fi

sbatch -J "$(basename "$(pwd)")" "${mailargs[@]}" -- artis/scripts/artis-leonardo.sh

# Add a line like this to your .bashrc to get the notifications:
# export EMAIL=your_email_address
