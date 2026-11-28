#!/usr/bin/env bash
# Remove the files that a simulation writes into the current model folder. Then copy input-newrun.txt to input.txt.

if [ $# -gt 1 ] || { [ $# -eq 1 ] && [ "$1" != -y ]; }; then
  echo "Usage: artis/scripts/clean.sh [-y], in the model folder. With -y, the script asks no question." >&2
  exit 1
fi

# for the +([0-9]) patterns that only match generated filenames (e.g. output_0-0.txt but not output_0-0_backup.txt)
shopt -s extglob

paths="input.txt.tmp input-newrun.txt.tmp gridsave*.tmp packets*.tmp vspecpol*.tmp vpackets*.tmp vpkt_grid*.tmp *.out *.out.* out.txt?(.zst) output_+([0-9])-+([0-9]).txt?(.zst|.gz|.xz) exspec*.txt* machine.file.* core.* job_from_ts+([0-9]) inputcheck packets vspecpol vpackets vpkt_grid speclc_angle_res bflist.dat ratecoeff.dat line_list.txt syn_dir.txt logfiles.tar*"

# Older job scripts made a job folder with the name of the Slurm job ID, e.g. 12411848.slurm. The loop adds only
# these folders and symlinks of these names, also a symlink with a missing target. A file such as submit.slurm stays.
# The name has no trailing slash, so rm removes a symlink of this name but not the files of its target. The pattern is in a string, as in paths, so bash -n can parse the script.
legacy_slurm_jobfolder_pattern="+([0-9]).slurm"
for legacy_slurm_jobfolder in $legacy_slurm_jobfolder_pattern; do
  if [[ -L "$legacy_slurm_jobfolder" || -d "$legacy_slurm_jobfolder" ]]; then
    paths+=" $legacy_slurm_jobfolder"
  fi
done

if [ 0 -lt $(ls -d -- $paths 2>/dev/null | wc -w) ]; then
  echo "The following ARTIS run files will be deleted:"
  ls -d -- $paths 2>/dev/null

  if [[ "$1" == "-y" ]]; then
    REPLY="y"
  else
    read -p "Are you sure you want to delete these ARTIS run files? [y/n]" -n 1 -r
    echo    # new line
  fi

  if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "Deleting:"
    rm -v -rf -- $paths 2>/dev/null
    if [[ -f "input-newrun.txt" ]]; then
      echo "cp input-newrun.txt input.txt"
      cp input-newrun.txt input.txt
    fi
  fi
else
  echo "No ARTIS run files to delete"
  if [[ -f "input-newrun.txt" ]] && ! cmp -s input-newrun.txt input.txt; then
    if [[ "$1" == "-y" ]]; then
      cp -v input-newrun.txt input.txt
    else
      cp -i -v input-newrun.txt input.txt
    fi
  elif [[ -f "input-newrun.txt" ]]; then
    echo "input-newrun.txt is the same as input.txt"
  fi
fi
