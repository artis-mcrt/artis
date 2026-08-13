#!/usr/bin/env bash

# for the +([0-9]) patterns that only match generated filenames (e.g. output_0-0.txt but not output_0-0_backup.txt)
shopt -s extglob

paths="gridsave*.tmp packets*.tmp vspecpol*.tmp vpackets*.tmp *.out *.out.* out.txt output_+([0-9])-+([0-9]).txt?(.zst|.gz|.xz) exspec*.txt* machine.file.* core.* *.slurm packets vspecpol vpackets speclc_angle_res bflist.dat ratecoeff.dat line_list.txt logfiles.tar*"

if [ 0 -lt $(ls $paths 2>/dev/null | wc -w) ]; then
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
  else
    echo "input-newrun.txt is the same as input.txt"
  fi
fi
