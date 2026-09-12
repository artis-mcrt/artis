#!/usr/bin/env bash

# only change the run folder if sn3d finished cleanly after the last timestep of the model. A run that
# stopped at an earlier timestep_finish needs its restart files to continue.
if grep -qs "No need for restart" output_0-0.txt && grep -qs "sn3d finished" output_0-0.txt; then
  rm -f packets_*.tmp gridsave_*.tmp vspecpol_*.tmp vpkt_grid_*.tmp

  # join 3D direction files, if they exist
  python3 ./artis/scripts/mergeangleres.py

  mkdir -p packets
  mv packets*.out* packets/ || true

  mkdir -p vpackets
  mv vpackets*.out* vpackets/ || true

  mkdir -p vspecpol
  mv vspecpol*.out* vspecpol/ || true

  mkdir -p vpkt_grid
  mv vpkt_grid*.out* vpkt_grid/ || true

  # remove empty directories, but keep the artis folder
  find . -maxdepth 1 -type d -empty ! -name artis -delete

  # One zstd process compresses its files one after the other. At level 13, -T0 divides a file
  # into jobs of 16 MiB, so a packet file of 44 MiB uses a maximum of three cores. xargs therefore
  # starts one zstd for each file, and runs one zstd on each CPU of the job. Without -n1, xargs
  # gives all the file names to one zstd. -T1 keeps the memory of each zstd small and gives the
  # same output as -T0.
  ncpus=$(nproc)

  # 3D kilonova model.txt and abundances.txt can be huge, so compress txt files
  # do maxdepth 1 first in case job gets killed during run folder compression
  echo "$(date): zstd compresses the .txt files of the run folder"
  find . -maxdepth 1 -name '*.txt' ! -name "output_0-0.txt" -size +200k -print0 | sort -z | xargs -r0 -n1 -P"$ncpus" zstd -T1 -13 -v --rm -f
  echo "$(date): zstd compresses the .out files of the run folder"
  find . -maxdepth 1 -name '*.out' ! -name "slurm-*.out" -size +200k -print0 | sort -z | xargs -r0 -n1 -P"$ncpus" zstd -T1 -13 -v --rm -f

  echo "$(date): zstd compresses the packet files"
  find packets/ -name 'packets*.out' -size +200k -print0 | sort -z | xargs -r0 -n1 -P"$ncpus" zstd -T1 -13 -v --rm -f

  # the artis folder holds the code and the data files of the code, so do not change them
  echo "$(date): zstd compresses the .txt files of the subfolders"
  find . -path ./artis -prune -o -name '*.txt' ! -name "output_0-0.txt" -size +200k -print0 | sort -z | xargs -r0 -n1 -P"$ncpus" zstd -T1 -13 -v --rm -f
  echo "$(date): zstd compresses the .out files of the subfolders"
  find . -path ./artis -prune -o -name '*.out' ! -name "slurm-*.out" -size +200k -print0 | sort -z | xargs -r0 -n1 -P"$ncpus" zstd -T1 -13 -v --rm -f

  echo "$(date): tar_rm_logs.sh archives the log files"
  ./artis/scripts/tar_rm_logs.sh

  echo "$(date): artistools converts the output files to parquet"
  export PATH="$(pwd)/../uv/bin:$PATH"
  export PATH="$HOME/.local/bin/:$PATH"
  if ! command -v uv >/dev/null 2>&1
  then
    # curl -LsSf https://astral.sh/uv/install.sh | sh
    # cosma disallows curl, so install uv with pip
    python3 -m ensurepip --upgrade
    python3 -m pip install --upgrade uv --target "$(pwd)/../uv"
  fi
  uv tool install -U artistools@latest

  # convert packets to parquet for fast reading
  uvx artistools lc --frompackets || true

  if [ -f vpkt.txt ]; then
    # convert virtual packets to parquet
    uvx artistools lc --frompackets -plotvspecpol 0 || true
  fi

  # convert estimators to parquet. On JUWELS, you might need to limit the number of processes to 16 in artistools/artistools/configuration.py
  uvx --from artistools -- python3 -c 'import artistools as at; at.estimators.scan_estimators()' || true

else
  echo "sn3d did not finish cleanly after the last timestep of the model, so exspec-after.sh changes no file"
fi
