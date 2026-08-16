#!/usr/bin/env bash

# only compress the files if we successfully ran exspec
if [[ -f emission.out || -f emission.out.zst || -f emissionpol.out ]]; then
  # keep the restart files unless the simulation has completed all of its timesteps: a run stopped
  # at an earlier timestep_finish needs them to continue
  if grep -qs "No need for restart" output_0-0.txt; then
    rm -f packets_*.tmp gridsave_*.tmp vspecpol_*.tmp vpkt_grid_*.tmp
  else
    echo "The simulation has not completed all of its timesteps, so keeping the restart files"
  fi

  mkdir -p speclc_angle_res
  mv *_res_*.out* speclc_angle_res/ || true
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

  # remove empty directories
  find . -maxdepth 1 -type d -empty -delete

  # 3D kilonova model.txt and abundances.txt can be huge, so compress txt files
  # do maxdepth 1 first in case job gets killed during run folder compression
  find . -maxdepth 1 -name '*.txt' ! -name "output_0-0.txt" -size +200k -print0 | sort -z | xargs -r0 -P8 zstd -T0 -13 -v --rm -f
  find . -maxdepth 1 -name '*.out' ! -name "slurm-*.out" -size +200k -print0 | sort -z | xargs -r0 -P8 zstd -T0 -13 -v --rm -f

  find packets/ -name 'packets*.out' -size +200k -print0 | sort -z | xargs -r0 -P8 zstd -T0 -13 -v --rm -f

  find . -name '*.txt' ! -name "output_0-0.txt" -size +200k -print0 | sort -z | xargs -r0 -P8 zstd -T0 -13 -v --rm -f
  find . -name '*.out' ! -name "slurm-*.out" -size +200k -print0 | sort -z | xargs -r0 -P8 zstd -T0 -13 -v --rm -f

  ./artis/scripts/tar_rm_logs.sh

  export PATH="$(pwd)/../uv/bin:$PATH"
  export PATH="$HOME/.local/bin/:$PATH"
  if ! command -v uv >/dev/null 2>&1
  then
    # curl -LsSf https://astral.sh/uv/install.sh | sh
    # cosma disallows curl, so install uv with pip
    python3 -m ensurepip --upgrade
    python3 -m pip install --upgrade uv --target "$(pwd)/../uv"
  fi
  uv tool install -U --no-cache artistools@latest

  # convert packets to parquet for fast reading
  uvx artistools lc --frompackets || true

  if [ -f vpkt.txt ]; then
    # convert virtual packets to parquet
    uvx artistools lc --frompackets -plotvspecpol 0 || true
  fi

  # convert estimators to parquet. On JUWELS, you might need to limit the number of processes to 16 in artistools/artistools/configuration.py
  uvx --from artistools -- python3 -c 'import artistools as at; at.estimators.scan_estimators()' || true

fi
