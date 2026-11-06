#!/usr/bin/env bash
# Make a model folder for a new simulation.
# The folder gets the ARTIS source in the subfolder artis, and an input.txt from scripts/input-template.txt.

set -euo pipefail

symlink_source=false
if [ "${1:-}" = "--symlink" ]; then
    symlink_source=true
    shift
fi

if [ $# -ne 1 ] || [ "${1#-}" != "$1" ]; then
    echo "Usage: $0 [--symlink] <model folder>" >&2
    echo "  The script clones the ARTIS repository into <model folder>/artis." >&2
    echo "  With --symlink, it makes a symlink to this copy of the source instead." >&2
    exit 1
fi

model_folder=$1
artis_source=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)

if [ -e "$model_folder/artis" ] || [ -L "$model_folder/artis" ]; then
    echo "$0: '$model_folder/artis' exists, so the script changes nothing." >&2
    exit 1
fi

mkdir -p "$model_folder"

model_folder_absolute=$(cd "$model_folder" && pwd -P)
if [ "$model_folder_absolute" = "$artis_source" ]; then
    echo "$0: '$model_folder' is the folder of the ARTIS source, so the script changes nothing." >&2
    exit 1
fi

if $symlink_source; then
    ln -s "$artis_source" "$model_folder/artis"
else
    # git runs in the source folder, so that a relative URL of origin gives the same repository as there
    git -C "$artis_source" clone "$(git -C "$artis_source" remote get-url origin)" "$model_folder_absolute/artis"
fi

if [ -e "$model_folder/input.txt" ] || [ -L "$model_folder/input.txt" ]; then
    echo "The script keeps the input.txt of '$model_folder'."
else
    cp "$model_folder/artis/scripts/input-template.txt" "$model_folder/input.txt"
fi

if [ -e "$model_folder/artis/artisoptions.h" ] || [ -L "$model_folder/artis/artisoptions.h" ]; then
    build_step="  1. The source already has an artisoptions.h. Build with these options:
       cd '$model_folder/artis' && make"
else
    build_step="  1. Select a preset and build, e.g. with the classic preset:
       cd '$model_folder/artis' && ln -s artisoptions_classic.h artisoptions.h && make"
fi

cat << EOF

The model folder '$model_folder' is ready. The next steps:
$build_step
  2. Put model.txt, abundances.txt, and an atomic database into '$model_folder'.
  3. Set lines 1 to 3, 18, and 19 of '$model_folder/input.txt' for your model.
     They hold the timesteps, the times, and the LTE and grey timesteps.
EOF
