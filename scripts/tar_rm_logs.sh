#!/usr/bin/env bash

# a failure of find or tar must stop the script before it removes the log files
set -o pipefail

# The archive gets its final name only when it is complete. The partial archive is in the model folder,
# because /tmp of a compute node can be too small for the log files. A killed run can leave a partial archive.
partial_archive="logfiles.tar.zst.partial"
rm -f "$partial_archive"

if compgen -G 'logfiles.tar*' > /dev/null; then
    echo "logfiles.tar* already exists! Not overwriting"
    exit 1
fi

trap 'rm -f "$partial_archive"' EXIT

# the artis folder holds the code and the data files of the code, so do not change them
find . -mindepth 2 ! -path './artis/*' -name "output_*.txt*" -print0 \
    | tar -cv --null --files-from - -f - \
    | zstd -T0 -15 -f -o "$partial_archive" \
    && mv -v "$partial_archive" logfiles.tar.zst \
    && find . -mindepth 2 ! -path './artis/*' -name "output_*.txt*" ! -name "output_0-0.txt*" -delete
