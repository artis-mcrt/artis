#!/usr/bin/env bash

if compgen -G 'logfiles.tar*' > /dev/null; then
    echo "logfiles.tar* already exists! Not overwriting"
    exit 1
fi

# The archive gets its final name only when it is complete. The partial archive is in the model folder,
# because /tmp of a compute node can be too small for the log files.
partial_archive="partial_logfiles.tar"
trap 'rm -f "$partial_archive"' EXIT

# exspec-after.sh compresses each large log file with zstd, so tar adds no second compression.
# the artis folder holds the code and the data files of the code, so do not change them
find . -mindepth 2 ! -path './artis/*' -name "output_*.txt*" -print0 > "$partial_archive.filelist" \
    && tar -cv --null --files-from "$partial_archive.filelist" -f "$partial_archive" \
    && mv -v "$partial_archive" logfiles.tar \
    && find . -mindepth 2 ! -path './artis/*' -name "output_*.txt*" ! -name "output_0-0.txt*" -delete
rm -f "$partial_archive.filelist"
