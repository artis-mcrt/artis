# Helpers for the setup scripts of the end-to-end tests. The setup scripts source this file.
# CI sources the setup scripts from bash, and a local run uses zsh, so the file must work in both.

# The release that supplies the atomic data archives. The cache keys of ci.yml hold
# the hash of this file, so a new tag gives a new cache.
ATOMICDATA_RELEASE=v2026.5.15

getatomicdata() {
    if [ ! -f "$1" ]; then
        curl -fL --retry 3 -o "$1.part" "https://github.com/artis-mcrt/artis/releases/download/${ATOMICDATA_RELEASE}/$1" && mv "$1.part" "$1"
    fi
}

# Replace each match of the pattern $1 in artisoptions.h with $2. Neither can contain a | character.
# sed gives no error for a pattern that matches nothing. The test then runs with the default
# value of the preset.
sedopt() {
    if ! grep -q -e "$1" artisoptions.h; then
        echo "[error] this pattern matches no line of artisoptions.h: $1" >&2
        exit 1
    fi
    sed -i.bak -e "s|$1|$2|g" artisoptions.h
}

# Set the maximum number of levels of each ion (column 5 of compositiondata.txt) to $1 for every element.
# The test then reads only the lowest $1 levels of each ion from the atomic data.
setnlevelsmax() {
    if ! awk -v n="$1" 'NR == 1 {nelements = $1} NR > 3 && NF == 7 {$5 = n; count++} {print} END {if (count != nelements) exit 1}' compositiondata.txt > compositiondata.txt.new; then
        echo "[error] compositiondata.txt does not have one line of seven columns for each element" >&2
        exit 1
    fi
    mv compositiondata.txt.new compositiondata.txt
}
