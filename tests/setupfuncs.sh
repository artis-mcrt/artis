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
