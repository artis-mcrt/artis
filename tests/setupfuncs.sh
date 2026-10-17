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

# Set the option $1 of artisoptions.h to $2. For a value option, $2 is the new value, e.g. 'true' or '1e5', and the
# option keeps its C++ type. For a function option, $2 is the new body with its braces, on one line or on more lines.
# The option name finds the statement. The statement must start with "constexpr" and have the name on its first line.
# Then it must have "= value;", only ";", or a parameter list and a body in braces. Stop the script unless:
# - $1 and $2 are given and $2 is not empty;
# - perl can read and change artisoptions.h;
# - the name matches exactly one statement;
# - $2 is a body in braces for a function option.
setoption() {
    if [ $# -ne 2 ] || [ -z "$2" ]; then
        echo "[error] setoption needs an option name and a value" >&2
        exit 1
    fi
    # perl -i only warns about a file that it cannot read or change, so END gives the error status for that case
    if ! OPTION_NAME="$1" OPTION_VALUE="$2" perl -0777 -i -pe '
        $file_read = 1;
        my $name = quotemeta $ENV{OPTION_NAME};
        my $value = $ENV{OPTION_VALUE};
        my $is_function;
        my $count = s/^(constexpr\b[^;{=\n]*?\b$name\b)(?:(\s*\([^)]*\)\s*)(\{(?:[^{}]++|(?3))*\})|\s*(?:=[^;]*)?;)/($is_function = defined $2) ? "$1$2$value" : "$1 = $value;"/gme || 0;
        die "[error] the option $ENV{OPTION_NAME} matches $count statements of artisoptions.h, not one\n" unless $count == 1;
        die "[error] the value of the function option $ENV{OPTION_NAME} must be a body in braces\n"
            if $is_function && $value !~ /^\s*\{.*\}\s*\z/s;
        END { $? = 1 unless $file_read; }
    ' artisoptions.h; then
        echo "[error] setoption could not set $1 in artisoptions.h" >&2
        exit 1
    fi
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
