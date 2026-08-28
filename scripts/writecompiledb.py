#!/usr/bin/env python3
"""Write the compilation database that clangd and clang-tidy read.

The Makefile calls this script. It supplies the compiler, the compile flags, and the
source files in the environment, because the Makefile already knows all of them.

Do not make this database with compiledb. compiledb does not recognise mpicxx as a
compiler, so it removes the compiler and the first flags from each command. The first
element of "arguments" then becomes a flag. clang-tidy takes that flag as the name of
the compiler, finds no resource directory, and finds no standard header. Each parse
error then gives many false diagnostics, e.g. cppcoreguidelines-pro-type-member-init on
a structure that has initialisers.
"""

import json
import os
import shlex
import sys

if len(sys.argv) < 3:
    sys.exit("usage: writecompiledb.py OUTFILE SOURCE...")

outfile = sys.argv[1]
sources = sys.argv[2:]

compiler = os.environ.get("CDB_CXX", "").strip()
if not compiler:
    sys.exit("writecompiledb.py: CDB_CXX is empty, so clang-tidy would find no header")

flags = shlex.split(os.environ.get("CDB_FLAGS", ""))
directory = os.getcwd()

entries = [
    {"directory": directory, "arguments": [compiler, *flags, "-c", source], "file": source}
    for source in sorted(set(sources))
]

with open(outfile, "w", encoding="utf-8") as file:
    json.dump(entries, file, indent=1)
    file.write("\n")

print(f"{outfile}: {len(entries)} entries for {compiler}")
