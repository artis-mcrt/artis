#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# ///
"""Write docs/options.md from artisoptions_doc.md and the artisoptions_*.h presets.

artisoptions_doc.md holds one fenced C++ block. The block has a group marker of the form "// --- Name ---", and then
for each option a comment block and a declaration. This script turns each group into a section of the page, and each
option into a subsection with its type, its description, and the value of each preset.

The script also checks that the presets and artisoptions_doc.md name the same options. A missing or an extra option
stops the script, because the page would then describe options that no preset sets.

Without arguments, the script writes the page. It gives the status 1 when the content of the page changes, so that a
pre-commit hook adds the new content to the commit. With --check, the script writes nothing and only compares.
"""

import argparse
import re
import sys
import typing as t
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DOC_PATH = REPO_ROOT / "artisoptions_doc.md"
PAGE_PATH = REPO_ROOT / "docs" / "options.md"

GROUP_PATTERN = re.compile(r"^// --- (.+?) ---$")
COMMENT_PATTERN = re.compile(r"^//(?: (.*))?$")
DECLARATION_PATTERN = re.compile(r"^constexpr\s+(.+?)\s*\b([A-Z][A-Z0-9_]*)\b\s*(\([^)]*\))?\s*;$")

PAGE_HEADER = """# Compile-time options

`artisoptions.h` holds the compile-time options of a run. It is a symlink to one of the `artisoptions_*.h` presets.
[Option presets](presets.md) says what each preset is for, and [Installation](installation.md) says how to select one.

A change of an option needs a new build. The log file of rank 0 holds the options of the run, so the log of a finished
simulation shows the values that gave its results.

This page lists every option of the presets, with the value of each preset. `std::optional` options with no value are
empty, e.g. `GRID_TYPE_OVERRIDE` keeps the grid type of `model.txt`.
"""


class Option(t.NamedTuple):
    """One compile-time option of artisoptions.h."""

    name: str
    type_name: str
    parameters: str
    description: list[str]


class Group(t.NamedTuple):
    """A group of options, as a group marker of artisoptions_doc.md names it."""

    name: str
    options: list[Option]


def comment_to_markdown(comment_lines: list[str]) -> list[str]:
    """Turn the comment lines of an option into markdown blocks.

    A line that starts with "- " starts a list item, an indented line continues the current item or paragraph, and an
    empty line starts a new paragraph.
    """
    blocks: list[str] = []
    for line in comment_lines:
        stripped = line.strip()
        if not stripped:
            blocks.append("")
        elif stripped.startswith("- ") or not blocks or not blocks[-1]:
            blocks.append(stripped)
        else:
            blocks[-1] = f"{blocks[-1]} {stripped}"
    return [block for block in blocks if block]


def parse_doc(text: str) -> list[Group]:
    """Read the groups and the options of artisoptions_doc.md."""
    groups: list[Group] = []
    comment_lines: list[str] = []
    for linenum, line in enumerate(text.splitlines(), start=1):
        if line.startswith("```"):
            continue
        group_match = GROUP_PATTERN.match(line)
        comment_match = COMMENT_PATTERN.match(line)
        declaration_match = DECLARATION_PATTERN.match(line)
        if group_match:
            groups.append(Group(group_match.group(1), []))
            comment_lines = []
        elif declaration_match:
            if not groups:
                message = f"{DOC_PATH.name} line {linenum}: an option comes before the first group marker"
                raise ValueError(message)
            groups[-1].options.append(
                Option(
                    name=declaration_match.group(2),
                    type_name=declaration_match.group(1),
                    parameters=declaration_match.group(3) or "",
                    description=comment_to_markdown(comment_lines),
                )
            )
        elif comment_match:
            comment_lines.append(comment_match.group(1) or "")
        else:
            comment_lines = []
    return groups


def find_preset_value(preset_text: str, name: str) -> str | None:
    """Give the value that a preset sets for the option, "" for a declaration with no value, or None when absent.

    The value of a function option is its body with the braces, and it keeps its line breaks. A body that holds a
    comment needs them, because one line would put the code after the "//" of the comment.
    """
    match = re.search(rf"^constexpr\b[^;{{=\n]*?\b{re.escape(name)}\b", preset_text, re.MULTILINE)
    if match is None:
        return None
    position = match.end()
    while preset_text[position].isspace():
        position += 1
    if preset_text[position] == ";":
        return ""
    if preset_text[position] == "(":
        position = preset_text.index(")", position) + 1
        while preset_text[position].isspace():
            position += 1
        if preset_text[position] != "{":
            message = f"the function option {name} has no body in braces"
            raise ValueError(message)
        depth = 0
        for end, character in enumerate(preset_text[position:], start=position):
            depth += {"{": 1, "}": -1}.get(character, 0)
            if depth == 0:
                return preset_text[position : end + 1]
        message = f"the body of the function option {name} has no closing brace"
        raise ValueError(message)
    if preset_text[position] == "=":
        end = preset_text.index(";", position)
        return " ".join(preset_text[position + 1 : end].split())
    message = f"the declaration of {name} has no value, no body, and no semicolon"
    raise ValueError(message)


def preset_values(presets: dict[str, str], option: Option) -> list[str]:
    """Give the markdown blocks that hold the value of the option in each preset."""
    values = {preset: find_preset_value(text, option.name) for preset, text in presets.items()}
    missing = sorted(preset for preset, value in values.items() if value is None)
    if missing:
        message = f"the presets {', '.join(missing)} declare no option {option.name}"
        raise ValueError(message)
    shown = {preset: value or "no value" for preset, value in values.items()}
    one_line = all("\n" not in value for value in shown.values())
    if len(set(shown.values())) == 1 and one_line:
        return [f"**Every preset:** `{next(iter(shown.values()))}`"]
    if one_line and all(len(value) <= 40 for value in shown.values()):
        return ["**Preset values:** " + ", ".join(f"{preset} `{value}`" for preset, value in shown.items())]
    blocks: list[str] = ["**Preset values:**"]
    blocks.extend(f"{preset}\n\n```c++\n{value}\n```" for preset, value in shown.items())
    return blocks


def make_page(groups: list[Group], presets: dict[str, str]) -> str:
    """Give the content of docs/options.md."""
    blocks: list[str] = [PAGE_HEADER.strip()]
    for group in groups:
        blocks.append(f"## {group.name}")
        for option in group.options:
            blocks.extend((f"### {option.name}", f"`{option.type_name} {option.name}{option.parameters}`"))
            blocks.extend(option.description)
            blocks.extend(preset_values(presets, option))
    return "\n\n".join(blocks) + "\n"


def main() -> int:
    """Write docs/options.md, or compare it with --check."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="compare the page and write nothing")
    args = parser.parse_args()

    presets = {
        path.stem.removeprefix("artisoptions_"): path.read_text() for path in sorted(REPO_ROOT.glob("artisoptions_*.h"))
    }
    if not presets:
        print(f"[error] {REPO_ROOT} holds no artisoptions_*.h preset", file=sys.stderr)
        return 1
    groups = parse_doc(DOC_PATH.read_text())
    documented = {option.name for group in groups for option in group.options}
    for preset, text in presets.items():
        declared = set(re.findall(r"^constexpr\s+.*?\b([A-Z][A-Z0-9_]*)\b\s*[(;=]", text, re.MULTILINE))
        undocumented = sorted(declared - documented)
        if undocumented:
            print(f"[error] {DOC_PATH.name} describes no option {', '.join(undocumented)} of {preset}", file=sys.stderr)
            return 1

    page = make_page(groups, presets)
    current = PAGE_PATH.read_text() if PAGE_PATH.exists() else ""
    if page == current:
        return 0
    if args.check:
        print(f"[error] {PAGE_PATH.name} is out of date. Run scripts/make_options_docs.py", file=sys.stderr)
        return 1
    PAGE_PATH.write_text(page)
    print(f"{PAGE_PATH.name}: {len(documented)} options of {len(presets)} presets")
    return 1


if __name__ == "__main__":
    sys.exit(main())
