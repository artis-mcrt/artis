#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# ///
"""Write docs/options.md from artisoptions_doc.md and the artisoptions_*.h presets.

artisoptions_doc.md holds one fenced C++ block. The block has a group marker of the form "// --- Name ---", and then
for each option a comment block and a declaration. This script turns each group into a section of the page, and each
option into a subsection with its type, its description, and the value of each preset. A comment block with no
declaration after it becomes a note at the top of its section.

The script reads the presets that the repository tracks, so an untracked preset of your own stays out of the page.
It also checks that the presets and artisoptions_doc.md name the same options, in both directions. An option that the
doc does not describe, and an option that a preset does not declare, both stop the script. Both checks cover every
tracked preset, also a preset that PRESETS_NOT_ON_PAGE keeps off the page.

The script gives the status 0 when it wrote the page or when the page was already current, and 1 for an error. The
pre-commit hook fails when the script changes the page, so that the commit holds the new content.
"""

import re
import subprocess  # ruff: ignore[suspicious-subprocess-import]
import sys
import typing as t
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DOC_PATH = REPO_ROOT / "artisoptions_doc.md"
PAGE_PATH = REPO_ROOT / "docs" / "options.md"

# the presets that the page leaves out, because they are work in progress and their values will change
PRESETS_NOT_ON_PAGE = ("kilonova_nlte",)

GROUP_PATTERN = re.compile(r"^// --- (.+?) ---$")
COMMENT_PATTERN = re.compile(r"^//(.*)$")
DECLARATION_PATTERN = re.compile(r"^constexpr\s+(.+?)\s*\b([A-Z][A-Z0-9_]*)\b\s*(\([^)]*\))?\s*;$")
# the name of an option of a preset, before its value, its body, or its semicolon
PRESET_OPTION_PATTERN = re.compile(r"^constexpr\s+.*?\b([A-Z][A-Z0-9_]*)\b\s*[(;={]", re.MULTILINE)

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
    notes: list[str]
    options: list[Option]


def escape_markdown(text: str) -> str:
    """Escape the characters of a comment that markdown reads as markup.

    Markdown reads "<rank>" as an HTML tag and shows nothing for it, so the text needs the backslash.
    """
    return text.replace("<", r"\<")


def comment_to_markdown(comment_lines: list[str]) -> list[str]:
    """Turn the comment lines of an option into markdown blocks.

    A line that starts with "- " starts a list item. An indented line continues the block above it, and so does a line
    that follows a paragraph. A line that follows a list item without an indent starts a new paragraph, and an empty
    line also starts one.
    """
    blocks: list[str] = []
    for line in comment_lines:
        stripped = escape_markdown(line.strip())
        previous = blocks[-1] if blocks else ""
        continues_block = bool(previous) and (line.startswith(" ") or not previous.startswith("- "))
        if not stripped:
            blocks.append("")
        elif stripped.startswith("- "):
            blocks.append(stripped)
        elif continues_block:
            blocks[-1] = f"{blocks[-1]} {stripped}"
        else:
            blocks.append(stripped)
    return [block for block in blocks if block]


def parse_doc(text: str) -> list[Group]:
    """Read the groups, the notes, and the options of artisoptions_doc.md."""
    groups: list[Group] = []
    comment_lines: list[str] = []
    comment_has_option = False
    for linenum, line in enumerate(text.splitlines(), start=1):
        if line.startswith("```"):
            continue
        group_match = GROUP_PATTERN.match(line)
        comment_match = COMMENT_PATTERN.match(line)
        declaration_match = DECLARATION_PATTERN.match(line)
        if declaration_match:
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
            comment_has_option = True
        elif comment_match and not group_match:
            text_of_line = comment_match.group(1)
            comment_lines.append(text_of_line.removeprefix(" "))
        else:
            # a group marker or a blank line ends the comment block. A block that no declaration used is a note
            if comment_lines and not comment_has_option and groups:
                groups[-1].notes.extend(comment_to_markdown(comment_lines))
            comment_lines = []
            comment_has_option = False
            if group_match:
                groups.append(Group(group_match.group(1), [], []))
    return groups


def skip_spaces_and_comments(text: str, position: int, name: str) -> int:
    """Give the index of the next character that is neither a space nor part of a // comment."""
    while position < len(text):
        if text[position].isspace():
            position += 1
        elif text.startswith("//", position):
            end_of_line = text.find("\n", position)
            position = len(text) if end_of_line < 0 else end_of_line
        else:
            return position
    message = f"the declaration of {name} ends without a value, a body, or a semicolon"
    raise ValueError(message)


def find_end_of_body(text: str, start: int, name: str) -> int:
    """Give the index after the closing brace of the body that starts at the brace at start."""
    depth = 0
    position = start
    while position < len(text):
        if text.startswith("//", position):
            end_of_line = text.find("\n", position)
            position = len(text) if end_of_line < 0 else end_of_line
            continue
        depth += {"{": 1, "}": -1}.get(text[position], 0)
        if depth == 0:
            return position + 1
        position += 1
    message = f"the body of the function option {name} has no closing brace"
    raise ValueError(message)


def find_end_of_value(text: str, start: int, name: str) -> int:
    """Give the index of the semicolon that ends the value which starts at start."""
    position = start
    while position < len(text):
        if text.startswith("//", position):
            end_of_line = text.find("\n", position)
            position = len(text) if end_of_line < 0 else end_of_line
            continue
        if text[position] == ";":
            return position
        position += 1
    message = f"the value of {name} has no semicolon"
    raise ValueError(message)


def find_preset_value(preset_text: str, name: str) -> str | None:
    """Give the value that a preset sets for the option, "" for a declaration with no value, or None when absent.

    The value of a function option is its body with the braces. A value keeps its line breaks when it holds a comment,
    because one line would put the code after the "//" of that comment.
    """
    match = re.search(rf"^constexpr\b[^;{{=\n]*?\b{re.escape(name)}\b", preset_text, re.MULTILINE)
    if match is None:
        return None
    position = skip_spaces_and_comments(preset_text, match.end(), name)
    if preset_text[position] == ";":
        return ""
    if preset_text[position] == "(":
        position = skip_spaces_and_comments(preset_text, preset_text.index(")", position) + 1, name)
        if preset_text[position] != "{":
            message = f"the function option {name} has no body in braces"
            raise ValueError(message)
        return preset_text[position : find_end_of_body(preset_text, position, name)]
    if preset_text[position] in "={":
        start = position + 1 if preset_text[position] == "=" else position
        value = preset_text[start : find_end_of_value(preset_text, start, name)]
        return value.strip("\n") if "//" in value else " ".join(value.split())
    message = f"the declaration of {name} has no value, no body, and no semicolon"
    raise ValueError(message)


def preset_values(presets: dict[str, str], option: Option) -> list[str]:
    """Give the markdown blocks that hold the value of the option in each preset of the page.

    Every tracked preset must declare the option, also a preset that PRESETS_NOT_ON_PAGE keeps off the page.
    """
    values = {preset: find_preset_value(text, option.name) for preset, text in presets.items()}
    missing = sorted(preset for preset, value in values.items() if value is None)
    if missing:
        message = f"the presets {', '.join(missing)} declare no option {option.name}"
        raise ValueError(message)
    shown = {preset: value or "no value" for preset, value in values.items() if preset not in PRESETS_NOT_ON_PAGE}
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
        blocks.extend(group.notes)
        for option in group.options:
            blocks.extend((f"### {option.name}", f"`{option.type_name} {option.name}{option.parameters}`"))
            blocks.extend(option.description)
            blocks.extend(preset_values(presets, option))
    return "\n\n".join(blocks) + "\n"


def read_presets() -> dict[str, str]:
    """Give the text of each preset that the repository tracks, by the name of the preset."""
    try:
        listing = subprocess.run(  # ruff: ignore[subprocess-without-shell-equals-true]
            ["git", "-C", str(REPO_ROOT), "ls-files", "artisoptions_*.h"],  # ruff: ignore[start-process-with-partial-path]
            capture_output=True,
            text=True,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        message = f"git cannot list the presets of {REPO_ROOT}: {exc}"
        raise ValueError(message) from exc
    presets = {
        Path(line).stem.removeprefix("artisoptions_"): (REPO_ROOT / line).read_text(encoding="utf-8")
        for line in sorted(listing.stdout.split())
    }
    if not presets:
        message = f"{REPO_ROOT} tracks no artisoptions_*.h preset"
        raise ValueError(message)
    return presets


def check_every_option_documented(presets: dict[str, str], documented: set[str]) -> None:
    """Stop when a preset holds an option that artisoptions_doc.md does not describe."""
    for preset, text in presets.items():
        undocumented = sorted(set(PRESET_OPTION_PATTERN.findall(text)) - documented)
        if undocumented:
            message = f"{DOC_PATH.name} describes no option {', '.join(undocumented)} of {preset}"
            raise ValueError(message)


def build_page() -> tuple[str, int, int]:
    """Give the content of the page, the number of options, and the number of presets."""
    presets = read_presets()
    groups = parse_doc(DOC_PATH.read_text(encoding="utf-8"))
    documented = {option.name for group in groups for option in group.options}
    check_every_option_documented(presets, documented)
    untracked = sorted(set(PRESETS_NOT_ON_PAGE) - presets.keys())
    if untracked:
        message = f"PRESETS_NOT_ON_PAGE names {', '.join(untracked)}, which the repository does not track"
        raise ValueError(message)
    return make_page(groups, presets), len(documented), len(presets) - len(PRESETS_NOT_ON_PAGE)


def main() -> int:
    """Write docs/options.md, and give the status 1 for an error."""
    try:
        page, noptions, npresets = build_page()
    except ValueError as exc:
        print(f"[error] {exc}", file=sys.stderr)
        return 1

    if page != (PAGE_PATH.read_text(encoding="utf-8") if PAGE_PATH.exists() else ""):
        PAGE_PATH.write_text(page, encoding="utf-8")
        print(f"{PAGE_PATH.name}: {noptions} options of {npresets} presets")
    return 0


if __name__ == "__main__":
    sys.exit(main())
