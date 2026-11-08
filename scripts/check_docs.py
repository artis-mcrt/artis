#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# ///
"""Check that the hand-written tables of docs/ agree with the source.

The script makes three checks:

- each option of the help text of the Makefile is on docs/makeoptions.md;
- each log message that docs/troubleshooting.md quotes is in the source;
- the table of docs/inputfiles.md has the lines of inputlinecomments in input.cc, with the same unused lines.

The script gives the status 0 when each check passes, and 1 otherwise.
"""

import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DOCS_PATH = REPO_ROOT / "docs"

# a shorter part of a quoted message, e.g. "cell", also occurs by chance in the source
MIN_QUOTED_FRAGMENT_LENGTH = 12


def find_make_options_without_docs() -> list[str]:
    """Give a message for each option of the help text of the Makefile that docs/makeoptions.md does not name."""
    makefile = (REPO_ROOT / "Makefile").read_text(encoding="utf-8")
    help_match = re.search(r"^define HELP_TEXT$(.*?)^endef$", makefile, re.MULTILINE | re.DOTALL)
    if help_match is None:
        return ["Makefile: the definition of HELP_TEXT is absent"]
    options_text = help_match.group(1).partition("Options:")[2]
    option_names = set(re.findall(r"\b([A-Z_]+)=", options_text))
    page = (DOCS_PATH / "makeoptions.md").read_text(encoding="utf-8")
    return [f"docs/makeoptions.md does not name {name}" for name in sorted(option_names) if f"`{name}=" not in page]


def read_source_with_joined_literals() -> str:
    """Give the text of the C++ source, with adjacent string literals joined."""
    source = "\n".join(
        path.read_text(encoding="utf-8") for pattern in ("*.cc", "*.h") for path in sorted(REPO_ROOT.glob(pattern))
    )
    return re.sub(r'"\s*\n\s*"', "", source)


def find_quoted_messages_not_in_source(source: str) -> list[str]:
    """Give a message for each part of a log message of docs/troubleshooting.md that the source does not hold."""
    page = (DOCS_PATH / "troubleshooting.md").read_text(encoding="utf-8")
    quoted_section = re.sub(r"\{term\}`[^`]*`", "", page[page.index("### Warnings that the tests also give") :])
    code_spans = (" ".join(code_span.split()) for code_span in re.findall(r"`([^`]*)`", quoted_section))
    errors: list[str] = []
    # a quoted log message has more than one word. A command also has more than one word, and it starts with artistools.
    for quoted_message in (span for span in code_spans if " " in span and not span.startswith("artistools ")):
        for fragment in quoted_message.split("..."):
            stripped_fragment = fragment.strip(" :;.")
            if len(stripped_fragment) >= MIN_QUOTED_FRAGMENT_LENGTH and stripped_fragment not in source:
                errors.append(f"docs/troubleshooting.md quotes text that the source does not hold: {stripped_fragment}")
    return errors


def find_input_line_differences(source: str) -> list[str]:
    """Give a message for each difference between inputlinecomments in input.cc and the table of docs/inputfiles.md."""
    array_match = re.search(r"inputlinecomments = std::array\{(.*?)\n\};", source, re.DOTALL)
    if array_match is None:
        return ["input.cc: the array inputlinecomments is absent"]
    source_comments = re.findall(r'"\s*(\d+): (.*?)",?\n', array_match.group(1) + "\n")
    page = (DOCS_PATH / "inputfiles.md").read_text(encoding="utf-8")
    input_section = page[page.index("## input.txt") : page.index("## model.txt")]
    table_rows = re.findall(r"^\| (\d+) \| (.*?) ?\| .* \|$", input_section, re.MULTILINE)
    if len(source_comments) != len(table_rows):
        return [f"input.cc has {len(source_comments)} input lines, and docs/inputfiles.md has {len(table_rows)}"]
    errors: list[str] = []
    for (source_line, comment), (page_line, value_names) in zip(source_comments, table_rows, strict=True):
        if source_line != page_line:
            errors.append(f"input line {source_line} of input.cc has the number {page_line} in docs/inputfiles.md")
        if comment.startswith("UNUSED") != (not value_names.strip()):
            errors.append(f"input line {source_line}: input.cc and docs/inputfiles.md differ on whether it is unused")
        if not comment.startswith("UNUSED"):
            names_in_source = set(comment.partition(":")[0].split())
            names_on_page = set(value_names.strip("` ").split())
            errors.extend(
                f"input line {source_line}: docs/inputfiles.md and input.cc differ on {name}"
                for name in sorted(names_in_source ^ names_on_page)
            )
    return errors


def main() -> int:
    """Run the checks, and give the status 1 for a difference."""
    source = read_source_with_joined_literals()
    errors = [
        *find_make_options_without_docs(),
        *find_quoted_messages_not_in_source(source),
        *find_input_line_differences(source),
    ]
    for error in errors:
        print(f"[error] {error}", file=sys.stderr)
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
