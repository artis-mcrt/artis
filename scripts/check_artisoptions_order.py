#!/usr/bin/env python3
"""Check that each artisoptions_*.h preset has the options and the headings of artisoptions_doc.md in the same order.

The script compares the sequence of the heading lines ("// --- ... ---") and of the option names. It prints the first
difference of each preset and returns 1 if a preset differs from artisoptions_doc.md.
"""

import re
import sys
from pathlib import Path

HEADING_PATTERN = re.compile(r"^// --- .* ---$")
OPTION_PATTERN = re.compile(r"^constexpr\s+.*?\b([A-Za-z_][A-Za-z_0-9]*)\s*(?:\(|\{|=|;)")


def get_headings_and_option_names(path: Path) -> list[str]:
    headings_and_option_names = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if HEADING_PATTERN.match(line):
            headings_and_option_names.append(line)
        elif match := OPTION_PATTERN.match(line):
            headings_and_option_names.append(match.group(1))
    return headings_and_option_names


def main() -> int:
    repo_path = Path(__file__).resolve().parent.parent
    doc_path = repo_path / "artisoptions_doc.md"
    doc_sequence = get_headings_and_option_names(doc_path)
    preset_paths = sorted(repo_path.glob("artisoptions_*.h"))
    presets_with_difference = 0
    for preset_path in preset_paths:
        preset_sequence = get_headings_and_option_names(preset_path)
        if preset_sequence == doc_sequence:
            continue
        presets_with_difference += 1
        index = next(
            (i for i, (a, b) in enumerate(zip(preset_sequence, doc_sequence, strict=False)) if a != b),
            min(len(preset_sequence), len(doc_sequence)),
        )
        preset_item = preset_sequence[index] if index < len(preset_sequence) else "(end of file)"
        doc_item = doc_sequence[index] if index < len(doc_sequence) else "(end of file)"
        print(
            f"[error] {preset_path.name} differs from {doc_path.name} at item {index + 1}: "
            f"the preset has '{preset_item}' and the doc file has '{doc_item}'"
        )
    if presets_with_difference:
        print("Give each preset the headings and the options of artisoptions_doc.md in the same order.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
