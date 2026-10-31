"""Sphinx configuration of the ARTIS documentation.

The pages are MyST Markdown. Read the Docs builds them with the settings of .readthedocs.yaml.
"""

project = "ARTIS"
author = "The ARTIS authors"
project_copyright = "The ARTIS authors. BSD 3-Clause license"

extensions = ["myst_parser", "sphinx_copybutton"]

# colon_fence gives a directive in ::: fences, which a Markdown editor reads more easily than a backtick fence.
myst_enable_extensions = ["colon_fence"]
# make an anchor for each heading down to level 3, so that a link can point at one option of options.md
myst_heading_anchors = 3

html_theme = "furo"
html_title = "ARTIS"
# the logo is the tracked image of data/, because a zip archive of the repository holds no symlink
html_logo = "../data/artislogo.png"
html_favicon = "../data/artislogo.png"
html_copy_source = False
html_show_sourcelink = False

# Furo takes these three values for the "Edit this page" link of each page
html_theme_options = {
    "source_repository": "https://github.com/artis-mcrt/artis/",
    "source_branch": "main",
    "source_directory": "docs/",
}
