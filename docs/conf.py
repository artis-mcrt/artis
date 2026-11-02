"""Sphinx configuration of the ARTIS documentation.

The pages are MyST Markdown. Read the Docs builds them with the settings of .readthedocs.yaml.
"""

import os

project = "ARTIS"
author = "The ARTIS authors"
project_copyright = "The ARTIS authors. BSD 3-Clause license"

extensions = ["myst_parser", "sphinx_copybutton"]

# colon_fence gives a directive in ::: fences, which a Markdown editor reads more easily than a backtick fence.
myst_enable_extensions = ["colon_fence"]
# make an anchor for each heading down to level 3, so that a link can point at one option of options.md
myst_heading_anchors = 3

templates_path = ["_templates"]
html_static_path = ["_static"]
html_css_files = ["version-selector.css"]
html_js_files = ["version-selector.js"]

html_theme = "furo"
html_title = "ARTIS"
# the logo is the tracked image of data/, because a zip archive of the repository holds no symlink
html_logo = "../data/artislogo.png"
html_favicon = "../data/artislogo.png"
html_copy_source = False
html_show_sourcelink = False

# the version selector comes under the name of the project. The other fragments are the default sidebar of Furo.
html_sidebars = {
    "**": [
        "sidebar/brand.html",
        "sidebar/version-selector.html",
        "sidebar/search.html",
        "sidebar/scroll-start.html",
        "sidebar/navigation.html",
        "sidebar/ethical-ads.html",
        "sidebar/scroll-end.html",
        "sidebar/variant-selector.html",
    ]
}

# the "Edit this page" link of each page points at the branch or the tag that Read the Docs builds. A build of a pull
# request gives the number of that request, which names no branch, so the link then points at the main branch.
source_branch = os.environ.get("READTHEDOCS_GIT_IDENTIFIER", "main")
if os.environ.get("READTHEDOCS_VERSION_TYPE") == "external":
    source_branch = "main"

# Furo takes these three values for the "Edit this page" link of each page
html_theme_options = {
    "source_repository": "https://github.com/artis-mcrt/artis/",
    "source_branch": source_branch,
    "source_directory": "docs/",
}
