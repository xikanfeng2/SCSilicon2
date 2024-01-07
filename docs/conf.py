from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING

HERE = Path(__file__).parent
sys.path[:0] = [str(HERE.parent), str(HERE / "extensions")]



# -- General configuration ------------------------------------------------


nitpicky = True  # Warn about broken links. This is here for a reason: Do not change.
needs_sphinx = "4.0"  # Nicer param docs
suppress_warnings = [
    "ref.citation",
    "myst.header",  # https://github.com/executablebooks/MyST-Parser/issues/262
]

# General information
project = "SCSilicon2"
author = "Xikang Feng"
repository_url = "https://github.com/xikanfeng2/SCSilicon2"
copyright = f"{datetime.now():%Y}, Xikang Feng."
version = "1.0.1"



# default settings
templates_path = ["_templates"]
master_doc = "index"
default_role = "literal"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.doctest",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx.ext.extlinks",
    "sphinx.ext.linkcode",
]



typehints_defaults = "braces"

pygments_style = "default"
pygments_dark_style = "native"


# -- Options for HTML output ----------------------------------------------

html_theme = "sphinx_book_theme"
html_theme_options = {
    "repository_url": repository_url,
    "use_repository_button": True,
}
html_static_path = ["_static"]
html_css_files = ["css/override.css"]
html_show_sphinx = False
# html_logo = "_static/img/Scanpy_Logo_BrightFG.svg"
html_title = "SCSilicon2"