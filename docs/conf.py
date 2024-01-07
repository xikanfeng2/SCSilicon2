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
    "sphinx_copybutton",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.doctest",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx.ext.extlinks",
    "matplotlib.sphinxext.plot_directive",
    "sphinx_autodoc_typehints",  # needs to be after napoleon
    "git_ref",  # needs to be before scanpydoc.rtd_github_links
    "sphinx.ext.linkcode",
    "sphinx_design",
    "sphinxext.opengraph",
    *[p.stem for p in (HERE / "extensions").glob("*.py") if p.stem not in {"git_ref"}],
]

# Generate the API documentation when building
autosummary_generate = True
autodoc_member_order = "bysource"
# autodoc_default_flags = ['members']
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
napoleon_custom_sections = [("Params", "Parameters")]
todo_include_todos = False
api_dir = HERE / "api"  # function_images
myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_image",
    "html_admonition",
]
myst_url_schemes = ("http", "https", "mailto")
nb_output_stderr = "remove"
nb_execution_mode = "off"
nb_merge_streams = True


ogp_site_url = "https://scsilicon2.readthedocs.io/en/stable/"

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


def setup(app: Sphinx):
    """App setup hook."""
    app.add_config_value(
        "recommonmark_config",
        {
            "auto_toc_tree_section": "Contents",
            "enable_auto_toc_tree": True,
            "enable_math": True,
            "enable_inline_math": False,
            "enable_eval_rst": True,
        },
        True,
    )


# -- Options for other output formats ------------------------------------------

htmlhelp_basename = f"{project}doc"
doc_title = f"{project} Documentation"
latex_documents = [(master_doc, f"{project}.tex", doc_title, author, "manual")]
man_pages = [(master_doc, project, doc_title, [author], 1)]
texinfo_documents = [
    (
        master_doc,
        project,
        doc_title,
        author,
        project,
        "One line description of project.",
        "Miscellaneous",
    )
]
