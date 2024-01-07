# Configuration file for the Sphinx documentation builder.

# -- Project information
from datetime import datetime

project = "SCSilicon2"
author = "Xikang Feng"
repository_url = "https://github.com/xikanfeng2/SCSilicon2"
copyright = f"{datetime.now():%Y}, Xikang Feng."
version = "1.0.1"
release = "1.0.1"

# -- General configuration
extensions = [
    # read Markdown files
    "myst_nb",
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
    "sphinx.ext.linkcode",
    "sphinx_design",
    "sphinxext.opengraph",
]

html_theme = "sphinx_book_theme"
html_title = "SCSilicon2"
html_theme_options = {
    "home_page_in_toc": True,
    "github_url": "https://github.com/executablebooks/rst-to-myst",
    "repository_url": "https://github.com/executablebooks/rst-to-myst",
    "use_issues_button": True,
    "use_repository_button": True,
    "repository_branch": "main",
    "path_to_docs": "docs",
}

intersphinx_mapping = {
    "python": ("https://docs.python.org/3.8", None),
    "sphinx": ("https://www.sphinx-doc.org/en/master", None),
    "markdown_it": ("https://markdown-it-py.readthedocs.io/en/latest", None),
}

nitpick_ignore = [
    ("py:class", name) for name in ["IO", "_io.StringIO", "docutils.nodes.document"]
]
