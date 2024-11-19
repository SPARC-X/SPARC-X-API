# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "SPARC-X-API"
copyright = "2024, SPARC-X Developmers"
author = "Tian Tian, Lucas R Timmerman, Ben Comer"


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.githubpages",
    "sphinx.ext.coverage",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",  # For Google/NumPy style docstrings
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",  # Adds links to source code
    "myst_parser",
]

source_suffix = {
    ".md": "markdown",
    ".txt": "markdown",
    ".rst": "restructuredtext",
}

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_context = {
    "display_github": True,  # Integrate GitHub
    "github_user": "sparc-x",  # Username
    "github_repo": "SPARC-X-API",  # Repo name
    "github_version": "master",  # Version
    "conf_py_path": "/doc/",  # Path in the checkout to the docs root
}

myst_enable_extensions = [
    "html_admonition",
]

coverage_Show_missing_items = True
autosummary_generate = True
