# docs/conf.py
import os
import sys
from pathlib import Path

# Make your Rust extension importable
sys.path.insert(0, os.path.abspath(".."))

project = 'bioshell'
copyright = '2006 - 2026, Dominik Gront'
author = 'Dominik Gront'

extensions = [
    "myst_parser",           # ← Enables Markdown
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",   # Support for Google/NumPy docstrings
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx_autodoc_typehints",  # Type hints support
]

exclude_patterns = ['_build', '.DS_Store']

# MyST settings
myst_enable_extensions = [
    "colon_fence",           # ::: directive support (like mkdocstrings)
    "deflist",
    "html_image",
]

# TOC settings
toc_object_entries = False

# Autodoc settings
autodoc_member_order = 'bysource'
napoleon_google_docstring = True
napoleon_numpy_docstring = True

html_theme = "sphinx_rtd_theme"   # or "furo", "sphinx_book_theme", etc.