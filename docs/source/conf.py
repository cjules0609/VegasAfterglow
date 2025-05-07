# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

project = 'Afterglow'
copyright = '2024, Afterglow Team'
author = 'Afterglow Team'
release = '0.1.0'  # Update this with your actual version

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',  # Support for NumPy and Google style docstrings
    'sphinx.ext.viewcode',
    'sphinx.ext.autosummary',
    'sphinx.ext.graphviz',
    'breathe',
    'sphinx_rtd_theme',
]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_logo = '../../logo.svg'
html_theme_options = {
    'logo_only': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': True,
    'navigation_depth': 4,
    'collapse_navigation': False,
    'sticky_navigation': True,
    'includehidden': True,
}
html_css_files = [
    'css/custom.css',
]

# GitHub Pages settings
html_baseurl = 'https://yihanwangastro.github.io/Afterglow/docs/'
html_use_index = True
html_copy_source = False  # Don't copy rst files to the output

# -- Breathe configuration ---------------------------------------------------
breathe_projects = {
    "Afterglow": "../doxygen/xml/"
}
breathe_default_project = "Afterglow"
breathe_default_members = ('members', 'undoc-members', 'protected-members', 'private-members')
breathe_show_include = False
breathe_show_enumvalue_initializer = True
breathe_show_define_initializer = True
breathe_show_details = True
breathe_separate_parameterlist = True
breathe_use_project_refids = True
breathe_implementation_filename_extensions = ['.c', '.cc', '.cpp']

# Additional Breathe customization for better detailed documentation
breathe_domain_by_extension = {
    "h": "cpp",
    "hpp": "cpp",
    "c": "c",
    "cpp": "cpp",
    "cc": "cpp",
}
breathe_domain_by_file_pattern = {
    '*/include/*': 'cpp',
    '*/src/*': 'cpp',
}
breathe_ordered_classes = True
breathe_show_include_files = False

# Debug options to trace issues with Breathe and Doxygen
breathe_debug_trace_directives = True
breathe_debug_trace_doxygen_ids = True
breathe_debug_trace_qualification = True

# Enhanced options for template and inline function documentation
breathe_template_relations = True
breathe_inline_details = True
breathe_show_define_initializer = True
breathe_show_template_parameters = True

# -- intersphinx configuration -----------------------------------------------
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
}

# -- autodoc configuration ---------------------------------------------------
autodoc_member_order = 'bysource'
autodoc_typehints = 'description'
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
    'inherited-members': True
}

# -- napoleon configuration --------------------------------------------------
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = True
napoleon_include_special_with_doc = True 