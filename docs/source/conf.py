# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
import os
import sys
from unittest.mock import MagicMock

# Mock C++ extension modules for documentation building
class BetterMock(MagicMock):
    # Add a proper __all__ attribute
    __all__ = ['ModelParams', 'Setups', 'ObsData', 'VegasMC']
    
    # Add documentation strings
    @classmethod
    def __getattr__(cls, name):
        mock = MagicMock()
        # Add docstrings to make autodoc happy
        mock.__doc__ = f"Mocked {name} class for documentation."
        # Make the signature inspection work better
        if name == '__init__':
            mock.__signature__ = None
        return mock

# Create fake module structure
systems_mock = type('VegasAfterglowC', (), {
    '__all__': ['ModelParams', 'Setups', 'ObsData', 'VegasMC'],
    'ModelParams': type('ModelParams', (), {'__doc__': 'ModelParams class documentation.'}),
    'Setups': type('Setups', (), {'__doc__': 'Setups class documentation.'}),
    'ObsData': type('ObsData', (), {'__doc__': 'ObsData class documentation.'}),
    'VegasMC': type('VegasMC', (), {'__doc__': 'VegasMC class documentation.'})
})

sys.modules['VegasAfterglow.VegasAfterglowC'] = systems_mock
sys.path.insert(0, os.path.abspath('../../'))

project = 'VegasAfterglow'
copyright = '2024, VegasAfterglow Team'
author = 'VegasAfterglow Team'
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
    'sphinx.ext.doctest',
    'sphinx.ext.todo',
    'sphinx.ext.ifconfig',
    'sphinx.ext.extlinks',  # For easily linking to external sites
    'breathe',
    'sphinx_rtd_theme',
]

templates_path = ['_templates']
exclude_patterns = []

# Enable todos
todo_include_todos = True

# Default role for inline markup
default_role = 'any'

# Define common links
extlinks = {
    'doxygen': ('doxygen/%s', '%s'),
    'source': ('doxygen/files.html#%s', '%s')
}

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_logo = '../../assets/logo.svg'
html_favicon = '../../assets/logo.svg'
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

# Add javascript for source code toggling
html_js_files = [
    'js/custom.js',
]
#html_add_permalinks = None    # older Sphinx
#html_permalinks     = "" 
# Add syntax highlighting style
pygments_style = 'sphinx'

# GitHub Pages settings
html_baseurl = 'https://yihanwangastro.github.io/VegasAfterglow/docs/'

# Create a custom javascript file for source code toggling
import os
js_dir = os.path.join(os.path.dirname(__file__), '_static', 'js')
os.makedirs(js_dir, exist_ok=True)
with open(os.path.join(js_dir, 'custom.js'), 'w') as f:
    f.write("""
// Function to toggle source code visibility
document.addEventListener('DOMContentLoaded', function() {
    // Add toggle buttons for implementation sections
    const sections = document.querySelectorAll('.breathe-sectiondef');
    sections.forEach(function(section) {
        if (section.querySelector('.cpp-source')) {
            const btn = document.createElement('button');
            btn.textContent = 'Toggle Implementation';
            btn.className = 'toggle-impl-btn';
            btn.style.cssText = 'background: #2980b9; color: white; border: none; padding: 5px 10px; margin: 5px 0; cursor: pointer; border-radius: 3px;';
            btn.onclick = function() {
                const sources = section.querySelectorAll('.cpp-source');
                sources.forEach(function(src) {
                    src.style.display = src.style.display === 'none' ? 'block' : 'none';
                });
            };
            section.insertBefore(btn, section.firstChild);
        }
    });
    
    // Add a link to source browser in the navigation
    const nav = document.querySelector('.wy-nav-side .wy-menu-vertical');
    if (nav) {
        const sourceLi = document.createElement('li');
        sourceLi.className = 'toctree-l1';
        const sourceLink = document.createElement('a');
        sourceLink.href = '/source_browser.html';
        sourceLink.textContent = 'Source Code Browser';
        sourceLi.appendChild(sourceLink);
        nav.appendChild(sourceLi);
    }
});
""")

# -- Breathe configuration ---------------------------------------------------
breathe_projects = {
    "VegasAfterglow": "../doxygen/xml"
}
# -- Breathe project and member defaults ------------------------------------
breathe_default_project               = "VegasAfterglow"                      # valid config key  [oai_citation:0‡Breathe](https://breathe.readthedocs.io/en/latest/directives.html?utm_source=chatgpt.com)
breathe_default_members               = ('members', 'undoc-members')           # valid config key  [oai_citation:1‡Breathe](https://breathe.readthedocs.io/en/latest/directives.html?utm_source=chatgpt.com)

# -- What to show in the documentation --------------------------------------
breathe_show_include                  = True                                  # valid config key  [oai_citation:2‡Breathe](https://breathe.readthedocs.io/en/latest/directives.html?utm_source=chatgpt.com)
breathe_show_enumvalue_initializer    = True                                  # valid config key  [oai_citation:3‡Breathe](https://breathe.readthedocs.io/en/latest/directives.html?utm_source=chatgpt.com)
breathe_show_define_initializer       = True                                  # valid config key  [oai_citation:4‡Breathe](https://breathe.readthedocs.io/en/latest/directives.html?utm_source=chatgpt.com)

# -- Parameter list placement (replaces breathe_separate_parameterlist) ------
breathe_order_parameters_first        = False                                  # valid config key  [oai_citation:5‡GitHub](https://github.com/breathe-doc/breathe/blob/main/breathe/renderer/sphinxrenderer.py?utm_source=chatgpt.com)

# -- Project-wide linking and file handling ---------------------------------
breathe_use_project_refids            = True                                  # valid config key  [oai_citation:6‡Breathe](https://breathe.readthedocs.io/en/latest/directives.html?utm_source=chatgpt.com)
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
breathe_show_include_files = True
breathe_doxygen_mapping = {
    # Map implementation file elements to their header declarations
    'function': 'function',
    'define': 'define',
    'property': 'variable',
    'variable': 'variable',
    'enum': 'enum',
    'enumvalue': 'enumvalue',
    'method': 'method',
    'typedef': 'typedef',
    'class': 'class',
    'struct': 'struct',
}

# Improved debug options for troubleshooting
breathe_debug_trace_directives = True
breathe_debug_trace_doxygen_ids = True
breathe_debug_trace_qualification = True

# Enhanced options for template and inline function documentation
breathe_template_relations = True
breathe_inline_details = True
breathe_show_define_initializer = True
breathe_show_template_parameters = True
breathe_show_templateparams = True

# -- intersphinx configuration -----------------------------------------------
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
}

# -- autodoc configuration ---------------------------------------------------
autodoc_member_order = 'groupwise'  # Changed to groupwise for better organization
autodoc_typehints = 'both'  # Changed to 'both' to show in signature and description
autodoc_default_options = {
    'members': True,
    'member-order': 'groupwise',
    'undoc-members': True,
    'private-members': True,  # Show private members
    'special-members': True,  # Show special members
    'show-inheritance': True,
    'inherited-members': True
}

# -- napoleon configuration --------------------------------------------------
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = True
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = True
napoleon_use_ivar = True
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = True
napoleon_attr_annotations = True

# -- Source code link configuration ------------------------------------------
viewcode_follow_imported_members = True
viewcode_enable_epub = True

# Custom implementation to include source from cpp files
def setup(app):
    # Custom directive for including source code from implementation files
    from docutils.parsers.rst import Directive
    from docutils import nodes
    
    class ImplementationDirective(Directive):
        required_arguments = 1  # function name
        optional_arguments = 1  # filename
        has_content = False
        
        def run(self):
            function_name = self.arguments[0]
            filename = self.arguments[1] if len(self.arguments) > 1 else None
            
            para = nodes.paragraph()
            para += nodes.Text(f"Implementation details for {function_name} can be found in the source browser.")
            
            return [para]
    
    app.add_directive('implementation', ImplementationDirective)
    
    # Add custom CSS class for implementation details
    app.add_css_file('css/custom.css') 