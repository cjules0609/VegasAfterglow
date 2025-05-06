# VegasAfterglow Documentation

This directory contains the documentation for the VegasAfterglow project.

## Building the Documentation

To build the documentation locally, you need to have the following dependencies installed:

- Python 3.7+
- Sphinx
- sphinx_rtd_theme
- breathe
- Doxygen
- Graphviz (for diagrams)

### Installing Dependencies

```bash
# Python dependencies
pip install sphinx sphinx_rtd_theme breathe

# On Ubuntu/Debian
sudo apt-get install doxygen graphviz

# On macOS
brew install doxygen graphviz
```

### Building

```bash
cd docs
make all
```

This will:
1. Run Doxygen to generate XML documentation from the C++ source code
2. Run Sphinx to generate HTML documentation, incorporating both the C++ API documentation (via Breathe) and the Python API documentation

The generated documentation will be available in the `build/html/` directory.

## Documentation Structure

- `source/` - Contains Sphinx RST files
- `Doxyfile` - Doxygen configuration
- `Makefile` - Build script for documentation

## Online Documentation

The latest documentation is automatically built and deployed to GitHub Pages on each push to the main branch:

https://yourusername.github.io/VegasAfterglow/ 