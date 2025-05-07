# Documentation for VegasAfterglow

This directory contains the documentation for the VegasAfterglow project.

## Building the Documentation

### Prerequisites

- Python 3.9 or later
- Doxygen (for C++ API documentation)
- Graphviz (for diagrams)
- Sphinx and related packages

### Install Dependencies

```bash
pip install sphinx sphinx-rtd-theme breathe
```

### Build Process

You can build the documentation using the provided script:

```bash
cd docs
chmod +x build_docs.sh
./build_docs.sh
```

The built documentation will be in the `build/html` directory.

## Documentation Structure

- `source/` - Contains the RST files for Sphinx
- `doxygen/` - Output directory for Doxygen (generated)
- `build/` - Output directory for Sphinx (generated)
- `Doxyfile` - Configuration for Doxygen
- `source/conf.py` - Configuration for Sphinx

## Automatic Deployment

The documentation is automatically built and deployed to the GitHub Wiki whenever changes are pushed to the `main` branch and affect files in the `docs/`, `include/`, `src/`, or `pybind/` directories.

### Manual Deployment

To manually trigger the documentation build and deployment, you can run:

```bash
gh workflow run docs.yml
```

## Documentation Guidelines

### C++ Documentation

Use Doxygen-style comments for C++ code:

```cpp
/**
 * @brief Brief description of the function/class
 * @details Detailed description
 * 
 * @param param1 Description of first parameter
 * @param param2 Description of second parameter
 * @return Description of return value
 */
```

### Python Documentation

Use Google-style docstrings for Python code:

```python
def function(param1, param2):
    """Brief description of the function.
    
    Detailed description of the function.
    
    Args:
        param1: Description of first parameter
        param2: Description of second parameter
        
    Returns:
        Description of return value
    """
```

### Template and Inline Functions

For template and inline functions, ensure you document both template parameters and function parameters:

```cpp
/**
 * @brief Brief description of the template function
 * @details Detailed description
 *
 * @tparam T Description of template parameter
 * @param param1 Description of first parameter
 * @return Description of return value
 */
template<typename T>
inline T function(T param1) { return param1; }
```

## Customizing Documentation

- Edit `source/conf.py` to modify Sphinx configuration
- Edit `Doxyfile` to modify Doxygen configuration 