#!/bin/bash
set -e  # Exit on error

# Path setup
DOCS_DIR=$(pwd)
PROJECT_ROOT=$(dirname "$DOCS_DIR")

echo "=== Cleaning existing build files ==="
rm -rf doxygen build

# Create directory structure
mkdir -p source/_static/css
mkdir -p doxygen

echo "=== Running Doxygen to generate XML for Breathe ==="
doxygen Doxyfile

# Verify XML files exist
if [ ! -f "doxygen/xml/index.xml" ]; then
    echo "ERROR: Doxygen XML not generated correctly. Missing index.xml file."
    exit 1
fi

echo "=== Debug: Listing some key XML files ==="
ls -la doxygen/xml/index.xml
ls -la doxygen/xml/class_gaussian_jet.xml || echo "GaussianJet XML missing!"
ls -la doxygen/xml/class_tophat_jet.xml || echo "TophatJet XML missing!"

# Create a simple custom CSS file
cat > source/_static/css/custom.css << EOF
/* Basic styling for C++ documentation */
dl.cpp.function {
    margin-bottom: 15px;
    padding: 10px;
    border-radius: 5px;
    background-color: #f7f7f7;
}

dl.cpp.class {
    padding: 10px;
    margin: 10px 0;
    border: 1px solid #eee;
    border-radius: 5px;
}
EOF

echo "=== Building Sphinx documentation ==="
sphinx-build -b html source build/html

# Copy Doxygen HTML for reference
echo "=== Adding source code browsing capability ==="
cp -r doxygen/html build/html/doxygen

echo "=== Documentation build complete ==="

