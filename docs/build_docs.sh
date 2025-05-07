#!/bin/bash
set -e  # Exit on error

# Path setup
DOCS_DIR=$(pwd)
PROJECT_ROOT=$(dirname "$DOCS_DIR")

echo "=== Cleaning existing build files ==="
rm -rf doxygen build

echo "=== Running Doxygen to generate XML ==="
doxygen Doxyfile

# Verify XML files exist
if [ ! -f "doxygen/xml/index.xml" ]; then
    echo "ERROR: Doxygen XML not generated correctly. Missing index.xml file."
    exit 1
fi

echo "=== Debug: Listing top-level XML files ==="
ls -la doxygen/xml/

echo "=== Debug: Verify key class files ==="
ls -la doxygen/xml/class_gaussian_jet.xml || echo "GaussianJet XML missing!"
ls -la doxygen/xml/class_tophat_jet.xml || echo "TophatJet XML missing!"

echo "=== Running Sphinx to generate HTML (verbose mode) ==="
sphinx-build -b html -v source build/html

echo "=== Documentation build complete ==="
echo "HTML documentation is in build/html"
echo "Doxygen XML is in doxygen/xml" 