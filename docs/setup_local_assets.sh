#!/bin/bash
# Script to set up assets for local development

# Create required directories
mkdir -p build/html/_static
mkdir -p build/html/_images
mkdir -p build/html/_static/assets
mkdir -p source/_static

# Copy logo from assets to needed locations
echo "Copying logo files to build locations..."
cp ../assets/logo.svg build/html/_static/
cp ../assets/logo.svg build/html/_images/
cp ../assets/logo.svg source/_static/
cp -r ../assets/* build/html/_static/assets/

echo "Logo files copied successfully!"
echo "Use this script before running make html if you're having logo issues." 