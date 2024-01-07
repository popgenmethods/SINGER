#!/bin/bash

# Check if a version number is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: bash $0 <version_number>"
    exit 1
fi

# Version number from the first argument
VERSION=$1

# Compile the program with optimizations and debugging information
clang++ -std=c++17 -O3 -g -c *.cpp -o ../../releases/singer

# Compile the debug version of the program
clang++ -std=c++17 -g -c *.cpp -o ../../releases/singer_debug

# Copy additional files
cp singer_master ../../releases/singer_master
cp convert_to_tskit ../../releases/convert_to_tskit

# Change directory to releases
cd ../../releases

# Create a tarball with the version number in the name
# Rename to indicate ARM architecture (macOS M1)
tar -cvzf singer-$VERSION-beta-mac-arm64.tar.gz \
    singer_debug \
    singer \
    convert_to_tskit \
    singer_master
