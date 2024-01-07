#!/bin/bash

# Check if a version number is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: bash $0 <version_number>"
    exit 1
fi

# Version number from the first argument
VERSION=$1

# Compile the program with optimizations and debugging information
g++ -std=c++17 -O3 -g -static *.cpp -o ../../releases/singer

# Compile the debug version of the program
g++ -std=c++17 -g -static *.cpp -o ../../releases/singer_debug

# Copy additional files
cp singer_master ../../releases/singer_master
cp convert_to_tskit ../../releases/convert_to_tskit

# Create a tarball with the version number in the name
cd ../../releases
tar -cvzf singer-$VERSION-beta-linux-x86_64.tar.gz \
    singer_debug \
    singer \
    convert_to_tskit \
    singer_master
