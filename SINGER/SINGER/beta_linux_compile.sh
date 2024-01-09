#!/bin/bash

# Check if a version number is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: bash $0 <version_number>"
    exit 1
fi

# Version number from the first argument
VERSION=$1

# Directory for the release
RELEASE_DIR="../../releases"
VERSION_DIR="$RELEASE_DIR/singer-$VERSION-beta-linux-x86_64"

# Create version directory
mkdir -p $VERSION_DIR

# Compile the program with optimizations and debugging information
g++ -std=c++17 -O3 -g -static *.cpp -o $VERSION_DIR/singer

# Compile the debug version of the program
g++ -std=c++17 -g -static *.cpp -o $VERSION_DIR/singer_debug

# Copy additional files
cp singer_master $VERSION_DIR/singer_master
cp convert_to_tskit $VERSION_DIR/convert_to_tskit
cp parallel_singer $VERSION_DIR/parallel_singer
cp index_vcf.py $VERSION_DIR/index_vcf.py
cp convert_long_ARG.py $VERSION_DIR/convert_long_ARG.py

# Change directory to releases
cd $RELEASE_DIR

# Create a tarball with the version number in the name
tar -cvzf "singer-$VERSION-beta-linux-x86_64.tar.gz" "singer-$VERSION-beta-linux-x86_64" 
#rm -rf "singer-$VERSION-beta-linux-x86_64" 
