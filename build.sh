#!/bin/bash

# Create build directory if it doesn't exist
mkdir -p build

# Navigate to build directory
cd build

# Configure with CMake
cmake ..

# Build the project
make -j$(nproc)

echo "Build completed! CFD executables are in build/bin/" 