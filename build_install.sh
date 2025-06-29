#!/bin/bash

# Build script for hoomd-bussi-thermostat
# 
# Usage:
#   ./build.sh           # Build with GPU support (default)
#   ./build.sh --no-gpu  # Build without GPU support

# Parse command line arguments
ENABLE_GPU=ON
for arg in "$@"; do
    case $arg in
        --no-gpu)
        ENABLE_GPU=OFF
        shift
        ;;
        --help|-h)
        echo "Usage: $0 [--no-gpu] [--help]"
        echo "  --no-gpu    Disable GPU/CUDA support"
        echo "  --help      Show this help message"
        exit 0
        ;;
    esac
done

echo "Building with GPU support: $ENABLE_GPU"

rm -rf build
cmake -B build -S . -DENABLE_GPU=$ENABLE_GPU
cmake --build build
cmake --install build
