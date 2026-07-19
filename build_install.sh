#!/bin/bash

# Build script for hoomd-bussi-thermostat
# 
# Usage:
#   ./build.sh           # Build with GPU support (default)
#   ./build.sh --no-gpu  # Build without GPU support

# Check HOOMD-blue version compatibility
echo "Checking HOOMD-blue version compatibility..."
HOOMD_VERSION=$(python3 -c "import hoomd; print(hoomd.version.version)" 2>/dev/null)
if [ $? -ne 0 ]; then
    echo "ERROR: HOOMD-blue is not installed or not accessible."
    echo "Please install HOOMD-blue 5.2.0 first following the documentation."
    exit 1
fi

if [ "$HOOMD_VERSION" != "5.2.0" ] && [ "$HOOMD_VERSION" != "5.4.0" ] && [ "$HOOMD_VERSION" != "7.1.0" ]; then
    echo "ERROR: HOOMD-blue version mismatch."
    echo "Required: 5.2.0, 5.4.0, or 7.1.0"
    echo "Found: $HOOMD_VERSION"
    echo "Please install a supported HOOMD-blue version following the documentation."
    exit 1
fi

echo "HOOMD-blue version $HOOMD_VERSION confirmed"

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
