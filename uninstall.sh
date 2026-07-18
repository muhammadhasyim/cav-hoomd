#!/bin/bash

# Uninstall script for hoomd-bussi-thermostat and cavitymd plugins
# This script removes the installed files from the HOOMD package directory

echo "HOOMD Plugin Uninstaller"
echo "========================"

# Find HOOMD installation directory
HOOMD_DIR=$(python -c "import hoomd; import os; print(os.path.dirname(hoomd.__file__))" 2>/dev/null)

if [ $? -ne 0 ] || [ -z "$HOOMD_DIR" ]; then
    echo "Error: Could not find HOOMD installation directory."
    echo "Make sure HOOMD is installed and accessible from Python."
    exit 1
fi

echo "Found HOOMD installation at: $HOOMD_DIR"
echo

# Function to safely remove files/directories
safe_remove() {
    local target="$1"
    if [ -e "$target" ]; then
        echo "Removing: $target"
        rm -rf "$target"
        if [ $? -eq 0 ]; then
            echo "  ✓ Successfully removed"
        else
            echo "  ✗ Failed to remove"
            return 1
        fi
    else
        echo "Not found (already removed): $target"
    fi
    return 0
}

# Check if running with appropriate permissions
if [ ! -w "$HOOMD_DIR" ]; then
    echo "Warning: You may not have write permissions to $HOOMD_DIR"
    echo "You might need to run this script with sudo or as an administrator."
    echo
fi

echo "Removing bussi_reservoir plugin..."
echo "-----------------------------------"

# Remove bussi_reservoir directory and all its contents
safe_remove "$HOOMD_DIR/bussi_reservoir"

echo
echo "Removing cavitymd plugin..."
echo "---------------------------"

# Remove cavitymd directory and all its contents  
safe_remove "$HOOMD_DIR/cavitymd"

echo
echo "Cleanup complete!"
echo

# Verify removal
echo "Verifying removal..."
REMAINING_FILES=0

if [ -d "$HOOMD_DIR/bussi_reservoir" ]; then
    echo "Warning: bussi_reservoir directory still exists"
    REMAINING_FILES=$((REMAINING_FILES + 1))
fi

if [ -d "$HOOMD_DIR/cavitymd" ]; then
    echo "Warning: cavitymd directory still exists"
    REMAINING_FILES=$((REMAINING_FILES + 1))
fi

if [ $REMAINING_FILES -eq 0 ]; then
    echo "✓ All plugin files have been successfully removed."
    echo
    echo "The bussi_reservoir and cavitymd plugins have been uninstalled."
    echo "You may need to restart Python to clear any cached imports."
else
    echo "✗ Some files may not have been removed completely."
    echo "  This could be due to permission issues or files in use."
    echo "  You may need to remove them manually or with elevated privileges."
fi

echo
echo "Uninstall process completed." 