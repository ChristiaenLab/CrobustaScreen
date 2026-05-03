#!/bin/bash
# Fix Julia executable stack issue on SELinux/security-hardened systems

echo "=== Fixing Julia Executable Stack Issue ==="

# Check if we're in the conda environment
if [[ -z "$CONDA_PREFIX" ]]; then
    echo "Please activate the crobustascreen environment first:"
    echo "  micromamba activate crobustascreen"
    exit 1
fi

JULIA_LIB_DIR="$CONDA_PREFIX/lib/julia"

# Method 1: Try to fix with execstack (if available)
if command -v execstack &> /dev/null; then
    echo "Using execstack to clear executable stack flag..."
    execstack -c "$JULIA_LIB_DIR"/*.so
    echo "Fixed with execstack"
    exit 0
fi

# Method 2: Try patchelf (if available)
if command -v patchelf &> /dev/null; then
    echo "Using patchelf to clear executable stack flag..."
    for lib in "$JULIA_LIB_DIR"/*.so; do
        patchelf --remove-needed libexecstack.so "$lib" 2>/dev/null || true
    done
    echo "Fixed with patchelf"
    exit 0
fi

# Method 3: SELinux context (if SELinux is enabled)
if command -v getenforce &> /dev/null && [ "$(getenforce)" != "Disabled" ]; then
    echo "SELinux is enabled. Trying to set permissive context..."
    if command -v chcon &> /dev/null; then
        chcon -t textrel_shlib_t "$JULIA_LIB_DIR"/*.so 2>/dev/null || {
            echo "Could not change SELinux context. You may need to run:"
            echo "  sudo setsebool -P selinuxuser_execstack 1"
            echo "OR temporarily:"
            echo "  sudo setenforce 0  # Run your Julia code"
            echo "  sudo setenforce 1  # Re-enable after"
        }
    fi
fi

echo ""
echo "If the above didn't work, try one of these solutions:"
echo ""
echo "1. Install execstack and run this script again:"
echo "   # On Fedora/RHEL/CentOS:"
echo "   sudo dnf install prelink"
echo "   # On Ubuntu/Debian:"
echo "   sudo apt-get install execstack"
echo ""
echo "2. Use system Julia instead of conda Julia:"
echo "   # Remove Julia from environment.yml and install system-wide"
echo "   sudo dnf install julia  # or apt-get"
echo ""
echo "3. Temporarily disable SELinux (not recommended for production):"
echo "   sudo setenforce 0"
echo "   # Run your Julia commands"
echo "   sudo setenforce 1"
echo ""
echo "4. Use Docker/Podman container where this isn't an issue"