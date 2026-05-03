#!/bin/bash
# Safer setup script that won't crash your terminal
# Can be either sourced or executed

echo "=== CrobustaScreen Safe Setup ==="

# Function to safely handle errors without exiting shell
safe_error() {
    echo "ERROR: $1" >&2
    if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
        # Being sourced - return instead of exit
        return 1 2>/dev/null || true
    else
        # Being executed - safe to exit
        exit 1
    fi
}

# Detect conda implementation
detect_conda() {
    if command -v micromamba &> /dev/null; then
        echo "micromamba"
    elif command -v mamba &> /dev/null; then
        echo "mamba"
    elif command -v conda &> /dev/null; then
        echo "conda"
    else
        echo "none"
    fi
}

CONDA_IMPL=$(detect_conda)

if [[ "$CONDA_IMPL" == "none" ]]; then
    safe_error "No conda implementation found. Please install micromamba, mamba, or conda first."
fi

echo "Found: $CONDA_IMPL"

# Create or update environment
create_environment() {
    local impl=$1
    
    echo "Creating/updating environment..."
    
    if [[ "$impl" == "micromamba" ]]; then
        micromamba create -f environment.yml -y || \
        micromamba update -f environment.yml -y || \
        safe_error "Failed to create/update environment"
    else
        $impl env create -f environment.yml || \
        $impl env update -f environment.yml || \
        safe_error "Failed to create/update environment"
    fi
}

# Only create env if not already activated
if [[ "$CONDA_DEFAULT_ENV" != "crobustascreen" ]]; then
    create_environment "$CONDA_IMPL"
    
    echo ""
    echo "Environment created successfully!"
    echo ""
    echo "To activate, run:"
    if [[ "$CONDA_IMPL" == "micromamba" ]]; then
        echo "  micromamba activate crobustascreen"
    else
        echo "  conda activate crobustascreen"
    fi
    echo ""
    echo "Then run the configuration:"
    echo "  ./configure_julia_r.sh"
else
    echo "Already in crobustascreen environment"
fi

echo ""
echo "For quick activation in the future:"
echo "  source ./activate_crobustascreen.sh"