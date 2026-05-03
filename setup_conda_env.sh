#!/bin/bash
# Setup script for CrobustaScreen conda environment
# This configures PyCall.jl and RCall.jl to use the conda-managed Python and R
# Supports: micromamba, mamba, and conda
#
# USAGE: ./setup_conda_env.sh (do not source this script)

set -e

# Check if being sourced (will crash terminal with exit commands)
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
    echo "ERROR: This script should be executed, not sourced!"
    echo "Correct usage: ./setup_conda_env.sh"
    echo "             or: bash setup_conda_env.sh"
    return 1 2>/dev/null || true
fi

echo "=== CrobustaScreen Conda Environment Setup ==="

# Check which conda implementation is available
if command -v micromamba &> /dev/null; then
    CONDA_CMD="micromamba"
    CONDA_TYPE="micromamba"
    echo "Using micromamba (fastest option)"
elif command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
    CONDA_TYPE="conda"
    echo "Using mamba for faster environment creation"
elif command -v conda &> /dev/null; then
    CONDA_CMD="conda"
    CONDA_TYPE="conda"
    echo "Using conda"
else
    echo "Error: No conda implementation found. Please install one of:"
    echo "  - Micromamba (fastest): https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html"
    echo "  - Mambaforge: https://github.com/conda-forge/miniforge#mambaforge"
    echo "  - Miniconda: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# Create/update the environment
echo "Creating/updating environment from environment.yml..."
if [ "$CONDA_TYPE" = "micromamba" ]; then
    # Micromamba syntax
    $CONDA_CMD create -f environment.yml -y || $CONDA_CMD update -f environment.yml -y
else
    # Conda/mamba syntax
    $CONDA_CMD env create -f environment.yml || $CONDA_CMD env update -f environment.yml
fi

# Activate the environment
echo "Activating crobustascreen environment..."
if [ "$CONDA_TYPE" = "micromamba" ]; then
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate crobustascreen
    # Get environment path for micromamba
    CONDA_PREFIX=$MAMBA_ROOT_PREFIX/envs/crobustascreen
else
    eval "$(conda shell.bash hook)"
    conda activate crobustascreen
    # Get environment path for conda/mamba
    CONDA_PREFIX=$(conda info --base)/envs/crobustascreen
fi
PYTHON_PATH=$CONDA_PREFIX/bin/python
R_HOME=$CONDA_PREFIX/lib/R

echo "Conda environment prefix: $CONDA_PREFIX"
echo "Python path: $PYTHON_PATH"
echo "R home: $R_HOME"

# Configure Julia packages
echo "Configuring Julia packages..."
julia --project=. << EOF
using Pkg

# Add required packages if not already present
println("Adding Julia packages...")
packages = ["PyCall", "RCall", "Conda"]
for pkg in packages
    if !haskey(Pkg.project().dependencies, pkg)
        Pkg.add(pkg)
    end
end

# Configure PyCall to use conda Python
ENV["PYTHON"] = "$PYTHON_PATH"
Pkg.build("PyCall")

# Configure RCall to use conda R
ENV["R_HOME"] = "$R_HOME"
Pkg.build("RCall")

# Verify configurations
using PyCall
using RCall
println("PyCall Python: ", PyCall.pyprogramname)
println("RCall R: ", RCall.Rhome)

# Install other Julia dependencies from Project.toml
Pkg.instantiate()
EOF

# Install R packages not available in conda
echo "Installing additional R packages..."
Rscript - << EOF
# Check and install packages from CRAN
cran_packages <- c("leiden", "umap")
for (pkg in cran_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg, repos = "https://cloud.r-project.org")
    }
}

# Install STRINGdb from Bioconductor if needed
if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    BiocManager::install("STRINGdb", update = FALSE, ask = FALSE)
}

# Configure reticulate to use conda Python
library(reticulate)
use_python("$PYTHON_PATH", required = TRUE)
py_config()
EOF

# Create activation script for future use
cat > activate_crobustascreen.sh << ACTIVATION
#!/bin/bash
# Quick activation script for CrobustaScreen environment
# Auto-detects micromamba/mamba/conda

if command -v micromamba &> /dev/null; then
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate crobustascreen
    CONDA_PREFIX=\$MAMBA_ROOT_PREFIX/envs/crobustascreen
elif command -v mamba &> /dev/null; then
    eval "\$(conda shell.bash hook)"
    conda activate crobustascreen
    CONDA_PREFIX=\$(conda info --base)/envs/crobustascreen
elif command -v conda &> /dev/null; then
    eval "\$(conda shell.bash hook)"
    conda activate crobustascreen
    CONDA_PREFIX=\$(conda info --base)/envs/crobustascreen
else
    echo "Error: No conda implementation found"
    exit 1
fi

# Set environment variables
export PYTHON=\$CONDA_PREFIX/bin/python
export R_HOME=\$CONDA_PREFIX/lib/R
export JULIA_PROJECT=@.
export RETICULATE_PYTHON=\$CONDA_PREFIX/bin/python

# For CUDA support
export CUDA_HOME=\$CONDA_PREFIX
export LD_LIBRARY_PATH=\$CONDA_PREFIX/lib:\$LD_LIBRARY_PATH

echo "CrobustaScreen environment activated"
echo "Python: \$PYTHON"
echo "R: \$R_HOME"
echo "Julia project: \$JULIA_PROJECT"
ACTIVATION

chmod +x activate_crobustascreen.sh

echo ""
echo "=== Setup Complete ==="
echo ""
echo "To activate this environment in the future, run:"
echo "  source ./activate_crobustascreen.sh"
echo ""
echo "Or manually:"
echo "  conda activate crobustascreen"
echo "  export JULIA_PROJECT=@."
echo ""