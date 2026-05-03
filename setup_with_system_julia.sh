#!/bin/bash
# Alternative setup using system Julia instead of conda Julia
# This avoids the executable stack issue

echo "=== Setting up with System Julia ==="

# Check for system Julia
if ! command -v julia &> /dev/null; then
    echo "System Julia not found. Please install it first:"
    echo ""
    echo "Option 1: Using juliaup (recommended):"
    echo "  curl -fsSL https://install.julialang.org | sh"
    echo "  source ~/.bashrc"
    echo ""
    echo "Option 2: Package manager:"
    echo "  # Fedora/RHEL:"
    echo "  sudo dnf install julia"
    echo "  # Ubuntu/Debian:"
    echo "  sudo apt-get install julia"
    echo "  # Arch:"
    echo "  sudo pacman -S julia"
    exit 1
fi

echo "Found system Julia: $(which julia)"
echo "Julia version: $(julia --version)"

# Create conda environment WITHOUT Julia
echo "Creating conda environment without Julia..."
cat > environment_no_julia.yml << 'EOF'
name: crobustascreen
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  # Python 3.11 for better compatibility
  - python=3.11
  - pip
  
  # Python packages for clustering
  - python-igraph
  - leidenalg
  - umap-learn
  - numpy
  - scipy
  - matplotlib
  
  # R 4.3 and essential packages
  - r-base=4.3
  - r-essentials
  - r-devtools
  - r-biocmanager
  
  # R packages for the pipeline
  - r-optparse
  - r-purrr
  - r-ggplot2
  - r-ggpubr
  - r-class
  - r-cluster
  - r-igraph
  - r-reticulate
  - r-circlize
  - bioconductor-complexheatmap
  - bioconductor-biomart
  - bioconductor-fgsea
  
  # NO JULIA HERE - using system version
  
  # CUDA support (adjust version as needed)
  - cudatoolkit=11.8
  - cudnn=8.9
  
  # Build tools and libraries
  - gcc
  - gxx
  - make
  - cmake
  - libxml2
  - libcurl
  - openssl
  - libpng
  - libuv
  - gsl
  - bzip2
  - icu
  
  # Additional utilities
  - git
  - wget
  - curl
  
  # Pip packages
  - pip:
    - pandas
    - seaborn
EOF

# Create/update environment
if command -v micromamba &> /dev/null; then
    micromamba create -f environment_no_julia.yml -y
    echo "Environment created. Activate with:"
    echo "  micromamba activate crobustascreen"
elif command -v mamba &> /dev/null; then
    mamba env create -f environment_no_julia.yml
    echo "Environment created. Activate with:"
    echo "  conda activate crobustascreen"
else
    conda env create -f environment_no_julia.yml
    echo "Environment created. Activate with:"
    echo "  conda activate crobustascreen"
fi

echo ""
echo "After activating, configure Julia packages with:"
echo "  julia --project=. configure_julia_system.jl"