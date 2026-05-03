# Conda Environment Setup for CrobustaScreen

This setup provides better control over Python and R environments when using PyCall.jl and RCall.jl, addressing cross-language dependency management.

## Why Conda?

The CrobustaScreen pipeline requires tight integration between:
- **Julia** calling Python (via PyCall.jl) for leiden clustering
- **Julia** calling R (via RCall.jl) for statistical analysis  
- **R** calling Python (via reticulate) for leiden and UMAP

Conda provides a unified environment where all three languages share the same Python interpreter and R installation, preventing version conflicts.

## Quick Start

1. **Install a Conda Implementation** (choose one):
   
   **Micromamba** (Fastest, ~5MB, recommended):
   ```bash
   "${SHELL}" <(curl -L micro.mamba.pm/install.sh)
   ```
   
   **Mambaforge** (Fast, includes conda):
   ```bash
   wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
   bash Mambaforge-Linux-x86_64.sh
   ```
   
   **Miniconda** (Traditional):
   ```bash
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   ```

2. **Run Setup** (Two safer options):

   **Option A: Two-step process (Recommended)**
   ```bash
   # Step 1: Create environment (safe to source or execute)
   source safe_setup.sh
   
   # Step 2: Activate and configure
   micromamba activate crobustascreen  # or conda/mamba
   ./configure_julia_r.sh
   ```
   
   **Option B: Original script (DO NOT source!)**
   ```bash
   chmod +x setup_conda_env.sh
   ./setup_conda_env.sh  # Execute, don't source! Will crash terminal if sourced
   ```

3. **Activate Environment** (for future sessions):
   ```bash
   source ./activate_crobustascreen.sh
   ```

4. **Verify Setup**:
   ```bash
   julia verify_environment.jl    # Tests Julia→Python and Julia→R
   Rscript verify_environment.R    # Tests R→Python
   ```

## What the Setup Does

1. Creates a conda environment with Python 3.11, R 4.3, and Julia 1.10
2. Installs all required packages (leiden, igraph, umap, etc.)
3. Configures PyCall.jl to use conda's Python
4. Configures RCall.jl to use conda's R
5. Configures R's reticulate to use conda's Python
6. Sets up CUDA support if available

## Environment Structure

```
$CONDA_PREFIX/
├── bin/
│   ├── python          # Single Python for all languages
│   ├── R               # Single R installation
│   └── julia           # Julia executable
├── lib/
│   ├── python3.11/     # Python packages
│   │   └── site-packages/
│   │       ├── leidenalg/
│   │       ├── igraph/
│   │       └── umap/
│   └── R/              # R packages
│       └── library/
│           ├── leiden/
│           ├── reticulate/
│           └── igraph/
```

## Troubleshooting

### PyCall using wrong Python
```julia
ENV["PYTHON"] = "/path/to/conda/env/bin/python"
using Pkg
Pkg.build("PyCall")
```

### RCall using wrong R
```julia
ENV["R_HOME"] = "/path/to/conda/env/lib/R"
using Pkg
Pkg.build("RCall")
```

### Reticulate using wrong Python
```r
library(reticulate)
use_python("/path/to/conda/env/bin/python", required = TRUE)
```

## Files Created

- `environment.yml` - Conda environment specification
- `setup_conda_env.sh` - One-time setup script
- `activate_crobustascreen.sh` - Quick activation script
- `verify_environment.jl` - Julia verification tests
- `verify_environment.R` - R verification tests

## Comparison with Nix

| Feature | Conda | Nix |
|---------|-------|-----|
| Setup complexity | Simple | Complex |
| Reproducibility | Good | Excellent |
| Cross-language control | Excellent | Good |
| Community adoption | High | Low |
| Package availability | Excellent | Good |
| Binary caching | Yes | Yes |
| Rollback capability | Limited | Full |

For this project, Conda is recommended because:
1. Explicit control over Python/R versions used by all languages
2. Simpler debugging of cross-language issues
3. Wider adoption in scientific computing
4. Better documentation for multi-language setups