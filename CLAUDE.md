# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CrobustaScreen is a pipeline for self-supervised phenotype detection from confocal images of *Ciona robusta* embryos. It uses autoencoders for dimension reduction and clustering algorithms to identify phenotypic patterns in embryo segmentation data.

## Development Environment

### Setup with Conda/Mamba/Micromamba (Recommended)
```bash
# Install micromamba (fastest) or mamba/conda
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)

# Run the setup script (auto-detects micromamba/mamba/conda)
chmod +x setup_conda_env.sh
./setup_conda_env.sh

# Activate the environment
source ./activate_crobustascreen.sh

# Verify installation
julia verify_environment.jl
Rscript verify_environment.R
```

### Alternative: Setup with Nix
```bash
nix develop .
```

## Common Commands

### Full Pipeline
```bash
# Generate all prerequisites (interactions, normalized data)
make all

# Train autoencoder (outputs to data/)
julia autoencoder.jl --path "data/"

# Cluster and visualize
Rscript cluster.R
Rscript plot.clust.R
```

### Preprocessing (run only if data changes)
```bash
Rscript readEmbryos.R          # Parse segdat/ -> data/embryodat.csv
Rscript readPheno.R            # Parse imaris.csv -> data/params.csv
julia preprocess.jl            # Normalize -> data/X.csv
```

### Clustering Options
```bash
# Default: use autoencoder embeddings
Rscript cluster.R

# With custom hyperparameter ranges
Rscript cluster.R -k 5 -K 40 -G 2.0 -l 500

# Using PCA instead of autoencoder
Rscript cluster.R -e data/PCs.csv -o data/PCA
```

### Visualization Options
```bash
# Default: select clusters by ES (enrichment score)
Rscript plot.clust.R

# Select by combined_score (ES * recall)
Rscript plot.clust.R -s combined_score

# With custom paths
Rscript plot.clust.R -e data/E.csv -c data/ -o fig/
```

### DEWAK Analysis
```bash
julia dewak.jl
```

## Architecture

### Data Flow
```
Raw segmentation (segdat/)
    → readEmbryos.R → data/embryodat.csv
    → readPheno.R → data/params.csv, data/z_dat.csv
    → preprocess.jl → data/X.csv (116 features, scaled [-1,1])
    → autoencoder.jl → data/E.csv (14-dim embeddings)
    → cluster.R → data/k.csv, data/leiden.csv
    → plot.clust.R → fig/ (UMAP, heatmaps, networks)
```

### Autoencoder Architecture
- Input: 116 normalized features
- Hidden layers: 58 → 29 → 14 (bottleneck) → 29 → 58
- Activation: tanh
- Training: 10,000 epochs, AdamW (η=0.0001, λ=0.0001), 10% test split
- Also trains a Sparse Autoencoder (SAE) variant, output in `data/SAE/`

### Hyperparameter Optimization
- **k selection** (k-NN graph): k ∈ [3, 53], selected by GSEA enrichment score
- **γ selection** (Leiden resolution): γ ∈ [0.01, 3.0], 1000 random samples
- **Metrics**: ES (enrichment), recall (vs known interactions), silhouette width, k-NN classifier error

### Key R Modules (`R/`)
- `io.R`: Data loading and command-line parsing (`data.parser()`, `parse.env()`)
- `optimization.R`: Hyperparameter search (`get.k()`, `get.res.unif()`)
- `leiden.R`: Leiden clustering wrapper
- `knn.R`: k-NN graph construction
- `gene.network.R`: Infer gene network from k-NN enrichment
- `hyper.R`: Hypergeometric enrichment tests
- `clustplots.R`, `plotfns.R`: Visualization utilities
- `dirfns.R`: Output directory management (`dir.csv()`, `dir.pdf()`, `dir.plot()`)

### Julia Custom Packages
Required packages not in standard registries (installed via setup scripts or flake.nix):
- [Autoencoders.jl](https://github.com/kewiechecki/Autoencoders.jl)
- [TrainingIO.jl](https://github.com/kewiechecki/TrainingIO.jl)
- [DeePWAK.jl](https://github.com/kewiechecki/DeePWAK.jl) (for DEWAK analysis)

### Protein Interaction Network
Built from orthologs via ENSEMBL and STRINGdb:
```bash
Rscript cint.ensembl.R     # Fetch C. robusta orthologs
Rscript STRINGdb.R         # Download interactions
Rscript get.interactions.R # Merge into data/interactions.csv
```

## Key Data Files
- `data/X.csv`: Normalized input features (n × 116)
- `data/E.csv`: Autoencoder embeddings (n × 14)
- `data/groups.csv`: Sample metadata (Condition column)
- `data/interactions.csv`: Known protein interactions for validation
- `data/k.csv`: k optimization results
- `data/leiden.csv`: Clustering results with all metrics

## Troubleshooting

### PyCall/RCall Issues
Ensure Julia uses the conda environment's Python/R:
```julia
ENV["PYTHON"] = "/path/to/conda/env/bin/python"
ENV["R_HOME"] = "/path/to/conda/env/lib/R"
```
The setup scripts handle this automatically.

### Missing Julia Packages
Custom packages must be installed from GitHub:
```julia
using Pkg
Pkg.add(url="https://github.com/kewiechecki/Autoencoders.jl")
Pkg.add(url="https://github.com/kewiechecki/TrainingIO.jl")
```
