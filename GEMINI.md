# CrobustaScreen Project Context

## Project Overview
**CrobustaScreen** is a computational pipeline for self-supervised phenotype detection from confocal images of *Ciona robusta* embryos. It processes segmentation data to identify phenotypic patterns using autoencoders for dimension reduction and graph-based clustering algorithms.

The project integrates **Julia** (machine learning/autoencoders), **R** (preprocessing, clustering, visualization), and **Python** (underlying libraries for clustering/UMAP).

## Environment Setup
Due to the multi-language nature of the project (Julia, R, Python), environment management is critical.

### Recommended: Conda/Mamba
The project uses a Conda environment to ensure all languages share the same Python and R executables, preventing version conflicts.

1.  **Create Environment**:
    ```bash
    # Use the safe setup script (recommended)
    source safe_setup.sh
    micromamba activate crobustascreen
    ./configure_julia_r.sh
    ```
    *Alternatively, use `setup_conda_env.sh` (do not source).*

2.  **Activate**:
    ```bash
    source ./activate_crobustascreen.sh
    ```

3.  **Verify**:
    ```bash
    julia verify_environment.jl
    Rscript verify_environment.R
    ```

### Alternative: Nix
A `flake.nix` file is provided for Nix users.
```bash
nix develop .
```

## Key Workflows & Commands
The pipeline is automated via `Makefile`, but steps can be run individually.

### 1. Build All Prerequisites
Generates all necessary processed data files from raw inputs.
```bash
make all
```

### 2. Preprocessing
Parses raw segmentation data and experimental metadata.
-   **Embryo Data**: `Rscript readEmbryos.R` (Parses `segdat/` -> `data/embryodat.csv`)
-   **Phenotypes**: `Rscript readPheno.R` (Parses `imaris.csv` -> `data/params.csv`)
-   **Normalization**: `julia preprocess.jl` (Normalizes stats -> `data/X.csv`)

### 3. Dimensionality Reduction (Autoencoder)
Trains an autoencoder to compress phenotype data into embeddings.
```bash
# Trains model and saves embeddings to data/E.csv
julia autoencoder.jl --path "data/"
```

### 4. Clustering & Analysis
Clusters the embeddings and generates visualizations.
-   **Clustering**: `Rscript cluster.R` (Uses Leiden algorithm on `data/E.csv`)
-   **Visualization**: `Rscript plot.clust.R` (Generates figures in `fig/`)
-   **DEWAK**: `julia dewak.jl` (Deep weighted averaging kernel analysis)

## Architecture & Data Flow

1.  **Input**: Raw confocal image segmentation statistics (in `segdat/`) and phenotype labels (`imaris.csv`).
2.  **Preprocessing**: Data is aggregated per embryo, normalized (z-score), and scaled to [-1, 1].
3.  **Embedding**: `Autoencoders.jl` compresses the 116-dimensional feature space.
4.  **Graph Construction**: A k-NN graph is built from embeddings.
5.  **Clustering**: The Leiden algorithm partitions the graph into phenotypic clusters.
6.  **Validation**: Clusters are validated against known protein interactions (fetched via `STRINGdb` and `biomaRt`) and experimental perturbations.

## Key Files
-   `autoencoder.jl`: Main Julia script for training the autoencoder.
-   `cluster.R`: R script for graph-based clustering and hyperparameter tuning.
-   `preprocess.jl`: Julia script for data normalization.
-   `readEmbryos.R`, `readPheno.R`: R scripts for parsing raw inputs.
-   `get.interactions.R`: Fetches and processes protein interaction networks.
-   `Makefile`: Orchestrates the build process.
-   `Project.toml`: Julia dependencies.
-   `environment.yml` (generated): Conda environment spec.
-   `CLAUDE.md`: Quick reference for development commands.

## Troubleshooting
-   **PyCall/RCall Errors**: Ensure `ENV["PYTHON"]` and `ENV["R_HOME"]` in Julia point to the Conda environment paths. The `setup_conda_env.sh` script handles this automatically.
-   **Missing Dependencies**: If `Autoencoders` or `TrainingIO` are missing in Julia, they are likely installed from GitHub. Check `Manifest.toml` or the setup scripts.
