#!/usr/bin/env Rscript
# Setup script for renv package management
# This captures all R dependencies and creates a reproducible environment

# Install renv if not already installed
if (!requireNamespace("renv", quietly = TRUE)) {
  # Try to install in user library
  user_lib <- Sys.getenv("R_LIBS_USER")
  if (user_lib == "") {
    user_lib <- file.path(Sys.getenv("HOME"), "R", 
                          paste0(R.version$platform, "-library"), 
                          paste(R.version$major, strsplit(R.version$minor, "\\.")[[1]][1], sep = "."))
  }
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_lib, .libPaths()))
  install.packages("renv", repos = "https://cloud.r-project.org", lib = user_lib)
}

# Initialize renv (creates renv/ directory and .Rprofile)
cat("Initializing renv for CrobustaScreen project...\n")
renv::init(
  settings = list(
    # Use conda's R packages as external library
    external.libraries = file.path(Sys.getenv("CONDA_PREFIX"), "lib", "R", "library"),
    # Don't use a global cache to avoid conflicts with conda
    use.cache = FALSE,
    # Automatically snapshot after installs
    snapshot.type = "explicit"
  )
)

# List of required packages based on the codebase
required_packages <- c(
  # Core data manipulation and utilities
  "optparse",    # Command line options
  "purrr",       # Functional programming
  
  # Clustering and graph analysis
  "igraph",      # Graph/network analysis
  "leiden",      # Leiden clustering algorithm
  "umap",        # UMAP dimension reduction
  "cluster",     # Clustering algorithms
  "class",       # Classification methods (for k-NN)
  
  # Visualization
  "ggplot2",     # Advanced plotting
  "ggpubr",      # Publication ready plots
  "circlize",    # Circular visualizations
  "ComplexHeatmap", # Complex heatmaps
  
  # Bioinformatics
  "biomaRt",     # Access to BioMart databases
  "fgsea",       # Gene set enrichment analysis
  "STRINGdb",    # STRING protein interaction database
  "BiocManager", # Bioconductor package management
  
  # Python integration
  "reticulate",  # Python interoperability
  
  # Statistical methods
  "parallel"     # Parallel computing (usually part of base R)
)

# Check which packages need to be installed
cat("\nChecking package requirements...\n")
missing_packages <- c()
bioc_packages <- c("biomaRt", "fgsea", "STRINGdb", "ComplexHeatmap")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing_packages <- c(missing_packages, pkg)
    cat("  [ ] Missing:", pkg, "\n")
  } else {
    cat("  [✓] Found:", pkg, "\n")
  }
}

# Install missing packages
if (length(missing_packages) > 0) {
  cat("\nInstalling missing packages...\n")
  
  # Separate Bioconductor and CRAN packages
  bioc_to_install <- intersect(missing_packages, bioc_packages)
  cran_to_install <- setdiff(missing_packages, bioc_packages)
  
  # Install CRAN packages
  if (length(cran_to_install) > 0) {
    cat("Installing from CRAN:", paste(cran_to_install, collapse = ", "), "\n")
    install.packages(cran_to_install, repos = "https://cloud.r-project.org")
  }
  
  # Install Bioconductor packages
  if (length(bioc_to_install) > 0) {
    cat("Installing from Bioconductor:", paste(bioc_to_install, collapse = ", "), "\n")
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(bioc_to_install, update = FALSE, ask = FALSE)
  }
}

# Python packages via reticulate (for leiden, umap, etc.)
cat("\nConfiguring Python integration via reticulate...\n")
library(reticulate)

# Use conda Python if available
conda_python <- file.path(Sys.getenv("CONDA_PREFIX"), "bin", "python")
if (file.exists(conda_python)) {
  use_python(conda_python, required = TRUE)
  cat("Using conda Python:", conda_python, "\n")
} else {
  cat("Note: Conda Python not found, using system Python\n")
}

# Check Python packages
py_packages <- c("igraph", "leidenalg", "umap-learn", "numpy", "scipy")
cat("\nChecking Python packages for R integration:\n")
for (pkg in py_packages) {
  if (py_module_available(pkg)) {
    cat("  [✓] Found:", pkg, "\n")
  } else {
    cat("  [ ] Missing:", pkg, "(install with: pip install", pkg, ")\n")
  }
}

# Create a snapshot of the current state
cat("\nCreating renv snapshot...\n")
renv::snapshot(prompt = FALSE)

cat("\n=== renv Setup Complete ===\n")
cat("\nrenv has been initialized with the following features:\n")
cat("  • Project-local R library in renv/library\n")
cat("  • Lockfile created at renv.lock\n")
cat("  • Integration with conda R packages\n")
cat("\nUseful renv commands:\n")
cat("  renv::snapshot()  # Save current package state to renv.lock\n")
cat("  renv::restore()   # Restore packages from renv.lock\n")
cat("  renv::install()   # Install new packages\n")
cat("  renv::update()    # Update packages\n")
cat("  renv::status()    # Check environment status\n")
cat("\nThe project will automatically activate renv when you start R in this directory.\n")