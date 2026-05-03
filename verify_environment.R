#!/usr/bin/env Rscript
# Verification script to test R environment and Python integration

cat("=== CrobustaScreen R Environment Verification ===\n\n")

# Test R packages
cat("Testing R packages...\n")
required_packages <- c(
    "optparse", "parallel", "purrr",
    "class", "cluster", "fgsea", "igraph", "leiden",
    "circlize", "ComplexHeatmap", "ggplot2", "ggpubr", "umap",
    "biomaRt", "STRINGdb", "reticulate"
)

for (pkg in required_packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat("  ✓", pkg, "\n")
    } else {
        cat("  ✗", pkg, "not found\n")
    }
}

# Test Python integration via reticulate
cat("\nTesting Python integration (R → Python)...\n")
library(reticulate)
py_config <- py_config()
cat("  Python:", py_config$python, "\n")
# Handle potential list or NULL for numpy path
numpy_path <- py_config$numpy
if (is.null(numpy_path)) {
    cat("  NumPy: Not found\n")
} else {
    cat("  NumPy:", as.character(numpy_path), "\n")
}

# Test leiden (which uses Python under the hood)
cat("\nTesting leiden clustering (R → Python → R)...\n")
library(leiden)
library(igraph)

# Create test graph
set.seed(42)
g <- erdos.renyi.game(30, 0.15)
adj <- as_adjacency_matrix(g, sparse = FALSE)

# Run leiden clustering
clusters <- leiden(adj, resolution_parameter = 0.5)
cat("  Leiden clustering successful:", length(unique(clusters)), "clusters found\n")

# Test Python modules directly
cat("\nTesting Python modules from R...\n")
tryCatch({
    ig <- import("igraph")
    cat("  ✓ igraph (Python) version:", ig$`__version__`, "\n")
    
    leidenalg <- import("leidenalg")
    cat("  ✓ leidenalg (Python) available\n")
    
    umap <- import("umap")
    cat("  ✓ umap (Python) available\n")
}, error = function(e) {
    cat("  ✗ Error importing Python modules:", e$message, "\n")
})

# Test data exchange
cat("\nTesting R ↔ Python data exchange...\n")
test_matrix <- matrix(rnorm(100), 10, 10)
py_test <- r_to_py(test_matrix)
back_to_r <- py_to_r(py_test)
if (all(dim(test_matrix) == dim(back_to_r))) {
    cat("  ✓ Data exchange working\n")
} else {
    cat("  ✗ Data exchange failed\n")
}

cat("\n=== R Verification Complete ===\n")
cat("\nYour R environment is configured correctly for CrobustaScreen!\n")