#!/bin/bash
# Configure Julia and R after activating conda environment
# Run this AFTER activating the crobustascreen environment

echo "=== Configuring Julia and R Integration ==="

# Check if in correct environment
# This check is being removed as `micromamba run` ensures the correct environment context.
# if [[ "$CONDA_DEFAULT_ENV" != "crobustascreen" ]]; then
#     echo "ERROR: Please activate the crobustascreen environment first:"
#     echo "  conda activate crobustascreen"
#     echo "  OR"
#     echo "  micromamba activate crobustascreen"
#     exit 1
# fi

# Get paths
PYTHON_PATH=$CONDA_PREFIX/bin/python
R_HOME=$CONDA_PREFIX/lib/R

echo "Python: $PYTHON_PATH"
echo "R: $R_HOME"

# Configure Julia packages
echo ""
echo "Configuring Julia packages..."
julia --project=. << 'EOF'
using Pkg

# Add required packages if not already present
println("Checking Julia packages...")
packages = ["PyCall", "RCall", "Conda"]
for pkg in packages
    if !haskey(Pkg.project().dependencies, pkg)
        println("Adding $pkg...")
        Pkg.add(pkg)
    end
end

# Configure PyCall to use conda Python
println("\nConfiguring PyCall...")
ENV["PYTHON"] = ENV["CONDA_PREFIX"] * "/bin/python"
Pkg.build("PyCall")

# Configure RCall to use conda R
println("\nConfiguring RCall...")
ENV["R_HOME"] = ENV["CONDA_PREFIX"] * "/lib/R"
Pkg.build("RCall")

# Verify configurations
println("\nVerifying configuration:")
using PyCall
using RCall
println("✓ PyCall Python: ", PyCall.pyprogramname)
println("✓ RCall R: ", RCall.Rhome)

# Install other Julia dependencies
println("\nInstalling project dependencies...")
Pkg.instantiate()
println("✓ Julia configuration complete")
EOF

# Configure R packages
echo ""
echo "Configuring R packages..."
Rscript - << 'EOF'
# Check and install packages from CRAN
cat("Installing R packages if needed...\n")
cran_packages <- c("leiden", "umap")
for (pkg in cran_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat("Installing", pkg, "...\n")
        install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
    } else {
        cat("✓", pkg, "already installed\n")
    }
}

# Install STRINGdb from Bioconductor if needed
if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    cat("Installing STRINGdb from Bioconductor...\n")
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", quiet = TRUE)
    BiocManager::install("STRINGdb", update = FALSE, ask = FALSE)
} else {
    cat("✓ STRINGdb already installed\n")
}

# Configure reticulate to use conda Python
library(reticulate)
use_python(Sys.getenv("CONDA_PREFIX"), required = TRUE)
cat("✓ R configuration complete\n")
cat("  Reticulate Python:", py_config()$python, "\n")
EOF

echo ""
echo "=== Configuration Complete ==="
echo ""
echo "You can now run:"
echo "  julia verify_environment.jl"
echo "  Rscript verify_environment.R"
echo ""
echo "To start working:"
echo "  julia --project=."
echo "  R"