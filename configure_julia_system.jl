#!/usr/bin/env julia
# Configure Julia to use conda Python and R (when using system Julia)

println("=== Configuring System Julia with Conda Environment ===")

# Check if in conda environment
if !haskey(ENV, "CONDA_PREFIX")
    error("Please activate the crobustascreen conda environment first!")
end

using Pkg

# Add required packages
println("Adding Julia packages...")
packages = ["PyCall", "RCall", "Conda"]
for pkg in packages
    if !haskey(Pkg.project().dependencies, pkg)
        println("  Adding $pkg...")
        Pkg.add(pkg)
    end
end

# Configure PyCall to use conda Python
println("\nConfiguring PyCall to use conda Python...")
ENV["PYTHON"] = ENV["CONDA_PREFIX"] * "/bin/python"
Pkg.build("PyCall")

# Configure RCall to use conda R  
println("\nConfiguring RCall to use conda R...")
ENV["R_HOME"] = ENV["CONDA_PREFIX"] * "/lib/R"
Pkg.build("RCall")

# Verify
println("\nVerifying configuration:")
using PyCall
using RCall

println("✓ PyCall Python: ", PyCall.pyprogramname)
println("✓ RCall R: ", RCall.Rhome)

# Test Python packages
py"""
import igraph
import leidenalg
import umap
print("✓ Python packages loaded successfully")
"""

# Test R
R"""
cat("✓ R loaded successfully\n")
cat("  R version:", R.version$version.string, "\n")
"""

# Install other project dependencies
println("\nInstalling project dependencies...")
Pkg.instantiate()

println("\n=== Configuration Complete ===")
println("You can now use Julia with the conda Python and R environments!")