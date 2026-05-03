#!/usr/bin/env julia
# Verification script to test cross-language environment setup
# Tests PyCall (Julia→Python) and RCall (Julia→R) integration

using Pkg
Pkg.instantiate()

using PyCall
using RCall
try
    using CUDA
catch
    # CUDA might not be installed or configured
end

println("=== CrobustaScreen Environment Verification ===\n")

# Test PyCall
println("Testing PyCall (Julia → Python)...")
try
    println("✓ PyCall loaded")
    println("  Python executable: $(PyCall.pyprogramname)")
    
    # Test Python packages
    py"""
    import sys
    import igraph
    import leidenalg
    import umap
    print(f"  Python version: {sys.version}")
    print(f"  igraph version: {igraph.__version__}")
    print(f"  leidenalg version: {leidenalg.__version__}")
    print(f"  umap version: {umap.__version__}")
    """
    
    # Test actual leiden clustering
    py"""
    import igraph as ig
    import leidenalg
    # Create small test graph
    g = ig.Graph.Erdos_Renyi(n=10, p=0.3)
    partition = leidenalg.find_partition(g, leidenalg.ModularityVertexPartition)
    print(f"  Leiden test clustering: {len(set(partition.membership))} clusters found")
    """
    println("✓ Python packages working\n")
catch e
    println("✗ PyCall error: $e\n")
end

# Test RCall
println("Testing RCall (Julia → R)...")
try
    println("✓ RCall loaded")
    println("  R home: $(RCall.Rhome)")
    
    # Test R version and packages
    R"""
    cat("  R version:", R.version$version.string, "\n")
    
    # Test required packages
    packages <- c("leiden", "igraph", "reticulate", "cluster", "fgsea")
    for (pkg in packages) {
        if (require(pkg, character.only = TRUE, quietly = TRUE)) {
            cat("  ✓", pkg, "loaded\n")
        } else {
            cat("  ✗", pkg, "not found\n")
        }
    }
    """
    
    # Test R calling Python via reticulate
    R"""
    library(reticulate)
    py_config <- py_config()
    cat("  Reticulate Python:", py_config$python, "\n")
    
    # Test leiden through R
    library(leiden)
    library(igraph)
    g <- erdos.renyi.game(20, 0.2)
    adj <- as_adjacency_matrix(g, sparse = FALSE)
    clusters <- leiden(adj, resolution_parameter = 0.5)
    cat("  R leiden test: ", length(unique(clusters)), " clusters found\n")
    """
    println("✓ R packages working\n")
catch e
    println("✗ RCall error: $e\n")
end

# Test data flow between languages
println("Testing cross-language data flow...")
try
    # Julia → Python → R → Julia
    julia_data = rand(5, 3)
    println("  Julia matrix: $(size(julia_data))")
    
    # Send to Python
    py"""
    import numpy as np
    python_data = $julia_data
    python_processed = np.mean(python_data, axis=0)
    print(f"  Python processed: shape {python_processed.shape}")
    """
    python_result = py"python_processed"
    
    # Send to R
    @rput python_result
    R"""
    r_processed <- mean(python_result)
    cat("  R processed: mean =", r_processed, "\n")
    """
    @rget r_processed
    
    println("  Final result back in Julia: $r_processed")
    println("✓ Cross-language data flow working\n")
catch e
    println("✗ Data flow error: $e\n")
end

# Test CUDA if available
println("Testing CUDA availability...")
try
    # using CUDA (already loaded at top)
    if isdefined(Main, :CUDA) && CUDA.functional()
        println("✓ CUDA is functional")
        println("  GPU: $(CUDA.name(CUDA.device()))")
        println("  CUDA version: $(CUDA.version())")
    else
        println("○ CUDA not functional (CPU mode will be used)")
    end
catch e
    println("○ CUDA not available: $e")
    println("  (This is OK - the pipeline can run on CPU)")
end

println("\n=== Verification Complete ===")
println("\nIf all tests passed, your environment is ready for CrobustaScreen!")