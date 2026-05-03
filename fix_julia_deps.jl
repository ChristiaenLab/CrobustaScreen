#!/usr/bin/env julia
# Fix Julia dependencies and configuration issues

using Pkg

println("=== Fixing Julia Dependencies ===")

# 1. Remove problematic CUDA packages if they exist
println("1. Removing problematic CUDA packages...")
try
    Pkg.rm(["CUDA", "cuDNN", "Atomix", "CUDNN_jll"], mode=Pkg.PKGMODE_MANIFEST)
    println("   ✓ Removed CUDA packages")
catch e
    println("   ℹ CUDA packages not found or already removed")
end

# 2. Configure Python for PyCall
println("2. Configuring Python for PyCall...")
if haskey(ENV, "CONDA_PREFIX") && isfile(joinpath(ENV["CONDA_PREFIX"], "bin", "python"))
    ENV["PYTHON"] = joinpath(ENV["CONDA_PREFIX"], "bin", "python")
    println("   ✓ Using conda Python: $(ENV["PYTHON"])")
elseif isfile("/usr/bin/python3")
    ENV["PYTHON"] = "/usr/bin/python3"
    println("   ✓ Using system Python: $(ENV["PYTHON"])")
else
    ENV["PYTHON"] = ""
    println("   ℹ Using Julia's built-in conda")
end

# 3. Rebuild PyCall and Conda
println("3. Rebuilding PyCall and Conda...")
try
    Pkg.build(["Conda", "PyCall"])
    println("   ✓ PyCall and Conda rebuilt successfully")
catch e
    println("   ⚠ Error rebuilding PyCall/Conda: $e")
    println("   Trying alternative approach...")
    
    # Alternative: Force reinstall
    try
        Pkg.rm(["PyCall", "Conda"])
        Pkg.add(["PyCall", "Conda"])
        Pkg.build(["PyCall", "Conda"])
        println("   ✓ PyCall and Conda reinstalled successfully")
    catch e2
        println("   ✗ Failed to fix PyCall: $e2")
    end
end

# 4. Clean up and instantiate
println("4. Cleaning up and instantiating project...")
try
    Pkg.gc()
    Pkg.instantiate()
    println("   ✓ Project instantiated successfully")
catch e
    println("   ⚠ Warning during instantiation: $e")
end

# 5. Test key packages
println("5. Testing key package imports...")
test_packages = ["ArgParse", "JLD2", "PlotUtils", "Latexify", "DataFrames"]

for pkg in test_packages
    try
        eval(Meta.parse("using $pkg"))
        println("   ✓ $pkg loads successfully")
    catch e
        println("   ⚠ $pkg failed to load: $e")
    end
end

# 6. Test Python integration
println("6. Testing Python integration...")
try
    using PyCall
    pyimport("sys")
    println("   ✓ PyCall working with Python $(PyCall.pyprogramname)")
catch e
    println("   ⚠ PyCall not working: $e")
end

println("\n=== Julia Dependencies Fixed ===")
println("You can now run: julia autoencoder.jl")
println("If you still have issues, try: julia --project=. -e 'using Pkg; Pkg.precompile()'")