using Pkg
#Pkg.activate(".")
Pkg.add("cuDNN")
for (pkg, path) in [
    ("igraph_jll", "/nix/store/p89x11x3nb62b9qvd57rgyghpkvggnvl-source"),
    ("leiden_jll", "/nix/store/5hyzr4d5nj51mvii224n0z1dw14ywnir-source"),
    ("Leiden", "/nix/store/611wgj6ynw0knrk6jrnjdc4q935pl83j-source"),
    ("Autoencoders", "/nix/store/qm93fajgx1k0gxg6sfcya1kkwrf1ykis-source"),
    ("TrainingIO", "/nix/store/136777yvj47rrrxadiaxnianwx42r0md-source"),
    ("DictMap", "/nix/store/hsx3vah9niq2b2wqzccyxk1dyvvqxipb-source"),
    ("DeePWAK", "/nix/store/g4dhmj7xrxx3nbvhipp2v57j2dc9ihw0-source")
]
    try
        @eval import $(Symbol(pkg))
        println("Package ", pkg, " is already installed.")
    catch e
        println("Developing package ", pkg, " from ", path)
        try
            Pkg.develop(path=path)
            Pkg.precompile(only=[pkg])
        catch e
            println("Error precompiling ", pkg, ": ", e)
            #exit(1)
        end
    end
end
#Pkg.instantiate()
Pkg.precompile()

