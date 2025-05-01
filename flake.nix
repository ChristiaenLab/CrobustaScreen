{
  description = "Flake for CrobustaScreen using R, Python, and Julia";

  nixConfig = {
    bash-prompt = "\[CrobustaScreen$(__git_ps1 \" (%s)\")\]$ ";
  };

  inputs = {
    utils.url = "github:numtide/flake-utils";
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

	REPLVim = {
	  url = "github:kewiechecki/julia-repl-vim";
	  flake = false;
	};
    igraph_jll = {
      url = "github:fcdimitr/igraph_jll.jl";
      flake = false;
    };
    leiden_jll = {
      url = "github:fcdimitr/leiden_jll.jl";
      flake = false;
    };
    Leiden = {
      url = "github:pitsianis/Leiden.jl";
      flake = false;
    };

    Autoencoders = {
      url = "github:kewiechecki/Autoencoders.jl";
      flake = true;
    };
    DictMap = {
      url = "github:kewiechecki/DictMap.jl";
      flake = true;
    };
    TrainingIO = {
      url = "github:kewiechecki/TrainingIO.jl";
      flake = true;
    };
    DeePWAK = {
      url = "github:kewiechecki/DeePWAK.jl";
      flake = true;
    };
  };

  outputs = { self, nixpkgs, utils, REPLVim, igraph_jll, leiden_jll, Leiden,
              Autoencoders, TrainingIO, DictMap, DeePWAK }:
    utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
          config.allowUnfree = true;
          config.cudaSupport = system == "x86_64-linux";
        };

        # Get library paths from the stdenv compiler and from gfortran.
        gccPath = toString pkgs.stdenv.cc.cc.lib;
        gfortranPath = toString pkgs.gfortran;
		rPath = toString pkgs.R;

        # Define the multi-line Julia script.
        # NOTE: The closing delimiter (two single quotes) MUST be flush with the left margin.
        juliaScript = ''
using Pkg
#Pkg.activate(".")
#Pkg.add(url="__REPLVim__")
Pkg.instantiate()

Pkg.add("cuDNN")
Pkg.add("PyCall")
Pkg.build("RCall")
Pkg.add("StructArrays")

for (pkg, path) in [
    ("REPLVim", "__REPLVim__"),
    ("igraph_jll", "__IGRAPH_JLL__"),
    ("leiden_jll", "__LEIDEN_JLL__"),
    ("Leiden", "__LEIDEN__"),
    ("Autoencoders", "__AUTOENCODERS__"),
    ("TrainingIO", "__TRAININGIO__"),
    ("DictMap", "__DICTMAP__"),
    ("DeePWAK", "__DEEPWAK__")
]
    try
        @eval import __DOLLAR_PLACEHOLDER__(Symbol(pkg))
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
Pkg.precompile()
using DeePWAK, PyCall, RCall, cuDNN
'';

      in {
        defaultPackage = pkgs.stdenv.mkDerivation {
          name = "CrobustaScreen";
          src = ./.;
          buildInputs = [];
        };

        devShell = with pkgs; mkShell {
          name = "crobusta-screen-shell";
          buildInputs = [
            R
            pkgs.rPackages.optparse
            pkgs.rPackages.purrr
			# to fetch interactions
            pkgs.rPackages.biomaRt
            pkgs.rPackages.STRINGdb
			# clustering
            pkgs.rPackages.class
            pkgs.rPackages.cluster
            pkgs.rPackages.fgsea
            pkgs.rPackages.igraph
            pkgs.rPackages.leiden
			# visualization
			pkgs.rPackages.circlize
			pkgs.rPackages.ComplexHeatmap
			pkgs.rPackages.ggplot2
			pkgs.rPackages.ggpubr
			pkgs.rPackages.umap
            # Python environment with selected packages.
            (python3.withPackages (ps: with ps; [ umap-learn leidenalg igraph ]))
            python3Packages.virtualenv
            julia
            git      # for git prompt support
            pkgs.stdenv.cc
            pkgs.gfortran
          ];
          shellHook = "
source ${git}/share/bash-completion/completions/git-prompt.sh
export DEVICE=\"cuda:0\"

# Write the Julia script to a file.
cat > julia_deps.jl <<'EOF'
${juliaScript}
EOF

# Replace placeholders with actual paths.
sed -i \"s|__REPLVim__|${toString REPLVim}|g\" julia_deps.jl
sed -i \"s|__IGRAPH_JLL__|${toString igraph_jll}|g\" julia_deps.jl
sed -i \"s|__LEIDEN_JLL__|${toString leiden_jll}|g\" julia_deps.jl
sed -i \"s|__LEIDEN__|${toString Leiden}|g\" julia_deps.jl
sed -i \"s|__AUTOENCODERS__|${toString Autoencoders}|g\" julia_deps.jl
sed -i \"s|__TRAININGIO__|${toString TrainingIO}|g\" julia_deps.jl
sed -i \"s|__DICTMAP__|${toString DictMap}|g\" julia_deps.jl
sed -i \"s|__DEEPWAK__|${toString DeePWAK}|g\" julia_deps.jl

# Replace the dollar placeholder with a literal dollar sign.
sed -i \"s|__DOLLAR_PLACEHOLDER__|\\\$|g\" julia_deps.jl

# Run the Julia script with LD_LIBRARY_PATH set to include libraries from gfortran and the stdenv compiler.
export R_HOME=${rPath}/lib/R
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${gfortranPath}/lib:${gccPath}/lib:${gccPath}/lib64:${rPath}/lib/R/lib
echo $LD_LIBRARY_PATH
julia --project=. julia_deps.jl
";
        };
      });
}
