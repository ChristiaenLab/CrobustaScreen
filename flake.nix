{
  description = "Flake for CrobustaScreen using R, Python, and Julia";

  nixConfig = {
    bash-prompt = "\[CrobustaScreen$(__git_ps1 \" (%s)\")\]$ ";
  };

  inputs = {
    utils.url = "github:numtide/flake-utils";
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

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
      flake = false;
    };
    DictMap = {
      url = "github:kewiechecki/DictMap.jl";
      flake = false;
    };
    TrainingIO = {
      url = "github:kewiechecki/TrainingIO.jl";
      flake = false;
    };
    DeePWAK = {
      url = "github:kewiechecki/DeePWAK.jl";
      flake = false;
    };
  };

  outputs = { self, nixpkgs, utils, igraph_jll, leiden_jll, Leiden,
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

        # Define the multi-line Julia script.
        # NOTE: The closing delimiter (two single quotes) MUST be flush with the left margin.
        juliaScript = ''
using Pkg
Pkg.add("cuDNN")
for (pkg, path) in [
    ("igraph_jll", "__IGRAPH_JLL__"),
    ("leiden_jll", "__LEIDEN_JLL__"),
    ("Leiden", "__LEIDEN__"),
    ("Autoencoders", "__AUTOENCODERS__"),
    ("TrainingIO", "__TRAININGIO__"),
    ("DictMap", "__DICTMAP__"),
    ("DeePWAK", "__DEE_PWAK__")
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
            exit(1)
        end
    end
end
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
            python3
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
sed -i \"s|__IGRAPH_JLL__|${toString igraph_jll}|g\" julia_deps.jl
sed -i \"s|__LEIDEN_JLL__|${toString leiden_jll}|g\" julia_deps.jl
sed -i \"s|__LEIDEN__|${toString Leiden}|g\" julia_deps.jl
sed -i \"s|__AUTOENCODERS__|${toString Autoencoders}|g\" julia_deps.jl
sed -i \"s|__TRAININGIO__|${toString TrainingIO}|g\" julia_deps.jl
sed -i \"s|__DICTMAP__|${toString DictMap}|g\" julia_deps.jl
sed -i \"s|__DEE_PWAK__|${toString DeePWAK}|g\" julia_deps.jl

# Replace the dollar placeholder with a literal dollar sign.
sed -i \"s|__DOLLAR_PLACEHOLDER__|\\\$|g\" julia_deps.jl

# Run the Julia script with LD_LIBRARY_PATH set to include libraries from gfortran and the stdenv compiler.
env LD_LIBRARY_PATH=${gfortranPath}/lib:${gccPath}/lib:${gccPath}/lib64:/usr/lib julia --project=. julia_deps.jl
";
        };
      });
}
