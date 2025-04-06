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
            git        # for git prompt support
            gcc
            gfortran
          ];
          shellHook = ''
source ${git}/share/bash-completion/completions/git-prompt.sh
export DEVICE="cuda:0"

# Write the Julia script to a file named julia_deps.jl.
cat > julia_deps.jl <<EOF
using Pkg
for (pkg, path) in [
    ("igraph_jll", "${toString igraph_jll}"),
    ("leiden_jll", "${toString leiden_jll}"),
    ("Leiden", "${toString Leiden}"),
    ("Autoencoders", "${toString Autoencoders}"),
    ("TrainingIO", "${toString TrainingIO}"),
    ("DictMap", "${toString DictMap}"),
    ("DeePWAK", "${toString DeePWAK}")
]
    try
        @eval import __DOLLAR__(Symbol(pkg))
        println("Package ", pkg, " is already installed.")
    catch e
        println("Developing package ", pkg, " from ", path)
        try
            Pkg.develop(path=path)
            Pkg.precompile(only=[pkg])
        catch e
            println("Error precompiling ", pkg, ": ", e)
            #break
			exit(1)
        end
    end
end
EOF

# Replace the placeholder __DOLLAR__ with a literal dollar sign.
sed -i 's/__DOLLAR__/\$/g' julia_deps.jl

# Run the Julia script with a temporary LD_LIBRARY_PATH.
env LD_LIBRARY_PATH=${gfortran.libc}/lib:${gcc.libc}/lib:${gcc.libc}/lib64:/usr/lib \
  julia --project=. julia_deps.jl
'';
        };
      });
}
