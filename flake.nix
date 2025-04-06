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

  outputs = { 
    self, nixpkgs, utils, 
	igraph_jll, leiden_jll, Leiden,
	Autoencoders, TrainingIO, DictMap, DeePWAK 
  }:
    utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
          config.allowUnfree = true;
          config.cudaSupport = system == "x86_64-linux";
        };
        ghcWithHasktorch = pkgs.haskellPackages.ghcWithPackages (pkgs: with pkgs; [
          hasktorch
          haskell-language-server
        ]);
      in {
        # A default package derivation that simply points to the repo source.
        defaultPackage = pkgs.stdenv.mkDerivation {
          name = "CrobustaScreen";
          src = ./.;
          # No build inputs if the project is interpreted
          buildInputs = [];
        };

        # A development shell with R, Python, and Julia available.
        devShell = with pkgs; mkShell {
          name = "crobusta-screen-shell";
          buildInputs = [
            R
            python3
            python3Packages.virtualenv
            julia
            git  # for git prompt support
            gcc
			gfortran
          ];
          shellHook = ''
            source ${git}/share/bash-completion/completions/git-prompt.sh
	        export DEVICE="cuda:0"

  			cat > julia_deps.jl <<EOF
			  using Pkg; 
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
                       @eval import \\$(Symbol(pkg))
                       println("Package ", pkg, " is already installed.")
                   catch e
                       println("Developing package ", pkg, " from ", path)
                       try
					       Pkg.develop(path=path)
					   catch e
					       println("Error precompiling ", pkg, ": ", e)
						   break
                   end
               end
			   EOF
            env LD_LIBRARY_PATH=${gfortran.libc}/lib:${gcc.libc}/lib:${gcc.libc}/lib64:/usr/lib \
			julia --project=. julia_deps.jl
          '';
        };
      });
}
