# crobustascreen/flake.nix (Refactored)
{
  description = "Flake for CrobustaScreen using R, Python, and Julia";

  nixConfig = {
    bash-prompt = "\[CrobustaScreen$(__git_ps1 \" (%s)\")\]$ ";
  };

  inputs = {
    utils.url = "github:numtide/flake-utils";
    # Point all flakes to the same nixpkgs instance for consistency
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

    # --- Julia Package Sources ---
    # Non-flake git repos (use .src suffix for clarity)
    REPLVimSrc = { url = "github:kewiechecki/julia-repl-vim"; flake = false; };
    igraph_jllSrc = { url = "github:fcdimitr/igraph_jll.jl"; flake = false; };
    leiden_jllSrc = { url = "github:fcdimitr/leiden_jll.jl"; flake = false; };
    LeidenSrc = { url = "github:pitsianis/Leiden.jl"; flake = false; };

    # Flake dependencies providing Nix packages (like the new Autoencoders)
    Autoencoders = {
      url = "github:kewiechecki/Autoencoders.jl";
      flake = true;
      inputs.nixpkgs.follows = "nixpkgs"; # <<< Ensure consistency
    };
    # !!! IMPORTANT: Assume these are ALSO updated to provide packages.default !!!
    # !!! If not, change flake=false and handle like REPLVimSrc above     !!!
    DictMap = {
      url = "github:kewiechecki/DictMap.jl";
      flake = true; # Assumed updated
      inputs.nixpkgs.follows = "nixpkgs";
    };
    TrainingIO = {
      url = "github:kewiechecki/TrainingIO.jl";
      flake = true; # Assumed updated
      inputs.nixpkgs.follows = "nixpkgs";
    };
    DeePWAK = {
      url = "github:kewiechecki/DeePWAK.jl";
      flake = true; # Assumed updated
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs = { 
    self, nixpkgs, utils, 
	REPLVimSrc,
    igraph_jllSrc, leiden_jllSrc, LeidenSrc,
    Autoencoders, DictMap, TrainingIO, DeePWAK
  }@inputs:
    utils.lib.eachDefaultSystem (system:
      let
        # --- Julia Package Overlay ---
        # Defines how Nix finds/builds all custom Julia deps
		/*
        juliaOverlay = final: prev: {
          juliaPackages = prev.juliaPackages // {
            # Build from source inputs using buildJuliaPackage
            REPLVim = final.juliaPackages.buildJuliaPackage {
              pname = "REPLVim"; version = "git"; src = inputs.REPLVimSrc;
            };
            igraph_jll = final.juliaPackages.buildJuliaPackage {
              pname = "igraph_jll"; version = "git"; src = inputs.igraph_jllSrc;
              # JLLs might need special care, check JuliaNix examples
            };
            leiden_jll = final.juliaPackages.buildJuliaPackage {
              pname = "leiden_jll"; version = "git"; src = inputs.leiden_jllSrc;
              # JLLs might need special care
            };
            Leiden = final.juliaPackages.buildJuliaPackage {
              pname = "Leiden"; version = "git"; src = inputs.LeidenSrc;
            };

            # Import directly from flake dependency outputs
            Autoencoders = inputs.Autoencoders.packages.${system}.default;
            # --- ASSUMPTIONS ---
            DictMap = inputs.DictMap.packages.${system}.default;
            TrainingIO = inputs.TrainingIO.packages.${system}.default;
            DeePWAK = inputs.DeePWAK.packages.${system}.default;
            # --- If assumptions are wrong, build from src like REPLVim ---
            # e.g., DictMap = final.juliaPackages.buildJuliaPackage { pname="DictMap"; ... src=inputs.DictMap; };
          };
        };

        # --- Base Nixpkgs with Overlay ---
        pkgs = import nixpkgs {
          inherit system;
          overlays = [ juliaOverlay ];
          config = { allowUnfree = true; cudaSupport = system == "x86_64-linux"; };
        };

        # --- Julia Environment ---
        # Builds an environment based on ./Project.toml in CrobustaScreen repo
        # Assumes Project.toml lists deps like "REPLVim", "Autoencoders", "cuDNN", "PyCall" etc.
        # The overlay ensures these names map to the Nix derivations.
        juliaEnv = pkgs.juliaPackages.buildJuliaApplication {
          name = "crobusta-screen-julia-env";
          src = ./.; # Needs Project.toml + optionally Manifest.toml here

          # You might need to propagate runtime libs here if julia packages need them
          # and the wrappers don't handle it automatically.
          # propagatedBuildInputs = [ pkgs.cudaPackages.cudnn pkgs.stdenv.cc.cc.lib ];
        };
*/
        pkgs = import nixpkgs {
          inherit system;
          config = { allowUnfree = true; cudaSupport = system == "x86_64-linux"; };
        };

        autoencodersPkg = Autoencoders.packages.${system}.default;
        dictMapPkg = DictMap.packages.${system}.default;
        trainingIOPkg = TrainingIO.packages.${system}.default;

        # --- Python Environment ---
        pythonEnv = pkgs.python3.withPackages (ps: with ps; [ umap-learn leidenalg igraph ]);

        # --- Shell Packages ---
        # List everything needed in the final shell
        shellPkgsNested = with pkgs; [
          # Language Environments
          #juliaEnv  # The derivation containing Julia and ALL specified packages
		  julia
          R         # Base R
          pythonEnv # Python with its packages
          # R Packages (consider managing via renv/Nix integration if complex)
          rPackages.optparse rPackages.purrr rPackages.biomaRt rPackages.STRINGdb
          rPackages.class rPackages.cluster rPackages.fgsea rPackages.igraph rPackages.leiden
          rPackages.circlize rPackages.ComplexHeatmap rPackages.ggplot2 rPackages.ggpubr rPackages.umap
          # System Tools & Runtime Libs
          git
          stdenv.cc.cc.lib # Core GCC libs (for R, Python, Julia deps)
          gfortran         # Needed if anything still links gfortran libs directly
          (lib.optional stdenv.isLinux cudaPackages.cudatoolkit) # Runtime CUDA
          (lib.optional stdenv.isLinux cudaPackages.cudnn)      # Runtime cuDNN
          python3Packages.virtualenv # Maybe remove if pythonEnv is sufficient?
        ];
        shellPkgs = pkgs.lib.flatten shellPkgsNested;

        # --- Julia setup script for shellHook ---
        juliaSetupScript = ''
          using Pkg
          # Use a temporary environment to avoid polluting user's default
          Pkg.activate(; temp=true)
          # Ensure current project is activated inside the temp env logic??
          # This part is tricky - maybe activate "." is better? Let's try activating "."
          Pkg.activate(".")

          let # Define local mapping from name to Nix store path
              dev_pkgs = [
                  ("igraph_jll", "${inputs.igraph_jllSrc}"),
                  ("leiden_jll", "${inputs.leiden_jllSrc}"),
                  ("Leiden", "${inputs.LeidenSrc}")
              ]
          end

          println("+++ Ensuring development packages (non-flake git repos) are linked +++")
          # Force develop for non-flake sources every time; Pkg handles idempotency
          for (pkg_name, pkg_path) in dev_pkgs
              println("Developing ", pkg_name, " from ", pkg_path)
              # Add try-catch for robustness
              try Pkg.develop(path=pkg_path) catch e; println("WARN: Pkg.develop failed for ", pkg_name, ": ", e) end
          end

          println("+++ Instantiating project environment (resolving registry pkgs, linking flake pkgs?) +++")
          try Pkg.instantiate() catch e; println("WARN: Pkg.instantiate failed: ", e) end

          # Optional: Precompile everything now to avoid delay later
          # println("+++ Precompiling project +++")
          # try Pkg.precompile() catch e; println("WARN: Pkg.precompile failed: ", e) end

          println("+++ Julia setup script finished +++")
        '';

      in {
        # defaultPackage = ??? # Define if CrobustaScreen itself builds something installable

        devShell = pkgs.mkShell {
          name = "crobusta-screen-shell";
          # Use the combined list of environments and tools
          buildInputs = shellPkgs;

          # Dramatically simplified shellHook
          shellHook = ''
            source ${pkgs.git}/share/bash-completion/completions/git-prompt.sh
            # Set env vars needed by tools/packages at runtime
            export DEVICE="cuda:0"
            export R_HOME="${pkgs.R}/lib/R" # For RCall.jl
            export JULIA_PROJECT="@."       # Tell Julia to use the project env

            # Set LD_LIBRARY_PATH from everything in buildInputs
            # Should include paths from juliaEnv, R libs, Python libs, CUDA, GCC libs etc.
            export LD_LIBRARY_PATH="${pkgs.lib.makeLibraryPath shellPkgs}";

            echo "CrobustaScreen dev shell activated. Nix handles R/Py/Jl environments."
            # No more julia_deps.jl !
          '';
        };
      }
    );
}
