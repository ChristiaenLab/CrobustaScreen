{
  description = "Flake for CrobustaScreen using R, Python, and Julia";

  nixConfig = {
    bash-prompt = "\[CrobustaScreen$(__git_ps1 \" (%s)\")\]$ ";
  };
  inputs = {
    utils.url = "github:numtide/flake-utils";
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs = { self, nixpkgs, utils  }:
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
          ];
          shellHook = ''
            source ${git}/share/bash-completion/completions/git-prompt.sh
	        export LD_LIBRARY_PATH=/usr/lib
	        export DEVICE="cuda:0"
          '';
        };
      });
}
