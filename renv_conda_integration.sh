#!/bin/bash
# Integration script for using renv with conda environment
# This ensures both conda and renv work together seamlessly

set -e

echo "=== CrobustaScreen renv + conda Integration ==="

# Check if conda environment is activated
if [ -z "$CONDA_PREFIX" ]; then
    echo "Warning: No conda environment detected."
    echo "Please activate the crobustascreen environment first:"
    echo "  source ./activate_crobustascreen.sh"
    echo ""
    echo "Or continue with system R (not recommended)"
    read -p "Continue anyway? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# Function to use conda's R if available
setup_r_environment() {
    if [ -n "$CONDA_PREFIX" ]; then
        export R_HOME=$CONDA_PREFIX/lib/R
        export PATH=$CONDA_PREFIX/bin:$PATH
        echo "Using conda R: $R_HOME"
    else
        echo "Using system R"
    fi
}

# Function to restore renv packages
restore_renv() {
    # Check if Rproject.toml exists and use TOML-based workflow
    if [ -f "Rproject.toml" ]; then
        echo "Syncing dependencies from Rproject.toml..."
        Rscript renv_toml.R sync
    else
        echo "Restoring R packages from renv.lock..."
        Rscript -e "
        if (!requireNamespace('renv', quietly = TRUE)) {
            install.packages('renv')
        }
        renv::restore(prompt = FALSE)
        cat('\nPackages restored successfully!\n')
        "
    fi
}

# Function to snapshot current state
snapshot_renv() {
    echo "Creating snapshot of current R packages..."
    Rscript -e "
    renv::snapshot(prompt = FALSE)
    cat('Snapshot saved to renv.lock\n')
    "
}

# Function to install a new package
install_package() {
    local package=$1
    local version=${2:-"*"}
    echo "Installing $package..."
    
    # Check if Rproject.toml exists and use TOML-based workflow
    if [ -f "Rproject.toml" ]; then
        echo "Using TOML-based dependency management..."
        Rscript renv_toml.R add "$package" "$version"
    else
        echo "Using direct renv installation..."
        Rscript -e "
        renv::install('$package')
        renv::snapshot(prompt = FALSE)
        cat('Package installed and lockfile updated\n')
        "
    fi
}

# Function to update packages
update_packages() {
    echo "Updating R packages..."
    Rscript -e "
    renv::update(prompt = FALSE)
    renv::snapshot(prompt = FALSE)
    cat('Packages updated and lockfile saved\n')
    "
}

# Function to check environment status
check_status() {
    echo "Checking renv status..."
    Rscript -e "
    renv::status()
    cat('\n')
    cat('Python config:\n')
    library(reticulate)
    py_config()
    "
}

# Main menu
case "${1:-}" in
    restore)
        setup_r_environment
        restore_renv
        ;;
    snapshot)
        setup_r_environment
        snapshot_renv
        ;;
    install)
        if [ -z "$2" ]; then
            echo "Usage: $0 install <package_name> [version]"
            exit 1
        fi
        setup_r_environment
        install_package "$2" "$3"
        ;;
    remove)
        if [ -z "$2" ]; then
            echo "Usage: $0 remove <package_name>"
            exit 1
        fi
        setup_r_environment
        if [ -f "Rproject.toml" ]; then
            Rscript renv_toml.R remove "$2"
        else
            echo "Rproject.toml not found. Cannot remove package."
            exit 1
        fi
        ;;
    list)
        if [ -f "Rproject.toml" ]; then
            Rscript renv_toml.R list "$2"
        else
            echo "Rproject.toml not found."
            exit 1
        fi
        ;;
    update)
        setup_r_environment
        update_packages
        ;;
    status)
        setup_r_environment
        check_status
        ;;
    *)
        echo "Usage: $0 {restore|snapshot|install|update|status|list|remove}"
        echo ""
        echo "Commands:"
        echo "  restore         - Sync from Rproject.toml or restore from renv.lock"
        echo "  snapshot        - Save current package state to renv.lock"
        echo "  install <pkg>   - Add package to Rproject.toml and install"
        echo "  remove <pkg>    - Remove package from Rproject.toml"
        echo "  update          - Update all packages and save state"
        echo "  status          - Check renv and Python integration status"
        echo "  list            - List dependencies from Rproject.toml"
        echo ""
        echo "Examples:"
        echo "  $0 restore              # Sync from Rproject.toml"
        echo "  $0 install ggplot2     # Add ggplot2 to dependencies"
        echo "  $0 install dplyr 1.1.0 # Add specific version"
        echo "  $0 remove ggplot2      # Remove ggplot2"
        echo "  $0 list                 # Show all dependencies"
        exit 1
        ;;
esac