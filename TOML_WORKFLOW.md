# TOML-Based R Dependency Management

## Rproject.toml Structure

Your dependencies are now managed in `Rproject.toml`, similar to Python's `pyproject.toml`:

```toml
[project]
name = "CrobustaScreen"
description = "Pipeline for self-supervised phenotype detection"
version = "1.0.0"
r-version = ">=4.3.0"

# Core dependencies (like Python's [project.dependencies])
[dependencies]
ggplot2 = "*"                    # Latest version
dplyr = "1.1.0"                  # Specific version
igraph = "*"
leiden = "*"

# Bioconductor packages
ComplexHeatmap = { source = "bioconductor" }
biomaRt = { source = "bioconductor" }

# Development dependencies (like Python's [project.optional-dependencies])
[dev-dependencies]
testthat = "*"
devtools = "*"
lintr = "*"

# Optional dependencies for enhanced features
[optional-dependencies]
plotly = "*"
gganimate = "*"

# GitHub packages
[github-dependencies]
"username/repo" = "main"
```

## Adding Dependencies

### Command Line (Recommended)

```bash
# Add to main dependencies
./renv_conda_integration.sh install ggplot2
./renv_conda_integration.sh install dplyr 1.1.0

# Add development dependency
Rscript renv_toml.R add testthat "*" dev-dependencies

# Add Bioconductor package
Rscript renv_toml.R add DESeq2 "*" dependencies bioconductor

# Add GitHub package
echo 'my-package = { source = "github", repo = "user/repo" }' >> Rproject.toml
```

### Manual Editing

Edit `Rproject.toml` directly and then sync:

```bash
# After editing Rproject.toml
./renv_conda_integration.sh restore
```

## Dependency Groups

### Main Dependencies
```toml
[dependencies]
# Packages required for the core functionality
ggplot2 = "*"
dplyr = "*"
```

### Development Dependencies
```toml
[dev-dependencies]
# Packages for development, testing, linting
testthat = "*"
devtools = "*"
lintr = "*"
```

### Optional Dependencies
```toml
[optional-dependencies]
# Packages for enhanced features
plotly = "*"
shiny = "*"
```

## Version Specifications

```toml
[dependencies]
# Any version
package1 = "*"

# Exact version
package2 = "1.0.0"

# Version range (if supported)
package3 = ">=1.0.0"

# Complex specification
package4 = { version = "1.0.0", source = "bioconductor" }

# GitHub package
package5 = { source = "github", repo = "user/repo", ref = "main" }
```

## Common Workflows

### Daily Development
```bash
# 1. Pull latest changes
git pull

# 2. Sync dependencies
./renv_conda_integration.sh restore

# 3. Start coding
R
```

### Adding New Features
```bash
# 1. Add required packages
./renv_conda_integration.sh install new_package

# 2. Develop feature
# 3. Test with all dependencies
./renv_conda_integration.sh restore
```

### Team Collaboration
```bash
# 1. Commit Rproject.toml and renv.lock
git add Rproject.toml renv.lock
git commit -m "Add new dependencies"

# 2. Team members sync
git pull
./renv_conda_integration.sh restore
```

## Migration from Traditional R

### From install.packages()
```r
# Old way
install.packages(c("ggplot2", "dplyr"))

# New way - add to Rproject.toml:
[dependencies]
ggplot2 = "*"
dplyr = "*"

# Then sync
```

### From BiocManager
```r
# Old way
BiocManager::install("DESeq2")

# New way - add to Rproject.toml:
[dependencies]
DESeq2 = { source = "bioconductor" }
```

### From GitHub
```r
# Old way
devtools::install_github("user/repo")

# New way - add to Rproject.toml:
[github-dependencies]
"user/repo" = "main"
```

## Commands Reference

| Action | Command |
|--------|---------|
| Add dependency | `./renv_conda_integration.sh install package` |
| Add with version | `./renv_conda_integration.sh install package 1.0.0` |
| Remove dependency | `./renv_conda_integration.sh remove package` |
| List dependencies | `./renv_conda_integration.sh list` |
| Sync from TOML | `./renv_conda_integration.sh restore` |
| Update all | `./renv_conda_integration.sh update` |
| Save lockfile | `./renv_conda_integration.sh snapshot` |

## Advanced Usage

### Direct R Commands
```r
# Read TOML configuration
source("renv_toml.R")
config <- read_rproject_toml()

# Install specific groups
install_from_toml(c("dependencies", "dev-dependencies"))

# Add dependency programmatically
add_dependency("new_package", version = "1.0.0", group = "dependencies")

# List dependencies
list_dependencies()
```

### Custom Groups
```toml
[plotting-dependencies]
ggplot2 = "*"
plotly = "*"
gganimate = "*"

[analysis-dependencies]  
dplyr = "*"
tidyr = "*"
broom = "*"
```

```bash
# Install specific group
Rscript renv_toml.R install plotting-dependencies
```

## Benefits Over Traditional R

1. **Declarative**: Dependencies are explicitly declared
2. **Reproducible**: Exact versions in lockfile
3. **Collaborative**: Easy to share and sync
4. **Organized**: Separate dev, optional, and core dependencies
5. **Familiar**: Same workflow as Python/Rust/Node.js
6. **Version Control**: Track dependency changes in git