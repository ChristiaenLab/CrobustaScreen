# renv + TOML Guide for CrobustaScreen

## Overview

This project uses `renv` with TOML-based dependency management, providing a `pyproject.toml`-like experience for R. Dependencies are declared in `Rproject.toml` and managed automatically, similar to Python's `uv` or Rust's Cargo.

## Quick Start

### Initial Setup (First Time)

```bash
# 1. Activate conda environment (if using)
source ./activate_crobustascreen.sh

# 2. Sync dependencies from Rproject.toml
./renv_conda_integration.sh restore

# 3. The project will automatically activate renv when you start R
```

### Daily Usage

```bash
# Activate conda environment (recommended)
source ./activate_crobustascreen.sh

# Sync dependencies (like uv sync or cargo build)
./renv_conda_integration.sh restore

# Start R - renv activates automatically
R
```

## Common Tasks

### TOML-Based Dependency Management

```bash
# List all dependencies from Rproject.toml
./renv_conda_integration.sh list

# Add a new dependency (like uv add or cargo add)
./renv_conda_integration.sh install ggplot2
./renv_conda_integration.sh install dplyr 1.1.0  # specific version

# Remove a dependency (like uv remove)
./renv_conda_integration.sh remove ggplot2

# Sync dependencies from Rproject.toml (like uv sync)
./renv_conda_integration.sh restore

# Check environment status
./renv_conda_integration.sh status

# Save current state to lockfile
./renv_conda_integration.sh snapshot
```

### Manual renv Commands in R

```r
# Install packages
renv::install("package_name")
renv::install(c("package1", "package2"))

# Install specific version
renv::install("package@1.0.0")

# Install from GitHub
renv::install("username/repo")

# Update packages
renv::update()              # Update all
renv::update("package")     # Update specific

# Save current state
renv::snapshot()

# Restore from lockfile
renv::restore()

# Check status
renv::status()

# Clean unused packages
renv::clean()
```

## File Structure

- `Rproject.toml` - Project configuration with dependencies (like pyproject.toml)
- `renv.lock` - Lockfile containing exact package versions (like uv.lock)
- `renv/` - Project-local R library
- `.Rprofile` - Auto-activates renv when R starts
- `renv_toml.R` - TOML parser and dependency manager
- `renv_conda_integration.sh` - CLI interface for common operations

## Integration with Conda

The setup is configured to work with conda's R installation while maintaining separate package management:

1. **Conda manages**: R binary, system libraries, Python packages
2. **renv manages**: R packages and their versions
3. **Both work together**: renv uses conda's R but installs packages locally

## Workflow Comparison

| Task | Python (uv) | R (TOML + renv) | Traditional R |
|------|-------------|-----------------|---------------|
| Config file | `pyproject.toml` | `Rproject.toml` | None |
| Add package | `uv add package` | `./renv_conda_integration.sh install package` | `install.packages()` |
| Sync from config | `uv sync` | `./renv_conda_integration.sh restore` | N/A |
| Save lockfile | `uv lock` | `./renv_conda_integration.sh snapshot` | N/A |
| Update packages | `uv update` | `./renv_conda_integration.sh update` | `update.packages()` |
| List packages | `uv list` | `./renv_conda_integration.sh list` | `installed.packages()` |
| Remove package | `uv remove package` | `./renv_conda_integration.sh remove package` | `remove.packages()` |

## Advantages

1. **Reproducibility**: Exact package versions recorded in `renv.lock`
2. **Isolation**: Project-specific library doesn't affect system R
3. **Collaboration**: Share `renv.lock` for identical environments
4. **Cache**: Packages cached locally for faster installs
5. **Integration**: Works with conda, Docker, CI/CD

## Troubleshooting

### Package Installation Issues

```r
# If a package fails to install, try:
options(repos = c(CRAN = "https://cloud.r-project.org"))
renv::install("package")

# For Bioconductor packages:
renv::install("bioc::package")
```

### Python Integration

```r
# Ensure reticulate uses conda Python
library(reticulate)
use_python(file.path(Sys.getenv("CONDA_PREFIX"), "bin", "python"))
```

### Reset renv

```r
# If things go wrong, reset:
renv::deactivate()
unlink("renv", recursive = TRUE)
unlink("renv.lock")
# Then reinitialize:
renv::init()
```

## Best Practices

1. **Always snapshot after installing**: `renv::snapshot()` after adding packages
2. **Commit renv.lock**: Always commit `renv.lock` to version control
3. **Don't commit renv/library**: Add `renv/library/` to `.gitignore`
4. **Document new packages**: Update this guide when adding major dependencies
5. **Test restore**: Periodically test `renv::restore()` works correctly

## Package Sources

renv can install from multiple sources:

- CRAN: `renv::install("package")`
- Bioconductor: `renv::install("bioc::package")`
- GitHub: `renv::install("user/repo")`
- GitLab: `renv::install("gitlab::user/repo")`
- Specific version: `renv::install("package@1.0.0")`
- Local: `renv::install("path/to/package")`

## For CI/CD

```yaml
# Example GitHub Actions
- name: Setup R
  uses: r-lib/actions/setup-r@v2
  
- name: Setup renv
  uses: r-lib/actions/setup-renv@v2
  
- name: Restore packages
  run: |
    Rscript -e 'renv::restore()'
```

## Additional Resources

- [renv documentation](https://rstudio.github.io/renv/)
- [renv cheatsheet](https://github.com/rstudio/cheatsheets/blob/master/renv.pdf)
- [Migration guide](https://rstudio.github.io/renv/articles/renv.html)