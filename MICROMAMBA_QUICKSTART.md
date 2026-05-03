# Micromamba Quick Start for CrobustaScreen

Micromamba is a lightweight, fast alternative to conda that's perfect for CI/CD and quick environment setup. It's 10-100x faster than conda for environment creation.

## Install Micromamba (2 minutes)

### Option 1: Quick Install (Linux/macOS)
```bash
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

### Option 2: Manual Install
```bash
# Download binary
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba

# Move to path
mkdir -p ~/.local/bin
mv bin/micromamba ~/.local/bin/

# Initialize shell
~/.local/bin/micromamba shell init -s bash -r ~/micromamba
source ~/.bashrc
```

### Option 3: Using Homebrew (macOS)
```bash
brew install micromamba
```

## Setup CrobustaScreen Environment (5 minutes)

```bash
# The setup script auto-detects micromamba
chmod +x setup_conda_env.sh
./setup_conda_env.sh

# Activate for future sessions
source ./activate_crobustascreen.sh
```

## Micromamba-Specific Commands

```bash
# Create environment directly (faster than script)
micromamba create -f environment.yml -y

# Activate
micromamba activate crobustascreen

# Install additional packages
micromamba install -c conda-forge some-package

# Update environment
micromamba update -f environment.yml

# List environments
micromamba env list

# Remove environment
micromamba env remove -n crobustascreen
```

## Why Micromamba?

| Feature | Micromamba | Conda | Mamba |
|---------|------------|-------|-------|
| Install size | ~5 MB | ~400 MB | ~100 MB |
| Env creation speed | Fastest | Slowest | Fast |
| Memory usage | Minimal | High | Medium |
| Solver | libmamba | conda | libmamba |
| Package manager | No | Yes | Yes |
| Drop-in replacement | Yes* | - | Yes |

*Micromamba can replace conda/mamba for most operations but doesn't include Python itself.

## Docker Integration

Micromamba is ideal for Docker containers:

```dockerfile
FROM mambaorg/micromamba:latest

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml && \
    micromamba clean --all --yes

# Set up Julia packages
USER root
RUN micromamba run -n base julia -e 'using Pkg; \
    ENV["PYTHON"]="/opt/conda/bin/python"; \
    Pkg.add("PyCall"); Pkg.build("PyCall"); \
    ENV["R_HOME"]="/opt/conda/lib/R"; \
    Pkg.add("RCall"); Pkg.build("RCall")'

USER $MAMBA_USER
WORKDIR /app
```

## CI/CD Example (GitHub Actions)

```yaml
name: Test
on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      
      - name: Setup Micromamba
        uses: mamba-org/setup-micromamba@v1
        with:
          environment-file: environment.yml
          cache-environment: true
          
      - name: Run tests
        shell: micromamba-shell {0}
        run: |
          julia verify_environment.jl
          Rscript verify_environment.R
```

## Troubleshooting

### Environment not found
```bash
# Check where micromamba stores environments
echo $MAMBA_ROOT_PREFIX
# Usually ~/micromamba or ~/.micromamba
```

### Julia can't find Python/R
```bash
# In Julia, after activating environment:
ENV["PYTHON"] = ENV["MAMBA_ROOT_PREFIX"] * "/envs/crobustascreen/bin/python"
ENV["R_HOME"] = ENV["MAMBA_ROOT_PREFIX"] * "/envs/crobustascreen/lib/R"
using Pkg
Pkg.build("PyCall")
Pkg.build("RCall")
```

### Activation issues
```bash
# Re-initialize shell integration
micromamba shell init -s bash
source ~/.bashrc
```

## Performance Tips

1. **Use libmamba solver** (default in micromamba)
2. **Pin major versions** in environment.yml to speed up solving
3. **Use `--strict-channel-priority`** to reduce solver complexity
4. **Cache environments** in CI/CD workflows
5. **Prefer conda-forge channel** for consistency

## Complete Setup Time

- Micromamba install: 1 minute
- Environment creation: 3-5 minutes (vs 15-30 with conda)
- Julia package config: 2 minutes
- Total: **~7 minutes**