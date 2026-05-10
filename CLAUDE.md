# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

CrobustaScreen is a self-supervised phenotype-detection pipeline for confocal images of *Ciona robusta* embryos. It feeds Imaris-derived segmentation statistics through an autoencoder, builds a kNN graph in the embedding, partitions it with Leiden, and validates clusters against a known protein–protein interaction network assembled from STRINGdb orthologs (mouse + human → C. robusta).

The pipeline is polyglot by design: **R** does data ingestion, interaction lookups, clustering, and plotting; **Julia** does autoencoder training and the DEWAK/DEPWAK models; **Python** is invoked from R/Julia for `umap`, `igraph`, and `leidenalg`.

## Environment

Use the Nix flake — it provisions R, Julia, Python (with `umap-learn`, `leidenalg`, `igraph`), CUDA/cuDNN, and the system libs that R packages link against:

```bash
nix develop .
```

The `shellHook` sets `R_HOME` (for `RCall.jl`), `JULIA_PROJECT=@.`, `DEVICE=cuda:0`, and an `LD_LIBRARY_PATH` that includes `libpng`, `icu75`, `bzip2`, `curl`, `libxml2`, `gsl`, `openssl`. If R or Julia packages fail to load with linker errors, suspect a missing entry there before anything else.

The custom Julia packages live in sibling repos and are wired in through flake inputs: `Autoencoders.jl`, `TrainingIO.jl`, `DictMap.jl`, `DeePWAK.jl` (flake-based), plus `igraph_jll`, `leiden_jll`, `Leiden.jl`, `julia-repl-vim` (non-flake git sources). The big commented `juliaOverlay` block in `flake.nix` is the in-progress "real" Nix-managed Julia env; the active path just exposes plain `julia` and relies on `Pkg` to resolve the project.

## Common commands

Preprocessing chain (Makefile target `all` builds everything from upstream sources):

```bash
make all                       # builds data/interactions.csv and data/X.csv
Rscript readEmbryos.R          # parse segdat/ → data/embryodat.csv
Rscript readPheno.R            # parse imaris.csv → data/params.csv
julia preprocess.jl            # data/z_dat.csv → data/X.csv (z-score, scale to [-1,1])
```

Interaction network (only needed to regenerate `data/interactions.csv`):

```bash
Rscript cint.ensembl.R         # ENSEMBL ortholog lookup (mouse, human)
Rscript STRINGdb.R             # download STRINGdb interactions
Rscript get.interactions.R     # merge across species
Rscript interactionGraph.R     # build the interaction graph
```

Training, clustering, plotting:

```bash
julia autoencoder.jl --path "data/"                     # writes E.csv (and SAE/E.csv)
Rscript cluster.R                                       # writes data/k.csv, data/leiden.csv
Rscript plot.clust.R                                    # all figures under the dated out dir

# Cluster on PCs instead of autoencoder embeddings:
Rscript cluster.R   -e data/PCs.csv -o data/PCA
Rscript plot.clust.R -e data/PCs.csv -o fig/PCA/ -c data/PCA/ -s combined_score
```

DEWAK / DEPWAK (require `DeePWAK.jl`):

```bash
julia dewak.jl       # PCA + autoencoder + SAE DEWAK passes, writes data/DEWAK/MSE/...
julia dewakES.jl     # second DEWAK pass optimizing NES instead of MSE → data/DEWAK/NES/
julia cluster.jl     # DEPWAK clustering on top of both DEWAK runs → data/DEPWAK/...
```

`autoencoder.jl` flags: `--path/-p` (output dir, defaults to today's date), `--savecheckpts/-c` (write a checkpoint each epoch).

`cluster.R` flags: `-k/--k_min` (default 3), `-K/--k_max` (default 53), `-G/--gamma_max` (default 3.0), `-l/--leiden_reps` (default 1000), plus `data.parser` options `-e` (embedding csv) and `-o` (output dir).

`plot.clust.R` flags: `-c/--clust_dir` (clustering output dir, default `data`), `-s/--clust_sel_method` (`combined_score`, `ES`, `log2error`, `mean_silhouette`, or `nclusts`).

## Architecture

### Data flow

1. **Imaris segmentations** in `segdat/` → `readEmbryos.R` → `data/embryodat.csv`. Per-cell statistics summarized to per-embryo features (min/max/mean/sd for TVCs, ATMs, and inter-cell distances/cosines). See `desc.txt` for the full parameter dictionary — TVC = Trunk Ventral Cells, ATM = Anterior Tail Muscles, suffixes `_Nucleus`/`_Cell`.
2. `readPheno.R` parses `imaris.csv` → `data/params.csv` (114 embryo-level parameters).
3. `preprocess.jl` z-scores then rescales to [-1, 1] → `data/X.csv` (the matrix the autoencoder sees).
4. `autoencoder.jl` trains a 5-hidden-layer MLP `m → 58 → 29 → 14 → 29 → 58 → m` with `tanh`, AdamW (η=1e-4, λ=1e-4), 10000 epochs, 10% test split. Writes `data/E.csv` (encoder output) and a sparse-autoencoder variant under `data/SAE/E.csv` (architecture: `Dense(m→m,tanh) → SAE(m, 5m, relu) → Dense(m→m,tanh)`, with L1 sparsity penalty α=1e-5).
5. `cluster.R` builds a kNN graph on `E.csv`, scans `k ∈ [k_min, k_max]` selecting by GSEA enrichment of known interactions (`get.k`), then samples `gamma` uniformly in `[0.05, gamma_max]` and runs Leiden `leiden_reps` times via the `leiden`/`leidenalg` bridge → `data/k.csv`, `data/leiden.csv`.
6. `plot.clust.R` reads cluster output back via `read.clusts`, computes UMAP coords, runs hypergeometric / t / u tests for cluster–condition and cluster–phenotype enrichment, and emits the figure tree described in `DESCRIPTION.md`.

### Hyperparameter selection metrics

The clustering selection logic in `R/optimization.R` is the conceptual core of the pipeline; four metrics are computed per (k, γ):

- **GSEA enrichment score** (`fgsea`) — interactions ranked by edge count between condition pairs.
- **Recall** — fraction of known interactions present in any cluster, after building a per-condition graph from the partial-modularity score `H_xy = e_xy − γ K_x K_y / (2M)`.
- **log2 error** — reduced kNN classifier trained on a subset of embeddings with cluster labels, evaluated on the rest, averaged over 1000 reps.
- **Mean silhouette width** — over the embedding-space distance matrix.

Optima are picked by maximizing the product (`combined_score`) or whichever single metric `--clust_sel_method` selects.

### Module layout

- **R/** — library modules `source()`d by the top-level R scripts. The boundary between "script" (top-level `.R`) and "library" (`R/*.R`) is strict: top-level scripts only define CLI parsing and the call sequence. Notable: `optimization.R` (k/γ search), `leiden.R` (Python leiden bridge), `gene.network.R`, `hyper.R` (hypergeometric tests), `clustplots.R`, `pois.R`, `wak.R`, `io.R` (the `data.parser`, `parse.env`, `read.clusts`, `dir.csv` helpers).
- **Top-level Julia** — `autoencoder.jl` (training entry point), `dewak.jl` / `dewakES.jl` / `cluster.jl` (DEWAK/DEPWAK pipelines), `clustfns.jl`, `dewak.jl`, `gsea.jl` (FFI wrapper around the R `fgsea` package via `RCall`), `readDEWAKloss.jl` (the `@readDEWAK` macro that loads `loss.csv`/`PCs.csv` and unpacks `d_pca`, `k_pca`, `d_pcaNES`, ... into the caller scope), `plotloss.jl`, `plotfns.jl`, `interactions.jl`, `preprocess.jl`.
- **data/** — checked-in intermediate state. Subdirs `DEWAK/`, `DEPWAK/`, `PCA/`, `SAE/`, `fgsea/` are created by the corresponding Julia/R passes. `interactions.csv`, `G_STRINGdb.csv`, `stringdblabels.csv` are the regulatory-network artifacts; the pipeline reuses these unless `make` rebuilds them.
- **R/** scripts assume a `data.parser("data")` working directory unless `-o` is passed; output dirs are created lazily by `dir.csv` / `dir.plot`.

### Cross-language bridges

- Julia → R: `RCall.jl` is in `Project.toml`; `gsea.jl` calls `fgsea` directly. `R_HOME` must point at the Nix R for this to load.
- Julia → Python: `PyCall.jl`, used in `cluster.jl` via `pyleiden` / `pygraph` to drive `leidenalg`.
- R → Python: the `leiden`/`umap`/`leidenalg` R packages spawn Python under the hood; the flake's `pythonEnv` is what they pick up.

### `.gitignore` conventions

The pipeline writes a *lot* of figures. Almost everything image-shaped (`*.pdf`, `*.svg`, `*.eps`, `*.dot`) is gitignored, as are dated output dirs (`20*-*-*/**`, `20*_*_*/**`), `out/**`, `fig/**`, `presentation*/**`, `paper/**`, `notes*/**`, and any path containing whitespace or shell metacharacters. When adding output, prefer one of these conventions instead of inventing a new top-level dir, otherwise it will get accidentally committed.
