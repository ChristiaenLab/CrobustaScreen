<script
  src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"
  type="text/javascript">
</script>

# CrobustaScreen
A pipeline for self-supervised phenotype detection from confocal images of *Ciona robusta* embryos.

![overview](https://github.com/ChristiaenLab/CrobustaScreen/blob/main/overview.png?raw=true)

## Dependencies

### R
`optparse`, `parallel`, `purrr`

To fetch interactions:
`biomaRt`, `STRINGdb`

Clustering:
`class`, `cluster`, `fgsea`, `igraph`, `leiden`, 

Visualization:
`circlize`, `ComplexHeatmap`, `ggplot2`, `ggpubr`, `umap` 

### Python 
`leidenalg`, `igraph`, `umap`

These should be automatically installed when loading the appropriate R packages.

### Julia 
[`Autoencoders.jl`](https://github.com/kewiechecki/Autoencoders.jl), [`TrainingIO.jl`](https://github.com/kewiechecki/TrainingIO.jl)

## Usage

### Internal prerequisites
All of the necessary files are included in `data/`. You should not need to run these; they are only included for completeness.

All missing files can be created with
```bash
make all
```

#### Preprocessing
The pipeline uses features extracted from segmentation of confocal images using Imaris. Summary statistics are extracted for segmened cells in each embryo. From the cell segmentation statistics 116 embryo-level parameters are computed. Parameters are normalized by z-score then scaled between -1 and 1.

![overview](https://github.com/ChristiaenLab/CrobustaScreen/blob/main/presentation/segmentation.png?raw=true)
Parse segmentation data from `segdat/`.
```bash
Rscript readEmbryos.R
```

Parse `imaris.csv`.
```bash
Rscript readPheno.R
```

Additional preprocessing to read z-scores into Julia.
```bash
julia preprocess.jl
```

#### Generate putative protein interactions
The final network combines known interactions from orthologous genes in *Ciona intestinalis*, mouse, and human.

Building the interaction table additionally requires `biomaRt` and `STRINGdb`.

##### Ortholog lookup
We use [STRINGdb](https://doi.org/10.1093/nar/gkq973) to construct a known protein interaction network of the perturbed genes. Because the *C. robusta* network is poorly characterized, we use ENSEMBL to obtain orthologs from *M. musculus* and *H. sapiens*. 

Fetch orthologs from ENSEMBL
```bash
Rscript cint.ensembl.R
```

##### Putative interaction network
Download interactions between orthologs from STRINGdb
```bash
Rscript STRINGdb.R
```

Merge interactions from all species
```bash
Rscript get.interactions.R
```

Generate graph of known target protein interactions 
```bash
Rscript interactionGraph.R
```

### Training
Train autoencoder. `--path` specifies where to save results.
By default, `cluster.R` and `plot.clust.R` assume this will be `data/`.
If unspecified, will be saved to `date +"%Y-%m-%d"`. 
```bash
julia autoencoder.jl --path "data/"
```

### Clustering
Generate and score clusters for randomized hyperparameters
```bash
Rscript cluster.R
```

### Visualization
Generate plots and characterize optimal clusterings
```bash
Rscript plot.clust.R
```

## Motivation
### Dimension Reduction

![alt text](https://github.com/ChristiaenLab/CrobustaScreen/blob/main/fig/encode.dot.svg?raw=true)

Sample parameters are often strongly correlated. This is undesirable for self-supervised learning because each parameter additively contributes to distance used for clustering, resulting in disproportionate weight being given to phenotypes captured by multiple parameters. Linear methods of dimenison reduction (e.g. PCA) assume that all variables are independent and can be linearly combined. We could not assume that all of our measured input parameters were independent, so we instead used an autoencoder for dimension reduction.

An autoencoder is a neural network architecture widely used for denoising and image recognition. It works by encoding the input data into a lower dimensional representation that can be decoded with minimal loss. By extracting this lower dimensional encoding (the "bottleneck" or "embedding" layer), an autoencoder [can be used for dimension reduction](https://doi.org/10.1016/j.neucom.2015.08.104).
This results in an embedding that corresponds to the information content of the input data rather than absolute distance in phenotype space.

`encode.R` trains four autoencoders using embedding layers of 2, 3, 7, and 14 dimensions. `cluster.R` selects the optimal embedding based on [Akaike Information Criterion](https://en.wikipedia.org/wiki/Akaike_information_criterion) defined as 

$$AIC = 2k - 2ln(\hat{L})$$ 

where $k$ is the number of parameters and $\hat{L}$ is a likelihood function, which we define as $1 - MSE$.

### Clustering Algorithm

![alt text](https://github.com/ChristiaenLab/CrobustaScreen/blob/main/fig/cluster.dot.svg?raw=true)

Euclidean distance between embeddings is used to compute a k-nearest neighbors graph. The graph is then partitioned into clusters by [modularity](https://en.wikipedia.org/wiki/Modularity_(networks)), which is defined as 

 $$\mathcal{H} = \frac{1}{2m}\sum_{c}( \,e_c - \gamma\frac{K_c^2}{2m}) $$

 where $m$ is the average degree of the graph, $e_c$ is the number of edges in cluster $c$, and $K_c$ is the number of nodes in cluster $c$. This equation has a nicely intuitive interpretation. Modularity $\mathcal{H}$ of a graph is given by the sum of how well-connected its clusters are, defined as the difference between the number of edges in the cluster and the expected number of edges given the number of nodes in the graph and average degree of a node.

Because optimizing modularity is NP-hard, we used the [leiden algorithm](https://arxiv.org/abs/1810.08473) to approximate an optimal solution.

Though modularity ensures that clusters are well-connected, the number of clusters returned is dependent on $\gamma$, which cannot be inferred from the data. A value of $k$ must also be selected for the input graph.

### Hyperparameter Selection

![alt text](https://github.com/ChristiaenLab/CrobustaScreen/blob/main/fig/sel.dot.svg?raw=true)

`cluster.R` performs clustering for 100 random $\gamma$ values between 0.01 and 3.0 for $k$ values ranging from 3 to 53.

We selected $k$ and $\gamma$ based on four validation metrics. Mean silhouette width was calculated from the euclidean distance etween embryos.

#### $k$ Selection

**GSEA**
The known protein interactions can be treated as a gene set for [GSEA](https://en.wikipedia.org/wiki/Gene_set_enrichment_analysis). Interactions can be ranked by edge count between embryos in two conditions. An enrichment score is calculated based on occurrence of known interactions near the top of the ranked list. An optimal $k$ can be selected by maximizing enrichment score.

**Gene Network**
A gene network can be created from the $k$-NN graph by drawing an edge between a pair of conditions if the $k$-NN graph is enriched in edges between embryos in that pair of conditions.
For each condition pair $(x,y)$, we use a hypergeometric test for enrichment of edges from embryos in $x$ to embryos in $y$.
We assume the null probability $p_{xy}$ to be given by

$$p_{xy}k = \frac{\binom{K_y}{k}\binom{M-K_y}{K_x-k}}{\binom{M}{K_x}}$$

where $K_y$ is the total degree of all nodes in $y$, $k$ is the number of edges from nodes in $x$ to nodes in $y$, $M$ is the total degree of all nodes in the graph, and $K_x$ is the total degree of all nodes in $x$. 
Effectively this means we consider all edges to be a population that the edges from nodes in $x$ are drawn from, and look for overrepresenation of edges connected to nodes in $y$. We consider a false disctovery rate of 0.05 to be significantly enriched. 
We define the odds ratio $OR_{xy}$ as

$$OR_{xy} = \frac{\frac{k}{K_x-k}}{\frac{K_y}{M-K_y}}$$

#### $\gamma$ Selection
After selecting a $k$-NN graph, clustering is performed for randomized $\gamma$ values. Four metrics are calculated for each clustering: $log2(error)$, enrichment score, recall, and mean silhouette width. $\gamma$ is selected by optimizing for the product of these values.

**Reduced $k$-NN classifier**
A reduced $k$-NN classifier is created from a subset of the embeddings using the clusters as labels. The remaining embeddings are used as a test set. This process is repeated 1000 times per clustering to obtain a mean error.

**GSEA**
Condition pairs for each clustering are ranked by the proportion of edges that are between embryos in the same cluster vs. between embryos in different clusters. An enrichment score can be calculates as with $k$ selection.

**Comparison to Known Protein Interactions**
A second gene network is constructed using partial modularity between pairs of conditions. We define the partial modularity $H_{xy}$ of a condition pair $(x,y)$ as 

$$ H_{xy} = \,e_{xy} - \gamma\frac{K_x\,K_y}{2M} $$

where $e_{xy}$ is the total number of edges from embryos of condition $x$ to embryos of condition $y$, $K_x$ is the total degree of all embryos of condition $x$, $K_y$ is the total degree of all embryos in condition $y$, and M is the total degree of all nodes in the graph. If $H_{xy}$ is positive, we draw an edge between genes $x$ and $y$. We then calculate a recall score by comparing this graph to the graph of known protein interactions.

**Mean Silhouette Width**
Pointwise [silhouette width](https://doi.org/10.1016/0377-0427(87)90125-7) $s(i)$ is given by 

$$s(i) = \frac{b(i) - a(i)}{max[a(i),b(i)]}$$

where $a(i)$ is the average distance between node $i$ and other nodes in the same cluster, and $b(i)$ is the average distance between $i$ and other nodes in the closest other cluster.

### Cluster Characterization

`plot.clusts.R` tests for enrichment of experimental perturbations and experimenter-labeled phenotypes in each cluster using a hypergeometric test. For each condition $c$ in each cluster $x$, we assume the probability $p_{xc}$ of the intersect between $c$ and $x$ is given by

$$p_{xc}k = \frac{\binom{n_c}{k}\binom{N-n_c}{n_x-k}}{\binom{N}{n_x}}$$

where $k$ is the number of embryos in both $x$ and $c$, $n_c$ is the number of embryos in $c$, $N$ is the total number of embryos, and $n_x$ is the total number of embryos in $x$.

We define the odds ratio $OR_{xc}$ as

$$OR_{xc} = \frac{\frac{k}{n_x-k}}{\frac{n_c}{N-n_c}}$$
