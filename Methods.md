<script
  src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"
  type="text/javascript">
</script>

# Methods
Code and its usage are available at [CrobustaScreen](https://github.com/ChristiaenLab/CrobustaScreen).

## Preprocessing

### Summary statistics 
The pipeline uses features extracted from segmentation of confocal images using Imaris. Summary statistics are extracted for segmened cells in each embryo. To generate embryo-level statistics, we compute the maximum, minimum, mean, and standard deviation for these statistics for TVCs and ATMs. We do the same for distances between TVCs and ATMs and the cosines of angles between TVCs. Angle cosines between ATMs were not computed because WT embryos have only two ATMs. 

**t.cond.param.pdf**
t-test for enrichment of summary statistics in each embryo treatment. Dot size indicates FDR-corrected p-values. Dot color indicates fold change relative to mean value in all samples. Only statistics with FDR > 0.05 and log2 fold change > 0.25 or < -0.25 in at least one condition are shown.

### Normalization 
From the cell segmentation statistics 114 embryo-level parameters are computed. Parameters are normalized by z-score then scaled between -1 and 1.

**z.pdf**
z-scores for all embryos.

## Autoencoder
Sample parameters are often strongly correlated. This is undesirable for self-supervised learning because each parameter additively contributes to distance used for clustering, resulting in disproportionate weight being given to phenotypes captured by multiple parameters. Linear methods of dimenison reduction (e.g. PCA) assume that all variables are independent and can be linearly combined. We could not assume that all of our measured input parameters were independent, so we instead used an autoencoder for dimension reduction.

An autoencoder is a neural network architecture widely used for denoising and image recognition. It works by encoding the input data into a lower dimensional representation that can be decoded with minimal loss. By extracting this lower dimensional encoding (the "bottleneck" or "embedding" layer), an autoencoder [can be used for dimension reduction](https://doi.org/10.1016/j.neucom.2015.08.104).
This results in an embedding that corresponds to the information content of the input data rather than absolute distance in phenotype space.

### Architecture
The normalized scaled parameters were passed to a multilayer perceptron with five hidden layers and a `tanh` activation function. Hidden layer sizes were $58 \to 29 \to 14 \to 29 \to 58$.

Training was done using [Flux.jl](https://github.com/FluxML) on an NVIDIA GeForce RTX 3060 GPU.

### Hyperparameters
The autoencoder was trained for 10000 epochs with learning rate 0.0001 and weight decay 0.0001. 10% of the data were held back as a test set.

### Dimension reduction
To use the autoencoder for dimension reduction, data are passed to only the first three hidden layers, yielding a 14-dimensional embedding.

## Putative regulatory network
We attempt to validate our phenotypically-derived gene regulatory network by comparing it to a putative network of known interactions.

### Protein-protein interactions
The putative network combines known interactions from orthologous genes in *Ciona intestinalis*, mouse, and human.

### Ortholog lookup
We use [STRINGdb](https://doi.org/10.1093/nar/gkq973) to construct a known protein interaction network of the perturbed genes. Because the *C. robusta* network is poorly characterized, we use ENSEMBL to obtain orthologs from *M. musculus* and *H. sapiens*. 

## $k$-NN graph
Distances were calculated between all samples in embedding space, which were used to select the $k$ nearest neighbors for a gene interaction graph.

### $k$ selection
We used a semisupervised method to select an optimal $k$ from all values between 3 and 53.

**GSEA**
The known protein interactions can be treated as a gene set for [GSEA](https://en.wikipedia.org/wiki/Gene_set_enrichment_analysis). Interactions can be ranked by edge count between embryos in two conditions. An enrichment score is calculated based on occurrence of known interactions near the top of the ranked list. An optimal $k$ can be selected by maximizing enrichment score. Our implementation uses the `fgsea` R library.

**knn.pdf**
Optimal $k$-NN graph ($k = 15$). Samples are plotted as vertices after UMAP dimension reduction with default parameters.

**Gene Network**
A gene network can be created from the $k$-NN graph by drawing an edge between a pair of conditions if the $k$-NN graph is enriched in edges between embryos in that pair of conditions.
For each condition pair $(x,y)$, we use a hypergeometric test for enrichment of edges from embryos in $x$ to embryos in $y$.
We assume the null probability $p_{xy}$ to be given by

$$p_{xy}k = \frac{\binom{K_y}{k}\binom{M-K_y}{K_x-k}}{\binom{M}{K_x}}$$

where $K_y$ is the total degree of all nodes in $y$, $k$ is the number of edges from nodes in $x$ to nodes in $y$, $M$ is the total degree of all nodes in the graph, and $K_x$ is the total degree of all nodes in $x$. 
Effectively this means we consider all edges to be a population that the edges from nodes in $x$ are drawn from, and look for overrepresenation of edges connected to nodes in $y$. We consider a false disctovery rate of 0.05 to be significantly enriched. 
We define the odds ratio $OR_{xy}$ as

$$OR_{xy} = \frac{\frac{k}{K_x-k}}{\frac{K_y}{M-K_y}}$$

**conditionEdgeNetworkUp.pdf**
Gene network as determined by enrichment of edges between samples in two treatments. Color of edges indicate log2 odds ratio of the condition pair as described above. Only significantly enriched edges are shown, not significantly depleted edges.

## Clustering
Clustering on the graph is performed using the leiden algorithm via [`leidenalg`](github.com/vtraag/leidenalg) with default parameters.

### $\gamma$ selection
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

**optimization.pdf**
top two rows: evaluation of randomly sampled resolution values.
bottom: GSEA enrichment score for each value of $k$.

top left: ES * recall; used for selecting optimal resolution.
top center: GSEA enrichment score for each resolution value.
top right: fraction of edges predicted by the putative regulatory network which connect samples in the same cluster.
middle left: log2 error of the reduced $k$-NN classifier as described above.
middle center: Mean silhouette width as described above.
middle left: Number of clusters for each resolution value.

### Cluster Characterization

`plot.clusts.R` tests for enrichment of experimental perturbations and experimenter-labeled phenotypes in each cluster using a hypergeometric test. For each condition $c$ in each cluster $x$, we assume the probability $p_{xc}$ of the intersect between $c$ and $x$ is given by

$$p_{xc}k = \frac{\binom{n_c}{k}\binom{N-n_c}{n_x-k}}{\binom{N}{n_x}}$$

where $k$ is the number of embryos in both $x$ and $c$, $n_c$ is the number of embryos in $c$, $N$ is the total number of embryos, and $n_x$ is the total number of embryos in $x$.

We define the odds ratio $OR_{xc}$ as

$$OR_{xc} = \frac{\frac{k}{n_x-k}}{\frac{n_c}{N-n_c}}$$

**knn.clust.pdf**
Same as `knn.pdf` with samples labeled by cluster.

**hyper.cond.pdf**
Hypergeometric test for enrichment of conditions in each cluster. Dot color indicates log2 odds ratio compared to background frequency of each condition. Dot size indicates false discovery rate. Note that a significant negative value indicates significant depletion.

**hyper.pheno.pdf**
Hypergeometric test for enrichment of labeled phenotypes in each cluster. 

**t.clust.param.pdf**
t-test for enrichment of each raw parameter in each cluser. Dot size indicates FDR-corrected p-values. Dot color indicates fold change relative to mean value in all samples. Only statistics with FDR > 0.05 and log2 fold change > 0.25 or < -0.25 in at least one condition are shown.
