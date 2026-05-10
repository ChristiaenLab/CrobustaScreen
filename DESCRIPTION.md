# Figures

## Parameters
**out/params.cond.pdf**
Raw parameter values split by condition.

**out/z.cond.pdf**
z scores split by condition.

**out/embedding.cond.pdf**
embedding values split by condition.

**out/params.clust.pdf**
Raw parameter values split by cluster.

**out/z.clust.pdf**
z scores split by cluster.

**out/embedding.clust.pdf**
embedding values split by cluster.

## Hyperparameter selection

### ES.pdf
GSEA enrichment score for known interactions by $k$ value.

### leiden.k15.optimization.pdf
Statistics for selecting optimal resolution score.

**combinned_score**
`ES * recall`. The value used to select an optimum.

**ES vs. resolution**
GSEA enrichment score for known interactions between embryos in the same cluster.

**recall**
Fraction of known interactions present within any cluster.

**log2error**
log2 classification error for a reduced $k$-NN classifier created from a subset of the embeddings using the clusters as labels. The remaining embeddings are used as a test set. This process is repeated 1000 times per clustering to obtain a mean error.

**mean_silhouette**
Pointwise [silhouette width](https://doi.org/10.1016/0377-0427(87)90125-7) $s(i)$ is given by 

$$s(i) = \frac{b(i) - a(i)}{max[a(i),b(i)]}$$

where $a(i)$ is the average distance between node $i$ and other nodes in the same cluster, and $b(i)$ is the average distance between $i$ and other nodes in the closest other cluster.

**nclust**
Number of clusters for each resolution value.

## kNN

**knn.pdf**
kNN of embryos.

**knn.cond.pdf**
As above with subset of conditions indicated.

**knn.clust.pdf**
As above but labeled by cluster.

### umap/point/
As `knn.clust.pdf`, but each plot shows black circles indicating embryos of one condition.

### umap/edge/
As above, but clusters are indicated by coloring edges.

## Gene network

**knn.network/conditionEdgePos.pdf**
Poisson test for enrichment of edges between each condition.

**knn.network/conditionEdgeNetworkUp.pdf**
Gene network omitting significant depletion of edges.

**knn.network.fr/conditionEdgeNetworkUp.pdf**
As above with an alternate node layout algorithm. The edges themselves are unchanged.

## Cluster characterization

### params/

**t.pdf**
t-test for enrichment of each raw parameter in each cluser. Dot size indicates FDR-corrected p-values. Dot color indicates fold change relative to mean value in all samples. Only statistics with FDR > 0.05 and log2 fold change > 0.25 or < -0.25 in at least one condition are shown.

**t.fc.pdf**
As above, but instead of mean value within a cluster, color indicates mean fold change from the sample mean.

**u.pdf**
As above, but with a u test rather than a t test.

### embeddings/
As figures in `params/` but testing embeddings rather than parameters.

### z/
As figures in `params/` but testing z score normalized parameters.

### condition/

**hyper.pdf**
Hypergeometric test for enrichment of conditions in each cluster.

**params/**
As above, but testing significance of parameters in each condition.

**embeddings/**
As above, but testing significance of embeddings in each condition.

**z/**
As above, but testing significance of z score normalized parameters in each condition.

### pheno/

**hyper.pdf**
Hypergeometric test for enrichment of labeled phenotypes in each cluster. 

