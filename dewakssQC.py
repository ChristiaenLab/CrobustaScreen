import os
from datetime import date
import argparse
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy
import statsmodels.api as sm
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from dewakss import denoise as dewakss
 
today = date.today()
fdate = date.today().strftime('%Y_%m_%d')
sc.settings.verbosity = 3
sc.settings.autoshow = False
sc.settings.autosave = True
sc.settings.figdir = os.path.join('out', fdate)
#sc.settings.format='eps'
os.makedirs(sc.settings.figdir, exist_ok=True)

def savePDF(fig, fname):
    fig.savefig(os.path.join(sc.settings.figdir, fname + '.pdf'), format='pdf')

# Load data:
df = pd.read_csv('screen/2020-06-24/dat.csv', index_col=0)
df = df[(df.nTVC==4) & (df.nATM==2) & (df['TVC.contiguous']==1) & (df['ATM.contiguous']==1)]

ann = df[['Condition','Phenotype']]

df = df.iloc[:,np.concatenate([np.arange(0,20), 
    np.arange(26,109), np.arange(111,114), np.arange(116,120)
    ])]

adata = sc.AnnData(df.values, obs=ann, var=df.columns)

adata.var_names = df.columns

#adata.obs = adata.obs.merge(df.reset_index(), left_on=0, right_on='index', how='left')

# Make sure samples are rows and features are columns. if not data = data.T
#adata = sc.AnnData(data, obs=sample_names, var=feature_annotations)

sc.pp.pca(adata, n_comps=40)
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=40)

# Run dewakss
neigbours = list(range(3, 20))
npcss = list(range(3, 40))

dewaxer = dewakss.DEWAKSS(adata, n_neighbors=neigbours, n_pcs=npcss, use_global_err=False)
dewaxer.fit(adata)
ddata = dewaxer.transform(adata, copy=True)

sc.tl.leiden(ddata, resolution = 1)
sc.tl.umap(ddata, min_dist = 0.01)

#ddata.obs.to_csv('out/dewakssdat.csv')
#print(ddata.uns['neighbors'])#.to_csv('out/neighbors.csv')
#print(ddata.uns['neighbors']['connectivities'])#.to_csv('out/neighbors.csv')
#ddata.obsm['X_umap'].to_csv('out/umap.csv')

fig, ax = dewaxer.plot_global_performance()
savePDF(fig,'global')
#fig.savefig(os.path.join(sc.settings.figdir, 'global.pdf'), format='pdf')
fig, ax = dewaxer.plot_local_performance()
savePDF(fig,'local')
#fig.savefig(os.path.join(sc.settings.figdir, 'local.pdf'), format='pdf')

#fig = plt.figure(figsize=(4, 4), constrained_layout=False)
#ax = fig.subplots(1, 1)

#sc.pl.embedding(ddata, basis='umap', color=['leiden', 'condition', 'phenotype'], edges=True)

ddata.write_csvs(os.path.join(sc.settings.figdir, 'dat'), skip_data=False)

sc.pl.embedding(ddata, basis='umap', color=ddata.obs.columns, edges=True)

sc.settings.figdir = os.path.join(sc.settings.figdir, 'feat')
sc.pl.embedding(ddata, basis='umap', color=ddata.var_names,edges=True)

def utest(x, y):
    y = df[y]
    x = y[ddata.obs.leiden == x]
    lfc = np.log2(scipy.mean(x) / scipy.mean(y))
    if lfc==0:
        p = 1
    else:
        p = scipy.stats.mannwhitneyu(x, y).pvalue
    res = (lfc, p)
    return res


feat = list(map(lambda x:
    list(map(lambda y:
        utest(x, y),
        df.columns)), 
    set(ddata.obs.leiden)))

fc = pd.DataFrame(map(lambda x: list(map(lambda y: y[0], x)), feat)).T
p = pd.DataFrame(map(lambda x: list(map(lambda y: y[1], x)), feat)).T

fc.to_csv(os.path.join(sc.settings.figdir, 'log2FoldChange.csv'))
p.to_csv(os.path.join(sc.settings.figdir, 'p.csv'))

map(lambda x:
        df[ddata.obs==x],
        set(ddata.obs['leiden']))
        
fdr = p.apply(lambda x: sm.stats.multipletests(x, method='fdr_bh')[1])
fdr.to_csv(os.path.join(sc.settings.figdir, 'FDR.csv'))
logfdr = fdr.apply(lambda x: -np.log10(x))

pvals = pd.melt(logfdr)['value']
plim = pvals.quantile([.05, .95])
pminmax = (logfdr - plim[.05])/(plim[.95] - plim[.05])

N = 256
vals = np.ones((N, 4))
vals[:, 0] = np.linspace(90/256, 1, N)
vals[:, 1] = np.linspace(39/256, 1, N)
vals[:, 2] = np.linspace(41/256, 1, N)
colscale = ListedColormap(vals)
colors = colscale(pminmax)

fig = sc.pl.stacked_violin(ddata, ddata.var_names, groupby='leiden', swap_axes=True, colors=colors)
savePDF(fig,'features')
type(fig)
fig

sc.tl.rank_genes_groups(ddata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(ddata, n_genes=25, sharey=False)

cols = map(np.log10, ddata.uns['rank_genes_groups']['pvals_adj'].flat)
print(df)

def logfn(x):
    return -np.log10(x)

log10p = pd.DataFrame(map(lambda x: map(logfn, x), ddata.uns['rank_genes_groups']['pvals_adj']))

#fig.savefig(os.path.join(sc.settings.figdir,'features.pdf'), format='pdf')

viridis = cm.get_cmap('viridis', 256)


np.log10(ddata.uns['rank_genes_groups']['pvals_adj'].flat)

newcolors = viridis(np.linspace(0, 1, 256))

fig, ax = plt.pcolor(log10p, cmap=viridis, rasterized=False, vmin=plim[.05], vmax=plim[.95])
fig.savefig(os.path.join(sc.settings.figdir, 'tmp.pdf'), format='pdf')

def plot_examples(cms):
    """
    helper function to plot two colormaps
    """
    np.random.seed(19680801)
    data = log10p

    fig, axs = plt.subplots(1, 2, figsize=(6, 3), constrained_layout=True)
    for [ax, cmap] in zip(axs, cms):
        psm = ax.pcolormesh(data, cmap=cmap, rasterized=True, vmin=plim[.05], vmax=plim[.95])
        fig.colorbar(psm, ax=ax)
    savePDF(fig,'tmp')

plot_examples([viridis, newcmp])

print(plim[.05])
print(np.linspace(plim[.05],plim[.95],12))
viridis(np.linspace(plim[.05],plim[.95],12))

#fig=plt.gcf()
#fig.savefig('clust.eps', format='eps')

#fig = plt.figure(figsize=(4,4), constrained_layout=False)
##ax = fig.subplots(1, 1)
#
#sc.settings.figdir = './out/cond'
#sc.pl.embedding(ddata, basis='umap', color='condition',edges=True)
##fig=plt.gcf()
##fig.savefig('cond.eps', format='eps')
#
#fig = plt.figure(figsize=(4,4), constrained_layout=False)
##ax = fig.subplots(1, 1)
#
#sc.settings.figdir = './out/phenotype'
#sc.pl.embedding(ddata, basis='umap', color='phenotype',edges=True)
##fig=plt.gcf()
##fig.savefig('phenotype.eps', format='eps')

