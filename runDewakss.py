#import matplotlib.pyplot as plt
#import numpy as np
#import scipy
#import statsmodels.api as sm
#from matplotlib import cm
#from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import os
import scanpy as sc
import pandas as pd
from dewakss import denoise as dewakss
from datetime import date
import argparse

parser = argparse.ArgumentParser(description='Cluster embryo data using DEWAKSS to find the optimal number of neighbors and PCs.')

parser.add_argument('--params', metavar='params', type=str, default='out/params.csv',
                    help='The input matrix with embryos as rows and parameters as columns.')
parser.add_argument('--groups', metavar='groups', type=str, default='out/groups.csv',
                    help='The annotation groups to be shown in the output figures.')
parser.add_argument('--out', metavar='groups', type=str, default='out/raw',
                    help='The output directory.')

args = parser.parse_args()
print(args.params)
print(args.groups)
print(args.out)

today = date.today()
fdate = date.today().strftime('%Y_%m_%d')
sc.settings.verbosity = 3
sc.settings.autoshow = False
sc.settings.autosave = True
sc.settings.figdir = os.path.join(args.out, fdate)
#sc.settings.format='eps'
os.makedirs(sc.settings.figdir, exist_ok=True)

def savePDF(fig, fname):
    fig.savefig(os.path.join(sc.settings.figdir, fname + '.pdf'), format='pdf')

# Load data:
df = pd.read_csv(args.params, index_col='Row.names')

ann = pd.read_csv(args.groups, index_col='Row.names')

adata = sc.AnnData(df.values, obs=ann, var=df.columns)

adata.var_names = df.columns

n_comps=min(adata.n_obs-1, adata.n_vars-1)
sc.pp.pca(adata, n_comps=n_comps)
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=n_comps)

# Run dewakss
neigbours = list(range(3, 20))
npcss = list(range(3, n_comps))

dewaxer = dewakss.DEWAKSS(adata, n_neighbors=neigbours, n_pcs=npcss, use_global_err=False)
dewaxer.fit(adata)
ddata = dewaxer.transform(adata, copy=True)

sc.tl.leiden(ddata, resolution = 1)
sc.tl.umap(ddata, min_dist = 0.01)

fig, ax = dewaxer.plot_global_performance()
savePDF(fig,'global')
fig, ax = dewaxer.plot_local_performance()
savePDF(fig,'local')

ddata.write_csvs(os.path.join(sc.settings.figdir, 'dat'), skip_data=False)

sc.pl.embedding(ddata, basis='umap', color=ddata.obs.columns, edges=True)

Rcmd = "Rscript dewakss.R --params " + args.params + " --clusts " + os.path.join(sc.settings.figdir, 'dat/obs.csv') + " --out " + os.path.join(args.out, fdate)
print(Rcmd)
os.system(Rcmd)

sc.settings.figdir = os.path.join(sc.settings.figdir, 'feat')
sc.pl.embedding(ddata, basis='umap', color=ddata.var_names,edges=True)
