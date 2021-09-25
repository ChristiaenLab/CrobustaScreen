#import scipy
#import statsmodels.api as sm
#from matplotlib import cm
#from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import os
from datetime import date
import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import operator as op
import functools as ft
import matplotlib.pyplot as plt
from dewakss import denoise as dewakss
from scipy.sparse import csr_matrix
from scipy.stats import hypergeom
from scipy.spatial import distance_matrix

parser = argparse.ArgumentParser(description='Cluster embryo data using DEWAKSS to find the optimal number of neighbors and PCs.')

parser.add_argument('--params', metavar='params', type=str, default='out/params.csv',
                    help='The input matrix with embryos as rows and parameters as columns.')
parser.add_argument('--groups', metavar='groups', type=str, default='out/groups.csv',
                    help='The annotation groups to be shown in the output figures.')
parser.add_argument('--out', metavar='groups', type=str, default='raw',
                    help='The output directory.')
parser.add_argument('--pcs', metavar='pcs', type=int, default=0,
                    help='The output directory.')
parser.add_argument('--neighbors', metavar='neighbors', type=int, default=20,
                    help='The output directory.')
#parser.add_argument('--res', metavar='res', type=float, default=1,
#                    help='The resolution parameter.')
parser.add_argument("-f", "--fff", help="a dummy argument to fool ipython", default="1")

args = parser.parse_args()
print(args.params)
print(args.groups)
print(args.out)
today = date.today()
fdate = date.today().strftime('%Y_%m_%d')
outdir = os.path.join(fdate, args.out)

sc.settings.verbosity = 3
sc.settings.autoshow = False
sc.settings.autosave = True
sc.settings.figdir = outdir

os.makedirs(sc.settings.figdir, exist_ok=True)

# Load data:
df = pd.read_csv(args.params, index_col='Row.names')

ann = pd.read_csv(args.groups, index_col='Row.names')

def savePDF(fig, fname):
    fig.savefig(os.path.join(sc.settings.figdir, fname + '.pdf'), format='pdf')

def addclust(dat,res):
    out = sc.tl.leiden(dat, resolution = res, neighbors_key = 'neighbors', copy = True)
    out.obs = out.obs.rename(columns = {'leiden' : 'leiden' + str(res)})

    out = sc.tl.leiden(out, resolution = res, neighbors_key = 'denoised', copy = True)
    out.obs = out.obs.rename(columns = {'leiden' : 'leiden_denoised' + str(res)})
    return out

def addclusts(dat,res=np.arange(0.25, 2.1, 0.25)):
    out = ft.reduce(addclust,res,dat)
    return out

def runDewakss(df, ann, n_neighbors=args.neighbors, n_pcs=0, seq=True):
    adata = sc.AnnData(df.values, obs=ann)
    adata.var_names = df.columns

    if n_pcs==0:
        n_pcs=min(adata.n_obs-1, adata.n_vars-1)

    adata = sc.pp.pca(adata, n_comps=n_pcs, copy=True)
    adata = sc.pp.neighbors(
            adata, 
            n_neighbors=n_neighbors,
            n_pcs=n_pcs, 
            copy=True)

    if seq:
        n_pcs = list(range(3, n_pcs))
        n_neighbors = list(range(3,n_neighbors))

    dewaxer = dewakss.DEWAKSS(
            adata, 
            n_neighbors=n_neighbors, 
            n_pcs=n_pcs, 
            use_global_err=False)
    dewaxer.fit(adata)

    if seq:
        fig, ax = dewaxer.plot_global_performance()
        savePDF(fig, 'global')
        fig, ax = dewaxer.plot_local_performance()
        savePDF(fig, 'local')
        fig, ax = dewaxer.plot_local_neighbour_hist(hspace=1)
        savePDF(fig, 'local_neighbors')
        return dewaxer

    else:
        ddata = dewaxer.transform(adata, copy=True)
        ddata = addclusts(ddata)
        return ddata

def clusthyper(cvclust,bg,clust):
    M = len(bg)
    n = len(clust)
    N = len(cvclust)
    x = len(np.intersect1d(cvclust,clust))
    p = hypergeom.cdf(x,M,n,N)
    fc = (x/n)/(N/M)

    dat = {"M": pd.Series([M]), 
           "n": pd.Series([n]), 
           "N": pd.Series([N]), 
           "x": pd.Series([x]), 
           "p": pd.Series([p]),
           "fc": pd.Series([fc])}
    return pd.DataFrame.from_dict(dat, orient='columns')

def clustweight(cvclust,bg,clusts):
    dat = list(map(lambda x: clusthyper(cvclust,bg,x),clusts))
    df = pd.concat(dat,ignore_index=True)
    df['weight'] = df['fc']/sum(df['fc'])
    return df

def assignclust(knn,train,clusts,weights):
    res = map(
            lambda x: list(map(
                lambda y: sum(weights.loc[x][y]),
                weights.columns.values)),knn)

    print(list(res))
    return res 

def clustcv(knn,test,train,bg,clusts='denoised_leiden'):
    print(clusts)
    trainclusts = train.obs.groupby(by=clusts).groups
    bgclusts = bg.obs.groupby(by=clusts).groups
    print(len(trainclusts))
    print(len(bgclusts))

    # enrichment of bg clusts in training clusts
    weights = [clustweight(
        i,
        bg.obs.index.values,
        bgclusts.values()) for i in trainclusts.values()]
#    weights = list(map(
#        lambda x: clustweight(
#            x,
#            bg.obs.index.values,
#            bgclusts.values()
#            ),
#        trainclusts.values()))

    # assign a weight vector for all bg clusts to each training clust
    train_bg = np.asmatrix([i.weight for i in weights])
    print(train_bg.shape)
    #train_bg = list(map(lambda x: x.weight,weights))

    # get training cluster for each neighbor of each test row
    knnclust = [train.obs.loc[i][clusts].astype('int').values for i in knn]
#    knnclust = list(map(lambda x: train.obs.loc[x][clusts].astype('int').values, knn))

    def f(x):
        w = [np.asarray(train_bg[i,:]).flatten() for i in x]
        w = np.asmatrix(w)
        res = w.sum(axis=0)
        clust = np.argmax(res)
        return clust

    testw = [f(i) for i in knnclust]
    print(testw)

    pred = bg.obs.loc[test]#[clusts]
    pred = pred.assign(predicted=testw)
    #print(pred)

    actual = list(pred.groupby(by=clusts))
    testp = list(pred.groupby(by='predicted'))
    #return {'pred':pred,'actual':actual,'test':testp}

    #print(actual)
    #print(testp)

    def pr(i, actual, test):
        actual = actual[i][1].index.values
        test = test[i][1].index.values
        print(actual[i])
        print(test[i])
        print(np.intersect1d(actual[i],test[i]))
        tp = len(np.intersect1d(actual[i],test[i]))
        fp = len(np.setdiff1d(test[i],actual[i]))
        fn = len(np.setdiff1d(actual[i],test[i]))
        tn = len(np.setdiff1d(np.setdiff1d(pred.index.values,actual[i]),test[i]))
#        precision = tp/(tp+fp)
#        recall = tp/(tp+fn)
#        f1 = 2*(precision*recall)/(precision+recall)
        return {'clust' : i, 'actual' : actual, 'test' : test, "TP" : tp, "FP" : fp, "FN" : fn, "TN" : tn}#, "F1" : f1}

    def pr2(actual, test):
        actual = actual[1].index.values
        test = test[1].index.values
        tp = len(np.intersect1d(actual,test))
        fp = len(np.setdiff1d(test,actual))
        fn = len(np.setdiff1d(actual,test))
        tn = len(np.setdiff1d(np.setdiff1d(pred.index.values,actual),test))
#        precision = tp/(tp+fp)
#        recall = tp/(tp+fn)
#        f1 = 2*(precision*recall)/(precision+recall)
        return {'actual' : actual, 'test' : test, "TP" : tp, "FP" : fp, "FN" : fn, "TN" : tn}#, "F1" : f1}

    #prc = [pr(i,actual,testp) for i in actual.index]
    prc = list(map(pr2,actual,testp))
    fields = ["TP","FP","FN"]
    res = [sum([i[j] for i in prc]) for j in fields]
    #res = {j : i for i in res for j in fields}
    precision = res[0]/(res[0]+res[1])
    recall = res[0]/(res[0]+res[2])
    f1 = 2*(precision*recall)/(precision+recall)
    f1adj = f1*len(actual)
    res = {"clusters" : clusts, "nclust" : len(actual),"TP" : res[0], "FP" : res[1], "FN" : res[2], 'precision' : precision, 'recall' : recall, 'F1' : f1, 'F1adj' : f1adj}
    return res


#    # assign weight vector to each test row 
#    trainw = list(map(
#        lambda x: train_bg[int(x)], 
#        map(int,train.obs[clusts])))
#
#    # get weight sum for neighbors of each test row for each bg clust
#    testw = map(
#            lambda x: assignclust(neighbors,train,trainclusts.values(),x),
#            weights)
#
#
#    cv = testpcs.apply(
#        lambda x: assignclust(
#            x,
#            trainpcs,
#            trainclusts.values(),
#            weights,
#            neighbors))
#
#    return cv

def dewakssCV(df,ann,bg,dwparams,frac=0.1):
    print(dwparams)
    test = df.sample(frac = frac)
    anntest = ann.loc[test.index.values]

    sel = np.setdiff1d(df.index.values, test.index.values)
    train = df.loc[sel]
    anntrain = ann.loc[sel]

    trainsel = np.isin(bg.obs.index, train.index.values)
    testsel = np.isin(bg.obs.index, test.index.values)

    dists = distance_matrix(
            bg.obsm['X_pca'][testsel,:int(dwparams.PCs)],
            bg.obsm['X_pca'][trainsel,:int(dwparams.PCs)])

    knn = map(np.argsort,dists)
    knn = list(map(lambda x: train.index.values[x[:int(dwparams.neighbors)]],knn))

    testdat = sc.AnnData(test.values, obs=anntest)

    traindat = runDewakss(
            train, anntrain, 
            int(dwparams.neighbors), 
            int(dwparams.PCs), 
            False)

    clusts = list(filter(lambda x: 'leiden' in x, bg.obs.columns))

    res = list(map(
        lambda x: clustcv(knn, test.index.values, traindat, bg, x),
        clusts))

    return res

dw = runDewakss(df, ann)

# Rerun with optimal params
dwparams = dw.global_performance.iloc[np.argmin(dw.global_performance.MSE)]

ann.iloc[:,2:5] = ann.iloc[:,2:5].astype('str')

clusts = runDewakss(
        df, ann, 
        int(dwparams.neighbors), 
        int(dwparams.PCs), 
        False)

sc.tl.umap(clusts, min_dist = 1.5, spread = 12, neighbors_key = 'denoised')
sc.pl.embedding(clusts, basis='umap', color=clusts.obs, edges=True)

dists = clusts.obsp['denoised_distances']
distout = os.path.join(sc.settings.figdir, 'dists.csv')
df = pd.DataFrame(csr_matrix.todense(dists))
df.to_csv(distout, index=False, header=None)

sc.settings.figdir = os.path.join(sc.settings.figdir, 'feat')
sc.pl.embedding(clusts, basis='umap', color=clusts.var_names, edges=True)

clusts.write_csvs(os.path.join(sc.settings.figdir, 'dat'), skip_data=False)

Rcmd = "Rscript dewakssPlots.R --params " + args.params + " --clusts " + outdir
print(Rcmd)
os.system(Rcmd)

#cv = dewakssCV(df,ann,clusts,dwparams)
#cv = pd.DataFrame.from_dict(cv, orient='columns')
#
#cv01 = dewakssCV(df,ann,clusts,dwparams,0.05)
#cv01 = pd.DataFrame.from_dict(cv, orient='columns')
#
#fig, ax = plt.barh(cv['clusters'],cv['F1adj'])
#savePDF(fig,'cv')
#
#plt.barh(cv['clusters'],cv['F1adj'])
#plt.xlabel('F1*nclust')
#plt.savefig(os.path.join(sc.settings.figdir, 'cv' + '.pdf'), format='pdf')
#
#f = plt.figure()
#savePDF(f,'cv')

