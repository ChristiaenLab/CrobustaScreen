using DeePWAK
include("clustfns.jl")

path = "data/"
method = :ES

meta = readcsv("data/groups.csv")
dat_k = readcsv(path * "k.csv")
dat = readcsv(path * "leiden.csv")

k = dat_k.k[argmax(dat_k.ES)]

clusts = Matrix(dat[:,Not(1:8)])
scores = dat[:,2:8]
sel = argmax(scores[:,method])
γ = scores.resolution[sel]
clust = clusts[sel,:]

E = readmat("data/E.csv")'
D = (zerodiag ∘ inveucl)(E)

clustix = clustdsort(clust, D)
centroids = centroidix(clust, D)
centroidIDs = meta.Name[centroids]

map(unique(clust)) do i
    ix = clustix[i]
    ids = meta.Name[ix]
    writecsv(ids, path * "sorted/", string(i) * ".csv")
end
