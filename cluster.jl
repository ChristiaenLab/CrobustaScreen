using DeePWAK

include("clustfns.jl")
include("readDEWAKloss.jl");


################################
# loss functions
loss_dewak = (dewak, G, E, Ê, X̂, X)->begin
    G = zerodiag(G)
    L = Flux.mse(X̂, X)
    ES, NES, precision, recall = loss_stringdb(G, P)
    L, ES, NES, precision, recall
end

loss_depwak = (depwak, G, E, Ê, X̂, X)->begin
    G = zerodiag(G)
    P = partition(depwak)
    Ĝ = G .* P
    sil = silhouette(dist(depwak), P)
    G_network = condnetwork(G, P, depwak.γ)
    
    prec_intx, rec_intx = cv_edge(G_network, G_cond)
    ES_cl, NES_cl, prec_cl, rec_cl = loss_stringdb(G, P)
    L, ES, NES, prec, rec = loss_dewak(depwak, Ĝ, E, Ê, X̂, X)

    1 / NES_cl, L, ES, NES, prec, rec,
    sil, ES_cl, NES_cl, prec_cl,
    rec_cl, prec_intx, rec_intx
end

losslabs = [:invNES_clust, :MSE, :ES, :NES, :precision, :recall, :silhouette, 
            :ES_clust, :NES_clust, :precision_clust, :recall_clust,
            :precision_network, :recall_network]
################################

    
################################
# load data
X = readmat("data/X.csv")'
E = readmat("data/E.csv")'
F = readmat("data/SAE/E.csv")'
################################


################################
# hyperparams
window_d = 5
window_k = 5
window_γ = 1.0
n_γ = 5

steps = 100
path = "data/DEPWAK/" 
################################


################################
# initialize models
dewak_pca = DEWAK(X; d_0 = d_pca, k_0 = k_pca,
                  lossfn = loss_depwak, losslabs = losslabs)
depwak_pca = DEPWAK(dewak_pca, pyleiden; graphfn = pygraph)

dewak_enc = DEWAK(X; d_0 = d_enc, k_0 = k_enc,
                  lossfn = loss_depwak, losslabs = losslabs)
depwak_enc = DEPWAK(dewak_enc, pyleiden; graphfn = pygraph)

dewak_sae = DEWAK(X; d_0 = d_sae, k_0 = k_sae,
                  lossfn = loss_depwak, losslabs = losslabs)
depwak_sae = DEPWAK(dewak_sae, pyleiden; graphfn = pygraph)
################################


################################
# optimization
L_pca = @showprogress mapreduce(vcat, 1:steps) do _
    update!(depwak_pca, G->loss(depwak_pca, G),
            0, 0, window_γ, n_γ)
end


################################
