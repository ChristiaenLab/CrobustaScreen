using DeePWAK

include("clustfns.jl")
include("readDEWAKloss.jl");

################################
# hyperparams
window_d = 5
window_k = 5
window_γ = 1.5
n_γ = 1000
steps = 100

path = "data/DEPWAK/" 
losslabs = [:combined, :MSE, :ES, :NES, :precision, :recall, :silhouette, 
            :ES_clust, :NES_clust, :precision_clust, :recall_clust,
            :precision_network, :recall_network]

################################


################################
# loss functions
loss_dewak = (dewak, G, E, Ê, X̂, X)->begin
    G = zerodiag(G)
    L = Flux.mse(X̂, X)
    ES, NES, precision, recall = loss_stringdb(G)
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

loss_depwak = (depwak, G)->begin
    X = data(depwak)
    E = encode(depwak, X)
    D = dist(depwak)
    G = zerodiag(G)
    P = partition(depwak)
    Ĝ = G .* P
    Ê = (wak(Ĝ .* D) * E')'
    X̂ = decode(depwak, Ê)

    sil = silhouette(dist(depwak), P)
    G_network = condnetwork(G, P, depwak.γ)
    
    prec_intx, rec_intx = cv_edge(G_network, G_cond)
    ES_cl, NES_cl, prec_cl, rec_cl = loss_stringdb(G, P)
    L, ES, NES, prec, rec = loss_dewak(depwak, Ĝ, E, Ê, X̂, X)
    combined = 1 / (NES_cl * rec_cl)

    combined, L, ES, NES, prec, rec,
    sil, ES_cl, NES_cl, prec_cl,
    rec_cl, prec_intx, rec_intx
end

writeloss! = (depwak, path)->begin
    L = update!(depwak, G->loss_depwak(depwak, G),
                0, 0, window_γ, n_γ)
    writecsv(losslog(depwak), path, "loss.csv")
    return L
end

# initialize models
f_init = (d_0, k_0)->begin
    dewak = DEWAK(X; d_0 = d_0, k_0 = k_0,
                  lossfn = loss_dewak, losslabs = losslabs)
    depwak = DEPWAK(dewak, pyleiden; graphfn = pygraph)
end

writeDEPWAK = path->begin
    depwak_pca = f_init(d_pca, k_pca)
    depwak_enc = f_init(d_enc, k_enc)
    depwak_sae = f_init(d_sae, k_sae)

    L_pca = writeloss!(depwak_pca, path * "/DEPWAK/PCA")
    L_enc = writeloss!(depwak_enc, path * "/DEPWAK/autoencoder")
    L_sae = writeloss!(depwak_sae, path * "/DEPWAK/SAE")
end
################################

    
################################
# load data
X = readmat("data/X.csv")'
E = readmat("data/E.csv")'
F = readmat("data/SAE/E.csv")'
################################


################################
@readDEWAK "data/DEWAK/MSE/"
writeDEPWAK("data/DEWAK/MSE/")

@readDEWAK "data/DEWAK/NES/"
writeDEPWAK("data/DEWAK/NES/")

depwak_pca = f_init(d_pca, k_pca)
depwak_enc = f_init(d_enc, k_enc)
depwak_sae = f_init(d_sae, k_sae)

depwak_pcaNES = f_init(d_pcaNES, k_pcaNES)
depwak_encNES = f_init(d_encNES, k_encNES)
depwak_saeNES = f_init(d_saeNES, k_saeNES)
################################


################################
# optimization
L_pca = writeloss!(depwak_pca, "PCA")
L_enc = writeloss!(depwak_enc, "autoencoder")
L_sae = writeloss!(depwak_sae, "SAE")

L_pcaNES = writeloss!(depwak_pcaNES, "PCA_NES")
L_encNES = writeloss!(depwak_encNES, "autoencoder_NES")
L_saeNES = writeloss!(depwak_saeNES, "SAE_NES")
################################
