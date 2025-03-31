using DeePWAK

include("gsea.jl")
include("readDEWAKloss.jl")

loss_dewak = (dewak, G, E, Ê, X̂, X)->begin
    G = zerodiag(G)
    L = Flux.mse(X̂, X)
    ES, NES, precision, recall = loss_stringdb(G)
    1 / NES, L, ES, NES, precision, recall
end

losslabs = [:invNES, :MSE, :ES, :NES, :precision, :recall]

path = "data/DEWAK/NES/" 
    
################################
# load data
X = readmat("data/X.csv")'
E = readmat("data/E.csv")'
F = readmat("data/SAE/E.csv")'
################################


################################
# initialize models
dewak_pca = DEWAK(X; d_0 = d_pcaNES, k_0 = k_pcaNES,
              lossfn = loss_dewak, losslabs = losslabs)

dewak_enc = DEWAK(E; d_0 = d_encNES, k_0 = k_encNES,
                lossfn = loss_dewak, losslabs = losslabs)

dewak_sae = DEWAK(F; d_0 = d_saeNES, k_0 = k_saeNES,
              lossfn = loss_dewak, losslabs = losslabs)
################################


################################
# optimization
L_pca = update!(dewak_pca, G_i->loss(dewak_pca, G_i),
            d_pcaNES - 1, k_pcaNES - 1)
writecsv(hcat(L_pca...), path * "PCA", "loss_dk.csv")

L_enc = update!(dewak_enc, G_i->loss(dewak_enc, G_i),
            d_encNES - 1, k_encNES - 1)
writecsv(hcat(L_enc...), path * "autoencoder", "loss_dk.csv")

L_sae = update!(dewak_sae, G_i->loss(dewak_sae, G_i),
            d_saeNES - 1, k_saeNES - 1)
writecsv(hcat(L_sae...), path * "SAE", "loss_dk.csv")
################################
