using DeePWAK

include("gsea.jl")
include("readDEWAKloss.jl")


################################
@readDEWAK "data/DEWAK/MSE/"

losslabs = [:invNES, :MSE, :ES, :NES, :precision, :recall]

window_d = 5
window_k = 5
steps = 100
path = "data/DEWAK/NES/" 

d_0 = max(d_pcaNES, d_saeNES, d_encNES)
k_0 = max(k_pcaNES, k_saeNES, k_encNES)
################################


################################
loss_dewak = (dewak, G, E, Ê, X̂, X)->begin
    G = zerodiag(G)
    L = Flux.mse(X̂, X)
    ES, NES, precision, recall = loss_stringdb(G)
    1 / NES, L, ES, NES, precision, recall
end

writeloss! = (dewak,path)->begin
    @showprogress mapreduce(vcat, 1:steps) do _
        update!(dewak, G_i->loss(dewak, G_i),
                window_d, window_k)
    end
    writecsv(losslog(dewak), path, "loss.csv")
    writecsv(dewak.pcs', path, "PCs.csv")
end
################################


################################
# load data
X = readmat("data/X.csv")'
E = readmat("data/E.csv")'
F = readmat("data/SAE/E.csv")'
################################


################################
# initialize models
dewak_pca = DEWAK(X; d_0 = d_0, k_0 = k_0,
              lossfn = loss_dewak, losslabs = losslabs)

dewak_enc = DEWAK(E; d_0 = d_encNES, k_0 = k_0,
                lossfn = loss_dewak, losslabs = losslabs)

dewak_sae = DEWAK(F; d_0 = d_0, k_0 = k_0,
              lossfn = loss_dewak, losslabs = losslabs)
################################


################################
# optimization
writeloss!(dewak_pca, path * "PCA/")
writeloss!(dewak_sae, path * "SAE/")
writeloss!(dewak_enc, path * "autoencoder/")
################################


################################
L_pca = update!(dewak_pca, G_i->loss(dewak_pca, G_i),
            d_0 - 1, k_0 - 1)
writecsv(hcat(L_pca...), path * "PCA", "loss_dk.csv")

L_enc = update!(dewak_enc, G_i->loss(dewak_enc, G_i),
            d_encNES - 1, k_0 - 1)
writecsv(hcat(L_enc...), path * "autoencoder", "loss_dk.csv")

L_sae = update!(dewak_sae, G_i->loss(dewak_sae, G_i),
            d_0 - 1, k_0 - 1)
writecsv(hcat(L_sae...), path * "SAE", "loss_dk.csv")
################################
