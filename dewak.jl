using DeePWAK

# FFI for `fgsea` R package
include("gsea.jl")


################################
# hyperparameters
window_d = 5
window_k = 5
window_γ = 1.0
n_γ = 5

steps = 100
path = "data/DEWAK/MSE/" 
################################


################################
# loss function
loss_dewak = (dewak, G, E, Ê, X̂, X)->begin
    G = zerodiag(G)
    L = Flux.mse(X̂, X)
    ES, NES, precision, recall = loss_stringdb(G)
    L, ES, NES, precision, recall
end

losslabs = [:MSE, :ES, :NES, :precision, :recall]
################################

    
################################
# load data
X = readmat("data/X.csv")'
m, n = size(X)

E = readmat("data/E.csv")'
ŋ, _ = size(E)

F = readmat("data/SAE/E.csv")'
################################


################################
# initialize model
dewak_pca = DEWAK(X; d_0 = window_d + 1, k_0 = window_k + 1,
              lossfn = loss_dewak, losslabs = losslabs)
################################


################################
# optimization
## initial pass iterates `steps` times with a small window to 
## efficiently search parameter space
L_pca = @showprogress mapreduce(vcat, 1:steps) do _
    update!(dewak_pca, G_i->loss(dewak_pca, G_i),
            window_d, window_k)
end
writecsv(losslog(dewak_pca), path * "PCA/", "loss.csv")
writecsv(dewak_pca.pcs', path, "PCs.csv")

## second pass is a single step testing all `d` values given optimal `k`
## and all `k` values given optimal `d`
L_pca_dk = update!(dewak_pca, G_i->loss(dewak_pca, G_i),
            dewak_pca.d - 1, dewak_pca.k - 1)
writecsv(hcat(L_pca_dk...), path * "PCA", "loss_dk.csv")
################################


################################
# use autoencoder embeddings rather than raw data
dewak_enc = DEWAK(E; d_0 = ŋ, k_0 = dewak_pca.k,
                  lossfn = loss_dewak, losslabs = losslabs)

## we don't want to reduce dimensions further, so d_window is set to 0
L_enc = update!(dewak_enc, G_i->loss(dewak_enc, G_i),
            0, dewak_pca.k - 1)
writecsv(hcat(L_enc...), path * "autoencoder", "loss_dk.csv")
################################


################################
# use embeddings from sparse autoencoder
dewak_sae = DEWAK(F; d_0 = window_d + 1, k_0 = window_k + 1,
              lossfn = loss_dewak, losslabs = losslabs)

L_sae = @showprogress mapreduce(vcat, 1:steps) do _
    update!(dewak_sae, G_i->loss(dewak_sae, G_i),
            window_d, window_k)
end
writecsv(losslog(dewak_sae), path * "SAE", "loss.csv")
writecsv(dewak_sae.pcs', path * "SAE", "PCs.csv")

L_sae_dk = update!(dewak_sae, G_i->loss(dewak_sae, G_i),
            dewak_sae.d - 1, dewak_sae.k - 1)
writecsv(hcat(L_sae_dk...), path * "SAE", "loss_dk.csv")
################################

