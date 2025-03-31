using TrainingIO

L_pca = readcsv("data/DEWAK/MSE/PCA/loss_dk.csv")
L_enc = readcsv("data/DEWAK/MSE/autoencoder/loss_dk.csv")
L_sae = readcsv("data/DEWAK/MSE/SAE/loss_dk.csv")

d_max_pca = max(L_pca.d...)
d_max_sae = max(L_sae.d...)
d_max = min(d_max_pca, d_max_sae)

k_max_pca = max(L_pca.k...)
k_max_enc = max(L_enc.k...)
k_max_sae = max(L_sae.k...)
k_max = min(k_max_pca, k_max_enc, k_max_sae)

L_d = (L_pca[1:d_max, :], L_sae[1:d_max, :])
L_k = (L_pca[(d_max_pca + 1):(d_max_pca + k_max), :],
       L_sae[(d_max_sae + 1):(d_max_sae + k_max), :],
       L_enc[2:(k_max + 1), :])

d_pca, k_pca = L_pca[argmin(L_pca.MSE),[:d, :k]]
d_sae, k_sae = L_sae[argmin(L_sae.MSE),[:d, :k]]
d_enc, k_enc = L_enc[argmin(L_enc.MSE),[:d, :k]]

d_pcaNES, k_pcaNES = L_pca[argmax(L_pca.NES),[:d, :k]]
d_saeNES, k_saeNES = L_sae[argmax(L_sae.NES),[:d, :k]]
d_encNES, k_encNES = L_enc[argmax(L_enc.NES),[:d, :k]]
