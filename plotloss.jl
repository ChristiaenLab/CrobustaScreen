using TrainingIO

include("plotfns.jl")

L = readcsv("data/loss.csv")
plotloss([eachcol(L)...], ["training", "test"], "MSE", "fig", "loss.pdf")

L1_L2 = readcsv("data/SAE/loss.csv")
plotloss([eachcol(L1_L2[:, [1, 4]])...], ["training", "test"],
         "L1 + MSE", "fig", "L1_L2.pdf")
plotloss([eachcol(L1_L2[:, [2, 5]])...], ["training", "test"],
         "L1", "fig", "L1.pdf")
plotloss([eachcol(L1_L2[:, [3, 6]])...], ["training", "test"],
         "MSE", "fig", "L2.pdf")

L_dewak = readcsv("data/fgsea/loss_dk.csv")
L_enc = readcsv("data/fgsea/loss_enc.csv")
L_sae = readcsv("data/SAE/loss_dk.csv")

d_max = min(max(L_dewak.d...), max(L_sae.d...))
k_max = min(max(L_dewak.k...), max(L_enc.k...), max(L_sae.k...))

L_d = (L_dewak[1:d_max, :], L_sae[1:d_max, :])
L_k = (L_dewak[(d_max + 1):end, :],
       L_sae[(max(L_sae.d...) + 1):end, :],
       L_enc[2:k_max, :])

plotloss([L_d[1].MSE, L_d[2].MSE], ["PCA (k = 78)", "SAE (k = 60)"],
         "MSE", "fig", "d_MSE.pdf", "d")
plotloss([L_d[1].ES, L_d[2].ES], ["PCA (k = 78)", "SAE (k = 60)"],
         "ES", "fig", "d_ES.pdf", "d")
plotloss([L_d[1].NES, L_d[2].NES], ["PCA (k = 78)", "SAE (k = 60)"],
         "NES", "fig", "d_NES.pdf", "d")
plotloss([L_d[1].precision, L_d[2].precision], ["PCA (k = 78)", "SAE (k = 60)"],
         "precision", "fig", "d_precision.pdf", "d")
plotloss([L_d[1].recall, L_d[2].recall], ["PCA (k = 78)", "SAE (k = 60)"],
         "recall", "fig", "d_recall.pdf", "d")

plotloss([L_k[1].MSE, L_k[2].MSE, L_k[3].MSE], ["PCA (d = 29)", "SAE (d = 41)", "autoencoder (d = 14)"],
         "MSE", "fig", "k_MSE.pdf", "k")
plotloss([L_k[1].ES, L_k[2].ES, L_k[3].ES], ["PCA (d = 29)", "SAE (d = 41)", "autoencoder (d = 14)"],
         "ES", "fig", "k_ES.pdf", "k")
plotloss([L_k[1].NES, L_k[2].NES, L_k[3].NES], ["PCA (d = 29)", "SAE (d = 41)", "autoencoder (d = 14)"],
         "NES", "fig", "k_NES.pdf", "k")
plotloss([L_k[1].precision, L_k[2].precision, L_k[3].precision], ["PCA (d = 29)", "SAE (d = 41)", "autoencoder (d = 14)"],
         "precision", "fig", "k_precision.pdf", "k")
plotloss([L_k[1].recall, L_k[2].recall, L_k[3].recall], ["PCA (d = 29)", "SAE (d = 41)", "autoencoder (d = 14)"],
         "recall", "fig", "k_recall.pdf", "k")
