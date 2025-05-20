using TrainingIO

include("plotfns.jl")
include("readDEWAKloss.jl")

L = readcsv("data/loss.csv")
plotloss([eachcol(L)...], ["training", "test"], "MSE";
         path = "fig", file = "loss.pdf")

L1_L2 = readcsv("data/SAE/loss.csv")
plotloss([eachcol(L1_L2[:, [1, 4]])...], ["training", "test"],
         "L1 + MSE"; path = "fig", file = "L1_L2.pdf")
p_L1 = plotloss([eachcol(L1_L2[:, [2, 5]])...], ["training", "test"],
                "L1"; legend = false);
p_L2 = plotloss([eachcol(L1_L2[:, [3, 6]])...], ["training", "test"],
         "MSE");
p_sae = plot(p_L1, p_L2; size = (600, 300), margin = 4mm);
savepath(p_sae, "fig", "loss_SAE.pdf")

markershape = :cross
markersize = 5

labs_d = ["PCA (k = " * string(k_pca) * ")",
          "SAE (k = " * string(k_sae) * ")"]
labs_k = ["PCA (d = " * string(d_pca) * ")",
          "SAE (d = " * string(d_sae) * ")",
          "autoencoder (d = " * string(d_enc) * ")"]

labs_dNES = ["PCA (k = " * string(k_pcaNES) * ")",
             "SAE (k = " * string(k_saeNES) * ")"]
labs_kNES = ["PCA (d = " * string(d_pcaNES) * ")",
             "SAE (d = " * string(d_saeNES) * ")",
             "autoencoder (d = " * string(d_encNES) * ")"]

p_dMSE = plotloss([L_d[1].MSE, L_d[2].MSE],
                  labs_d,
                  "MSE", "d";
                  markershape = markershape, markersize = markersize,
                  legend = false);
p_dES = plotloss([L_d[1].ES, L_d[2].ES],
                 labs_d,
                 "ES", "d";
                 markershape = markershape, markersize = markersize,
                 legend = false);
p_dNES = plotloss([L_d[1].NES, L_d[2].NES],
                  labs_d,
                  "NES", "d";
                  markershape = markershape, markersize = markersize,
                  legend = false);
p_dPrec = plotloss([L_d[1].precision, L_d[2].precision],
                   labs_d,
                   "precision", "d";
                   markershape = markershape, markersize = markersize,
                   legend = false);
p_dRec = plotloss([L_d[1].recall, L_d[2].recall],
                  labs_d,
                  "recall", "d";
                  markershape = markershape, markersize = markersize);

p_kMSE = plotloss([L_k[1].MSE, L_k[2].MSE, L_k[3].MSE],
                  labs_k,
                  "MSE", "k";
                  markershape = markershape, markersize = markersize, 
                  legend = false);
p_kES = plotloss([L_k[1].ES, L_k[2].ES, L_k[3].ES],
                 labs_k,
                 "ES", "k";
                 markershape = markershape, markersize = markersize, 
                 legend = false);
p_kNES = plotloss([L_k[1].NES, L_k[2].NES, L_k[3].NES],
                  labs_k,
                  "NES", "k";
                  markershape = markershape, markersize = markersize, 
                  legend = false);
p_kPrec = plotloss([L_k[1].precision, L_k[2].precision, L_k[3].precision],
                   labs_k,
                   "precision", "k";
                   markershape = markershape, markersize = markersize, 
                   legend = false);
p_kRec = plotloss([L_k[1].recall, L_k[2].recall, L_k[3].recall],
                  labs_k
                  "recall", "k";
                  markershape = markershape, markersize = markersize);

p = plot(p_dMSE, p_dES, p_dNES, p_dPrec, p_dRec,
         p_kMSE, p_kES, p_kNES, p_kPrec, p_kRec;
         layout = (2, 5), size = (1800, 600), margin = 10mm,
         legendfontsize = 10, tickfontsize = 10, axisfontsize = 16);
savepath(p, "fig", "loss_DEWAK.pdf")
