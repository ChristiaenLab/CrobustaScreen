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

p_dMSE = plotloss([L_d[1].MSE, L_d[2].MSE],
                  ["PCA (k = 78)", "SAE (k = 60)"],
                  "MSE", "d";
                  markershape = markershape, markersize = markersize,
                  legend = false);
p_dES = plotloss([L_d[1].ES, L_d[2].ES],
                 ["PCA (k = 78)", "SAE (k = 60)"],
                 "ES", "d";
                 markershape = markershape, markersize = markersize,
                 legend = false);
p_dNES = plotloss([L_d[1].NES, L_d[2].NES],
                  ["PCA (k = 78)", "SAE (k = 60)"],
                  "NES", "d";
                  markershape = markershape, markersize = markersize,
                  legend = false);
p_dPrec = plotloss([L_d[1].precision, L_d[2].precision],
                   ["PCA (k = 78)", "SAE (k = 60)"],
                   "precision", "d";
                   markershape = markershape, markersize = markersize,
                   legend = false);
p_dRec = plotloss([L_d[1].recall, L_d[2].recall],
                  ["PCA (k = 78)", "SAE (k = 60)"],
                  "recall", "d";
                  markershape = markershape, markersize = markersize);

p_kMSE = plotloss([L_k[1].MSE, L_k[2].MSE, L_k[3].MSE],
                  ["PCA (d = 29)", "SAE (d = 41)", "autoencoder (d = 14)"],
                  "MSE", "k";
                  markershape = markershape, markersize = markersize, 
                  legend = false);
p_kES = plotloss([L_k[1].ES, L_k[2].ES, L_k[3].ES],
                 ["PCA (d = 29)", "SAE (d = 41)", "autoencoder (d = 14)"],
                 "ES", "k";
                 markershape = markershape, markersize = markersize, 
                 legend = false);
p_kNES = plotloss([L_k[1].NES, L_k[2].NES, L_k[3].NES],
                  ["PCA (d = 29)", "SAE (d = 41)", "autoencoder (d = 14)"],
                  "NES", "k";
                  markershape = markershape, markersize = markersize, 
                  legend = false);
p_kPrec = plotloss([L_k[1].precision, L_k[2].precision, L_k[3].precision],
                   ["PCA (d = 29)", "SAE (d = 41)", "autoencoder (d = 14)"],
                   "precision", "k";
                   markershape = markershape, markersize = markersize, 
                   legend = false);
p_kRec = plotloss([L_k[1].recall, L_k[2].recall, L_k[3].recall],
                  ["PCA (d = 29)", "SAE (d = 41)", "autoencoder (d = 14)"],
                  "recall", "k";
                  markershape = markershape, markersize = markersize);

p = plot(p_dMSE, p_dES, p_dNES, p_dPrec, p_dRec,
         p_kMSE, p_kES, p_kNES, p_kPrec, p_kRec;
         layout = (2, 5), size = (1800, 600), margin = 10mm,
         legendfontsize = 10, tickfontsize = 10, axisfontsize = 16);
savepath(p, "fig", "loss_DEWAK.pdf")
