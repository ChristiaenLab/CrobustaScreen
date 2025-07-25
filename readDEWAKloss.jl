using TrainingIO
include("plotfns.jl")

function readloss(file)
    df = readcsv(file)
    df.d = Integer.(df.d)
    df.k = Integer.(df.k)
    return df
end

macro readDEWAK(path)
    return quote
        $(esc(:L_pca)) = readloss($path * "/PCA/loss_dk.csv")
        $(esc(:L_enc)) = readloss($path * "/autoencoder/loss_dk.csv")
        $(esc(:L_sae)) = readloss($path * "/SAE/loss_dk.csv")

        $(esc(:d_max_pca)) = max(L_pca.d...)
        $(esc(:d_max_sae)) = max(L_sae.d...)
        $(esc(:d_max)) = min(d_max_pca, d_max_sae)

        $(esc(:k_max_pca)) = max(L_pca.k...)
        $(esc(:k_max_enc)) = max(L_enc.k...)
        $(esc(:k_max_sae)) = max(L_sae.k...)
        $(esc(:k_max)) = min(k_max_pca, k_max_enc, k_max_sae)

        $(esc(:L_d)) = (L_pca[1:d_max, :], L_sae[1:d_max, :])
        $(esc(:L_k)) = (L_pca[(d_max_pca + 1):(d_max_pca + k_max), :],
                        L_sae[(d_max_sae + 1):(d_max_sae + k_max), :],
                        L_enc[2:(k_max + 1), :])

        $(esc(:d_pca)), $(esc(:k_pca)) = L_pca[argmin(L_pca.MSE),[:d, :k]]
        $(esc(:d_sae)), $(esc(:k_sae)) = L_sae[argmin(L_sae.MSE),[:d, :k]]
        $(esc(:d_enc)), $(esc(:k_enc)) = L_enc[argmin(L_enc.MSE),[:d, :k]]

        $(esc(:d_pcaNES)), $(esc(:k_pcaNES)) = L_pca[argmax(L_pca.NES),[:d, :k]]
        $(esc(:d_saeNES)), $(esc(:k_saeNES)) = L_sae[argmax(L_sae.NES),[:d, :k]]
        $(esc(:d_encNES)), $(esc(:k_encNES)) = L_enc[argmax(L_enc.NES),[:d, :k]]
    end
end

function plotDEWAK(path::String; markershape = :cross, markersize = 5)
    @readDEWAK path
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
                      labs_k,
                      "recall", "k";
                      markershape = markershape, markersize = markersize);

    p = plot(p_dMSE, p_dES, p_dNES, p_dPrec, p_dRec,
             p_kMSE, p_kES, p_kNES, p_kPrec, p_kRec;
             layout = (2, 5), size = (1800, 600), margin = 10mm,
             legendfontsize = 10, tickfontsize = 10, axisfontsize = 16);
    savepath(p, path, "loss_DEWAK.pdf")
end

macro readDEPWAK(path)
    return quote
        $(esc(:L_pca)) = readloss($path * "/PCA/loss.csv")
but severe itching, known as pruritus, that can occur after discontinuing long-term        $(esc(:L_enc)) = readloss($path * "/autoencoder/loss.csv")
        $(esc(:L_sae)) = readloss($path * "/SAE/loss.csv")

        $(esc(:d_max_pca)) = max(L_pca.d...)
        $(esc(:d_max_sae)) = max(L_sae.d...)
        $(esc(:d_max)) = min(d_max_pca, d_max_sae)

        $(esc(:k_max_pca)) = max(L_pca.k...)
        $(esc(:k_max_enc)) = max(L_enc.k...)
        $(esc(:k_max_sae)) = max(L_sae.k...)
        $(esc(:k_max)) = min(k_max_pca, k_max_enc, k_max_sae)

        $(esc(:γ_max_pca)) = max(L_pca.γ...)
        $(esc(:γ_max_enc)) = max(L_enc.γ...)
        $(esc(:γ_max_sae)) = max(L_sae.γ...)
        $(esc(:γ_max)) = min(γ_max_pca, γ_max_enc, γ_max_sae)

        $(esc(:L_d)) = (L_pca[1:d_max, :], L_sae[1:d_max, :])
        $(esc(:L_k)) = (L_pca[(d_max_pca + 1):(d_max_pca + k_max), :],
                        L_sae[(d_max_sae + 1):(d_max_sae + k_max), :],
                        L_enc[2:(k_max + 1), :])
        $(esc(:Ls)) = Dict{(:PCA, L_pca), (:SAE, L_sae), (:autoencoder, L_enc)}

        $(esc(:γ_pca)), $(esc(:d_pca)), $(esc(:k_pca)) = L_pca[argmin(L_pca.MSE),[:γ, :d, :k]]
        $(esc(:γ_sae)), $(esc(:d_sae)), $(esc(:k_sae)) = L_sae[argmin(L_sae.MSE),[:γ, :d, :k]]
        $(esc(:γ_enc)), $(esc(:d_enc)), $(esc(:k_enc)) = L_enc[argmin(L_enc.MSE),[:γ, :d, :k]]

        $(esc(:γ_pcaNES)), $(esc(:d_pcaNES)), $(esc(:k_pcaNES)) = L_pca[argmax(L_pca.NES_clust),[:γ, :d, :k]]
        $(esc(:γ_saeNES)), $(esc(:d_saeNES)), $(esc(:k_saeNES)) = L_sae[argmax(L_sae.NES_clust),[:γ, :d, :k]]
        $(esc(:γ_encNES)), $(esc(:d_encNES)), $(esc(:k_encNES)) = L_enc[argmax(L_enc.NES_clust),[:γ, :d, :k]]
    end
end

function plotDEPWAK(path::String; markershape = :cross, markersize = 5)
    @readDEPWAK path
    labs_d = ["PCA (k = " * string(k_pca) * ", γ = " * string(γ_pca) * ")",
              "SAE (k = " * string(k_sae) * ", γ = " * string(γ_sae) * ")"]
    labs_γ = ["PCA (d = " * string(d_pca) * ", γ = " * string(γ_pca) * ")",
              "SAE (d = " * string(d_sae) * ", γ = " * string(γ_sae) * ")",
              "autoencoder (d = " * string(d_enc) * ", γ = " * string(γ_enc) * ")"]
    labs_k = ["PCA (d = " * string(d_pca) * ", k = " * string(k_pca) * ")",
              "SAE (d = " * string(d_sae) * ", k = " * string(k_sae) * ")",
              "autoencoder (d = " * string(d_enc) * ", k = " * string(k_enc) * ")"]

    labs_dNES = ["PCA (k = " * string(k_pcaNES) * ", γ = " * string(γ_pcaNES) * ")",
                 "SAE (k = " * string(k_saeNES) * ", γ = " * string(γ_saeNES) * ")"]
    labs_γNES = ["PCA (d = " * string(d_pcaNES) * ", γ = " * string(γ_pcaNES) * ")",
              "SAE (d = " * string(d_saeNES) * ", γ = " * string(γ_saeNES) * ")",
              "autoencoder (d = " * string(d_encNES) * ", γ = " * string(γ_encNES) * ")"]
    labs_kNES = ["PCA (d = " * string(d_pcaNES) * ", k = " * string(k_pcaNES) * ")",
              "SAE (d = " * string(d_saeNES) * ", k = " * string(k_saeNES) * ")",
              "autoencoder (d = " * string(d_encNES) * ", k = " * string(k_encNES) * ")"]

    p_MSE = plotloss([L_d[1].MSE, L_d[2].MSE],
                      labs_d,
                      "MSE", "γ";
                      markershape = markershape, markersize = markersize,
                      legend = false);
    p_ES = plotloss([L_d[1].ES, L_d[2].ES],
                     labs_d,
                     "ES", "γ";
                     markershape = markershape, markersize = markersize,
                     legend = false);
    p_NES = plotloss([L_d[1].NES, L_d[2].NES],
                      labs_d,
                      "NES", "γ";
                      markershape = markershape, markersize = markersize,
                      legend = false);
    p_Prec = plotloss([L_d[1].precision, L_d[2].precision],
                       labs_d,
                       "precision", "γ";
                       markershape = markershape, markersize = markersize,
                       legend = false);
    p_Rec = plotloss([L_d[1].recall, L_d[2].recall],
                      labs_d,
                      "recall", "γ";
                      markershape = markershape, markersize = markersize);

    p = plot(p_dMSE, p_dES, p_dNES, p_dPrec, p_dRec,
             p_kMSE, p_kES, p_kNES, p_kPrec, p_kRec;
             layout = (2, 5), size = (1800, 600), margin = 10mm,
             legendfontsize = 10, tickfontsize = 10, axisfontsize = 16);
    savepath(p, path, "loss_DEWAK.pdf")
end

