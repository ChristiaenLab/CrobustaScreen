using DeePWAK

using PyCall
@pyimport leidenalg
@pyimport igraph

#using RCall
#@rimport fgsea

include("gsea.jl")

function pygraph(D, G)
    igraph.Graph[:Weighted_Adjacency](D .* G)
end

function pyleiden(graph, γ, args...; kwargs...)
    clusts =leidenalg.find_partition(graph,
                             leidenalg.RBConfigurationVertexPartition,
                             args...;
                             n_iterations = -1,
                             resolution_parameter = γ, kwargs...)
    clusts.membership .+ 1
end

function jleiden(G, γ)
    leiden(G, "cpm"; γ = γ)
end

loss_dewak = (dewak, G, E, Ê, X̂, X)->begin
    G = zerodiag(G)
    L = Flux.mse(X̂, X)
    ES, NES, precision, recall = loss_stringdb(G)
    L, ES, NES, precision, recall
end

losslabs = [:MSE, :ES, :NES, :precision, :recall]
    
################################
# load data
X = readmat("data/X.csv")'
m, n = size(X)

E = readmat("data/E.csv")'
ŋ, _ = size(E)

F = readmat("data/SAE/E.csv")'
################################

window_d = 5
window_k = 5
window_γ = 1.0
n_γ = 5

steps = 100
path = "data/" 

# initialize model
dewak = DEWAK(X; d_0 = window_d + 1, k_0 = window_k + 1,
              lossfn = loss_dewak, losslabs = losslabs)

L_dewak = @showprogress mapreduce(vcat, 1:steps) do _
    update!(dewak, G_i->loss(dewak, G_i),
            window_d, window_k)
end
writecsv(losslog(dewak), path*"fgsea", "loss_DEWAK.csv")
writecsv(dewak.pcs', path, "PCs.csv")

L_dk = update!(dewak, G_i->loss(dewak, G_i),
            dewak.d - 1, dewak.k - 1)
writecsv(hcat(L_dk...), path * "fgsea", "loss_dk.csv")

dewak_E = DEWAK(E; d_0 = ŋ, k_0 = dewak.k,
                lossfn = loss_dewak, losslabs = losslabs)
L_dewak_E = update!(dewak_E, G_i->loss(dewak_E, G_i),
            0, dewak.k - 1)
writecsv(hcat(L_dewak_E...), path * "fgsea", "loss_enc.csv")

dewak_sae = DEWAK(F; d_0 = window_d + 1, k_0 = window_k + 1,
              lossfn = loss_dewak, losslabs = losslabs)

L_sae = @showprogress mapreduce(vcat, 1:steps) do _
    update!(dewak_sae, G_i->loss(dewak, G_i),
            window_d, window_k)
end
writecsv(losslog(dewak_sae), path*"SAE", "loss_DEWAK.csv")
writecsv(dewak_sae.pcs', path * "SAE", "PCs.csv")

L_sae_dk = update!(dewak_sae, G_i->loss(dewak, G_i),
            dewak_sae.d - 1, dewak_sae.k - 1)
writecsv(hcat(L_sae_dk...), path * "SAE", "loss_dk.csv")


depwak = DEPWAK(dewak, pyleiden; graphfn = pygraph)

#
L_depwak = @showprogress mapreduce(vcat, 1:steps) do _
    update!(depwak, G->loss(depwak, G),
            window_d, window_k, window_γ, n_γ)
end

