using DeePWAK
using PyCall

include("gsea.jl")

# score for enrichment of interactions within clusters
# where `P` is an adjacency matrix partitioning `G` into clusters
function interactionscore(G::AbstractMatrix, G_inv::AbstractMatrix,
                          i::Integer, j::Integer; scale = false)
    s = (sum ∘ intxsub)(G, i, j)
    s_inv = (sum ∘ intxsub)(G_inv, i, j)
    if scale
        s = s / (sum(C[i, :]) * sum(C[j, :]))
        s_inv = s_inv / (sum(C[i, :]) * sum(C[j, :]))
    end
    s / s_inv
end

function interactionscore(G::AbstractMatrix, P::AbstractMatrix;
                          scale = false, weighted = false)
    if !weighted
        G = G .> 0
    end
    Ĝ = G .* P
    Ĝ_inv = G .* (1 .- P)
    score = mapintx(interactionscore, Ĝ, Ĝ_inv; scale = scale)
    score[isnan.(score)] .= 0
    score = reshape(score, c^2)
    sub = filter(isfinite, score)
    if length(sub) > 0
        score[!isfinite(score)] .= maximum(sub)
    else
        score .= 0
    end
    score
end

function fgsea(G::AbstractMatrix, P::AbstractMatrix)
    scores = interactionscore(G, P)
    fgsea(scores)
end

function loss_stringdb(G::AbstractMatrix, P::AbstractMatrix;
                       weighted = false)
    precision, recall = pr_stringdb(G .* P; weighted = weighted)
    #ES = es_stringdb(G)
    gsea = @suppress fgsea(G, P)
    ES = gsea.ES
    NES = gsea.NES
    ES, NES, precision, recall
end

@pyimport leidenalg
@pyimport igraph

function pygraph(G)
    (igraph.Graph[:Adjacency] ∘ zerodiag)(G)
end

function pygraph(D, G)
    (igraph.Graph[:Adjacency] ∘ zerodiag)(G)
end

function pygraphw(G)
    (igraph.Graph[:Weighted_Adjacency] ∘ zerodiag)(G)
end

function pygraphw(D, G)
    (igraph.Graph[:Weighted_Adjacency] ∘ zerodiag)(D .* G)
end


function pyleiden(graph, γ, args...; kwargs...)
    clusts = leidenalg.find_partition(graph,
                             leidenalg.RBConfigurationVertexPartition,
                             args...;
                             n_iterations = -1,
                             resolution_parameter = γ, kwargs...)
    clusts.membership .+ 1
end

function jleiden(G, γ)
    leiden(G, "cpm"; γ = γ)
end

function partitionloss(M::DEPWAK)
    X = data(M)
    D = dist(M)
    P = partitionmat(M)
    G = wak(D .* P)
    Flux.mse((G * X')', X)
end

function intxmod(G::AbstractMatrix, P::AbstractMatrix, γ::AbstractFloat,
                 i::Integer, j::Integer)
    rbpm(intxsub(G, i, j), intxsub(P, i, j), γ)
end

function intxmod(G::AbstractMatrix, P::AbstractMatrix, γ::AbstractFloat)
    mapintx(intxmod, G, P, γ)
end
    
function condnetwork(G::AbstractMatrix, P::AbstractMatrix, γ::AbstractFloat)
    Ĝ = intxmod(G, P, γ)
    Ĝ[Ĝ .< 0] .= 0
    Ĝ
end

