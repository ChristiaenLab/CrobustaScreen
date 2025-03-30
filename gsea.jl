using OneHotArrays, SparseArrays, Statistics
using Suppressor

delete!(ENV, "LD_LIBRARY_PATH")
using RCall

R"""
library(fgsea)
library(BiocParallel)

# Create a SerialParam object
serial_param <- SerialParam()
"""

conds = readcsv("data/conds.csv")
C = sparse(Matrix(conds)')
c,n = size(C)
G_STRINGdb = readmat("data/G_STRINGdb.csv")
n_STRINGdb = sum(G_STRINGdb)

labels = readmat("data/interactionlabels.csv")

ids = reshape(labels, c^2)
ids_sig = readcsv("data/stringdblabels.csv")[:, 1]
@rput ids ids_sig

i_stringdb = mapreduce(x->findall(ids .== x), vcat, ids_sig)

function interactionscore(G::AbstractMatrix, i::Integer, j::Integer;
                          weighted = false)
    if !weighted
        G = G .> 0
    end
    sum(G[C[i, :], C[j, :]])
end

function interactionscore(G::AbstractMatrix; scale = false)
    score = mapreduce(hcat, 1:c) do i
        mapreduce(vcat, 1:c) do j
            s = interactionscore(G, i, j)
            if scale
                s = s / (sum(C[i, :]) * sum(C[j, :]))
            end
            s
        end
    end
    reshape(score, c^2)
end

function fgsea(G::AbstractMatrix)
    S = interactionscore(G)
    @rput S
    R"""
names(S) <- ids
res <- fgsea(list(score=ids_sig),S,
             BPPARAM=serial_param,scoreType="pos")
"""
    @rget res
    return res
end

function loss_GSEA(loss::Function, G)
    gsea = @suppress fgsea(G)
    L = loss(G)
    return L, gsea.ES, gsea.NES
end

function loss_GSEA(M::AbstractDEWAK, G)
    gsea = @suppress fgsea(G)
    L = loss(M,G)
    return L, gsea.ES, gsea.NES
end

function loss_DEWAK(M, d::Integer, k::Integer)
    G = kern(M.cache, d, k)
    loss_GSEA(M, G)
end

function cv_edge(G, G_test)
    TP = sum(G .* G_test)
    FP = sum(G .* (1 .- G_test))
    TN = sum((1 .- G) .* (1 .- G_test))
    FN = sum((1 .- G) .* G_test)
    TP,FP,TN,FN
end

function pr_stringdb(G)
    TP = sum(G .* G_STRINGdb)
    precision = TP / sum(G)
    recall = TP / n_STRINGdb
    return precision, recall
end

function es(N::Integer, R::AbstractVector{<:Integer}, ties::AbstractVector)
    N_H = length(R)                    # Number of genes in the gene set
    P_hit = 1.0 / N_H                  # Increment for hits
    P_miss = 1.0 / (N - N_H)           # Decrement for misses

    # Initialize steps with P_miss decrements
    steps = fill(-P_miss, N)

    # Update steps at indices of hits with P_hit increments
    steps[R] .= P_hit
    map(i->begin
            P_tie = mean(steps[i])
            steps[i] .= P_tie
        end, ties)

    # Compute the running sum (cumulative sum)
    running_sum = cumsum(steps)

    # Find the maximum deviation from zero
    idx = findmax(abs.(running_sum))[2]
    ES = running_sum[idx]

    return ES
end

function es(S::AbstractVector{<:Integer}, i::AbstractVector)
    n = size(S, 1)
    R = sortperm(S, rev = true)
    S = S[R]
    ties = map(x->S .== x,unique(S))
    return es(n, R[i], ties)
end

function es(S::AbstractVector{<:Real}, i::AbstractVector)
    n = size(S, 1)
    R = sortperm(S, rev = true)
    S = S[R]
    ties = [S .== 0]
    return es(n, R[i], ties)
end

function es_stringdb(G)
    scores = interactionscore(G)
    es(scores, i_stringdb)
end

function loss_stringdb(G)
    precision,recall = pr_stringdb(G)
    #ES = es_stringdb(G)
    gsea = @suppress fgsea(G)
    ES = gsea.ES
    NES = gsea.NES
    ES, NES, precision, recall
end

