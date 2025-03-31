using OneHotArrays, SparseArrays, Statistics
using Suppressor # prevent R calls from printing to stdout

# fix for bug in RCall
delete!(ENV, "LD_LIBRARY_PATH")
using RCall

R"""
library(fgsea)
library(BiocParallel)

# Create a SerialParam object
serial_param <- SerialParam()
"""

# 1-hot encoding of embryo conditions
conds = readcsv("data/conds.csv")
C = sparse(Matrix(conds)')
c,n = size(C)

# adjacency matrix for whether a pair of embryos corresponds to a
# protein interaction in STRINGdb
G_STRINGdb = readmat("data/G_STRINGdb.csv")
n_STRINGdb = sum(G_STRINGdb)

# labels for interactions in "protein->protein" format
labels = readmat("data/interactionlabels.csv")
ids = reshape(labels, c^2)

# labels for interactions in STRINGdb
# `fgsea` requires scores to be named
ids_sig = readcsv("data/stringdblabels.csv")[:, 1]
@rput ids ids_sig

# indices of labels in STRINGdb
i_stringdb = mapreduce(x->findall(ids .== x), vcat, ids_sig)
G_cond = Matrix{Int64}(undef, c, c)
G_cond .= 0
G_cond[i_stringdb] .= 1

function intxsub(G::AbstractMatrix, i, j)
    G[C[i, :], C[j, :]]
end

function mapintx(f, args...; kwargs...)
    score = mapreduce(hcat, 1:c) do i
        mapreduce(vcat, 1:c) do j
            f(args..., i, j; kwargs...)
        end
    end
end

# count edges between `names(conds)[i]` and `names(conds)[j]`
# in adjacency matrix `G`
# if `weighted`, 
function interactionscore(G::AbstractMatrix, i, j; scale = false)
    (sum ∘ intxsub)(G, i, j; scale = scale)
end

# if indices are not specified, calculates `interactionscore` in `G`
# for each pair of proteins
function interactionscore(G::AbstractMatrix;
                          scale = false, weighted = false)
    if !weighted
        G = G .> 0
    end
    score = mapintx(interactionscore, G; scale = scale)
    reshape(score, c^2)
end
# FFI for `fgsea` R function
function fgsea(scores::AbstractVector)
    @rput scores
    R"""
#scores[!is.finite(scores)] <- max(Filter(is.finite, scores))
names(scores) <- ids
res <- fgsea(list(score=ids_sig),scores,
             BPPARAM=serial_param,scoreType="pos")
"""
    @rget res
    return res
end

function fgsea(G::AbstractMatrix)
    scores = interactionscore(G)
    fgsea(scores)
end

# crossvalidation for adjacency matrix compared to an expected
# adjacency matrix
function cv_edge(G, G_0)
    G = G .> 0
    TP = sum(G .* G_0)
    precision = TP / sum(G)
    recall = TP / sum(G_0)
    return precision, recall
    #FP = sum(G .* (1 .- G_test))
    #TN = sum((1 .- G) .* (1 .- G_test))
    #FN = sum((1 .- G) .* G_test)
    #TP, FP, TN, FN
end

# precision & recall for an adjacency matrix compared to one predicted
# from STRINGdb interactions
function pr_stringdb(G; weighted = false)
    if !weighted
        G = G .> 0
    end
    TP = sum(G .* G_STRINGdb) 
    precision = TP / sum(G)
    recall = TP / n_STRINGdb
    return precision, recall
end

# native julia GSEA implementation
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

# loss functions
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

function loss_stringdb(G)
    precision, recall = pr_stringdb(G)
    #ES = es_stringdb(G)
    gsea = @suppress fgsea(G)
    ES = gsea.ES
    NES = gsea.NES
    ES, NES, precision, recall
end
