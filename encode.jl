using Flux: train!
using CSV, Flux, Statistics, ProgressMeter, DataFrames, InvertedIndices
using StatsBase
using NearestNeighbors, TensorCast, TensorOperations

struct KStep
    Kn::Matrix
    MSE::Float64
end

mutable struct KSteps
    vals::Matrix
    K::Matrix
    steps::AbstractArray{KStep,1}

    function init(vals,K)
        E = matmul(K,vals)
        MSE = mse(vals,E)
        return KSteps(vals,K,[KStep(K,MSE)])
    end
end

function kinit(vals,K)
    E = matmul(K,vals)
    MSE = mse(vals,E)
    return KSteps(vals,K,[KStep(K,MSE)])
end


#function mse(X,E)
#    return mean((X-E).^2)
#end

function aic(nvars,loss)
    return 2 .* nvars .- 2 .* log2.(1 .- loss)
end

function matmul(K,X)
    @tensor E[k,j] := K[i,j] * X[k,i]
    return E
end

function colsum(X)
    @reduce W[j] := sum(i) X[i,j]
    return W |> cpu
end


function kstep(steps::KSteps)
    Kn = matmul(steps.steps[length(steps.steps)].Kn,steps.K)
    E = matmul(Kn,steps.vals)
    MSE = mse(steps.vals,E)
    return KSteps(steps.vals,K,cat(steps.steps,KStep(Kn,MSE),dims=1))
end

function nsteps(n)
    return eval(reduce((x,y) -> Expr(:call,:∘,x,y),repeat([:kstep],n)))
    #return eval(reduce((x,y) -> Expr(:call,:∘,x,y),repeat([:kstep],n)))(kinit(embedding,K))
end

function scaledat(X)
    Y = mapslices(maximum,abs.(X),dims=1)[1,:]
    X = eachslice(X,dims=2)./Y
    X = X[Y.!=0]
    X = transpose(hcat(X...))
    return X
end

function sampledat(X,frac)
    j = size(X)[2]
    sel = sample(1:j,j ÷ frac)
    test = X[:,sel]
    train = X[:,Not(sel)]
    return test,train
end

function encoder(X,n,epochs,batchsize,η=0.01)
    i = size(X)[1]
    model = Chain(
        Dense(i => n, relu),   # activation function inside layer
        #BatchNorm(n),
        Dense(n => i)) |> gpu        # move model to GPU, if available
    
    #out1 = model(X |> gpu) |> cpu
    loader = Flux.DataLoader((X, X) |> gpu, batchsize=batchsize, shuffle=true);
    optim = Flux.setup(Flux.AdamW(η), model)  # will store optimiser momentum, etc.

    # Training loop, using the whole data set 1000 times:
    losses = []
    @showprogress for epoch in 1:epochs
        for (x, y) in loader
            loss, grads = Flux.withgradient(model) do m
                # Evaluate model and loss inside gradient context:
                y_hat = m(x)
                Flux.mse(y_hat, y)
            end
            Flux.update!(optim, model, grads[1])
            push!(losses, loss)  # logging, outside gradient context
        end
    end
    return model
end

function knnmat(tree,embedding,k)
    j = length(tree.data)
    g = knn(tree,embedding, k)

    function f(i,d)
        x = repeat([0.0],j)
        x[i] = d
        return x
    end

    K = hcat(map(f,g[1],g[2])...)
    return K
end

function knnkern(K)
    K = -log.(K)
    K[Not(isfinite.(K))] .= 0

    #@reduce W[j] := sum(i) K[i,j]
    W = colsum(K)
    K = transpose(K) ./ W
    #@reduce V[i] := sum(j) K[i,j]
    V = colsum(K)
    all(isapprox.(V,1))
    return K
end

struct PrunableDense
  dense :: Dense
  mask:: BitMatrix
end

function (a:: PrunableDense)(x:: AbstractVecOrMat)
  W, b,sigma , M = a.dense.W, a.dense.b, a.dense.sigma , a.mask
  return sigma.((W.*M)*x .+ b)
end


dat = DataFrame(CSV.File("out/z_dat.csv",normalizenames=true))
X = Matrix(dat[:,2:end])
X = scaledat(X)
test,train = sampledat(X,10)

model = encoder(train,114,15000,128)

encoders = map(x -> encoder(train,x,1000,128),1:10)

encMSE = map(x -> mse(test,(x(test |> gpu) |> cpu)),encoders)
encAIC = aic(collect(1:10),encMSE)

embedding = encoders[1][1](X |> gpu) |> cpu

tree = KDTree(embedding)
knns = map(k -> knnmat(tree,embedding,k),3:100)
Ks = map(knnkern,knns)

kMSE = map(K -> mse(embedding,matmul(K,embedding)),Ks)   
sel = argmin(kMSE)
K = Ks[argmin(kMSE)]
steps = eval(reduce((x,y) -> Expr(:call,:∘,x,y),repeat([:kstep],100)))(kinit(embedding,K))
argmin(map((x->x.MSE),steps.steps))
argmin(map((n,x)->aic(n,x.MSE),1:length(steps.steps),steps.steps))

ksteps = map(K->nsteps(5)(kinit(embedding,K)),Ks)
kMSEs = map(J->map(x->x.MSE,J.steps),ksteps)
kMSEs = map(i->map(x->x[i],kMSEs),1:length(kMSEs[1]))
p = scatter(3:100,kMSEs,xlabel=L"k",ylabel=L"\mathrm{MSE}(K_{knn}X)")
savefig(p,"kMSE.pdf")

using Leiden, Distributions, LinearAlgebra

function clust(K,res)
    return Leiden.leiden((K.+transpose(K)) .> 0,resolution=res)
end

function clustmat(H)
    n = maximum(vcat(H.partition...))
    M = zeros(n,n)
    for P in H.partition
        for i in P
            M[i,P] .= 1
        end
    end
    M = M .* ((Matrix(I,n,n) .- 1) .* (-1))
    return M ./ colsum(M)
end

γs = rand(Uniform(0.001,0.2),100)
clusts = map((γ -> clust(knns[sel], γ)), γs)

Js = map(clustmat,clusts)

Jstep = map(J->nsteps(5)(kinit(embedding,J)),Js)
JMSE = map(J->map(x->x.MSE,J.steps),Jstep)

using Plots, LaTeXStrings

MSEs = map(i->map(x->x[i],JMSE),1:length(JMSE[1]))
labels = map(n->L"n="*string(n),1:6)
    
p = scatter(γs,MSEs,label=reshape(labels,1,length(MSEs)),xlabel=L"γ",ylabel=L"\mathrm{MSE}(M^{n}_{clust}X)")
savefig(p,"clustMSE.pdf")

