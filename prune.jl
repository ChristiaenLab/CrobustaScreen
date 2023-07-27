using NaiveNASflux, Flux, Test
using Flux: train!, mse
import Random

using CSV, Flux, Statistics, ProgressMeter, DataFrames, InvertedIndices
using StatsBase
using NearestNeighbors, TensorCast, TensorOperations

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

dat = DataFrame(CSV.File("out/z_dat.csv",normalizenames=true))
X = Matrix(dat[:,2:end])
X = scaledat(X)
test,train = sampledat(X,10)
i = size(X)[1]
n = i

θ = CompGraph(fluxvertex(Dense(i => n,relu)),
              fluxvertex(Dense(n => i))) |> gpu

Random.seed!(0)
niters = 50

densevertex(in, outsize, act) = fluxvertex(Dense(nout(in),outsize, act), in, layerfun=ActivationContribution)
invertex = denseinputvertex("input", i) |> gpu
layer1 = densevertex(invertex, n, relu)
layer2 = fluxvertex(Dense(n, i),layer1)
original = CompGraph(invertex, layer2) |> gpu

# Training params, nothing to see here
opt = AdamW(0.1) |> gpu
loss(g) = (x, y) -> mse(g(x), y)

# Training data: xor truth table: y = xor(x) just so we don't need to download a dataset.
x = Float32[0 0 1 1;
            0 1 0 1]
y = Float32[0 1 1 0]

# Train the model
η = 0.01
batchsize = 256
loader = Flux.DataLoader((X, X), batchsize=batchsize, shuffle=true) |> gpu
optim = Flux.setup(Flux.AdamW(η), original)  # will store optimiser momentum, etc.
train!(loss(original), Flux.params(original), loader, opt)
@test loss(original)(x, y) < 0.001

nprune = 16

pruned_least = deepcopy(original)
Δnout!(pruned_least[2] => -nprune)

pruned_most = deepcopy(original)
Δnout!(pruned_most[2] => -nprune) do v
    vals = NaiveNASlib.defaultutility(v)
    return 2*sum(vals) .- vals # Ensure all values are still > 0, even for last vertex
end

pruned_random = deepcopy(original)
Δnout!(v -> rand(nout(v)), pruned_random[2] => -nprune)

@test   loss(pruned_most)(x, y)   >
        loss(pruned_random)(x, y) >
        loss(pruned_least)(x, y)  >=
        loss(original)(x, y)
@test loss(pruned_least)(x, y) ≈ loss(original)(x, y) atol = 1e-5
