using Autoencoders, TrainingIO

###############################
# add command line arguments
using  ArgParse

s = ArgParseSettings()
@add_arg_table! s begin
    "--path", "-p"
        help = "Where to save the model"
        default = date()
    "--savecheckpts", "-c"
        help = "Export a model after each epoch"
        action = :store_true
end

args = parse_args(ARGS, s)

# where to save model
path = args["path"]
###############################


###############################
# define hyperparameters
epochs = 10000
batchsize=512

# learning rate
η = 0.0001

# weight decay rate
λ = 0.0001
optimiser = OptimiserChain(Flux.AdamW(η),Flux.WeightDecay(λ))

################################


################################
# load data
X = readmat("data/X.csv")'
m,n = size(X)
################################


################################
# define training set & test set; split training set into minibatches 
# fraction to use for test set
frac = 10
n_test = Integer(2^round(log2(n) - log2(frac)))
n_train = n - n_test

# split into training and test set
test,train = sampledat(X,n_test) |> gpu

# split into minibatches of batchsize
loader = Flux.DataLoader((train,train),batchsize=n_test,shuffle=true) 
################################


################################
# define loss function
# `loss(f)` converts a binary loss function `f::(X -> Y -> AbstractFloat)`
# to a function `g::((X -> Y) -> X -> Y -> AbstractFloat`
f = loss(Flux.mse)

# `loss_test(f, model, testdata, trainingdata)` returns a tuple `(L_train, L_test)`,
# where `L_train` is training set loss and `L_test` is test set loss.
lossfn = (M,x,y)->loss_test(f,M,test,x)
################################


################################
# define autoencoder

# build a dense MLP from a vector of layer widths
function mlp(l::AbstractVector{<: Integer},f::Function)
    θ = foldl(l[3:length(l)],
              init=Chain(Dense(l[1] => l[2],f))) do layers,d
        d_0 = size(layers[length(layers)].weight)[1]
        return Chain(layers...,Dense(d_0 => d,f))
    end
end

#layer widths
layers = [m, 58, 29, 14]

# define encoder & decoder; move to GPU
encoder = mlp(layers, tanh) |> gpu
decoder = mlp(reverse(layers), tanh) |> gpu

autoenc = Autoencoder(encoder,decoder)

################################


################################
# train autoencoder
# save model and loss table to path
# return loss table as `L`
L = train!(autoenc,
           path,
           lossfn,
           loader,
           optimiser,
           epochs,
           savecheckpts=args["savecheckpts"]) # if true, saves a separate model each epoch
################################


################################
# write embeddings to `data/E.csv`

# embeddings
E = encode(autoenc,gpu(X))

# move to CPU and write to file
writecsv(cpu(E'), path,"E.csv")
