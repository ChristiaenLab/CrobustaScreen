using Autoencoders, TrainingIO, DataFrames, CSV

dat = (DataFrame ∘ CSV.File)("data/z_dat.csv", normalizenames = true);
# dat = dat[:, Not(contains.(names(dat), "Threshold_"))]
dat = Matrix(dat[:, 2:end]);
dat = hcat(filter(x -> sum(x) != 0, eachslice(dat, dims=2))...);

X = scaledat(dat')
writecsv(X', "data/", "X.csv")
