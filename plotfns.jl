using Plots
using Plots.PlotMeasures

function newplot(ylab::String, xlab::String = "batch", args...; kwargs...)
    scatter(xlabel = xlab, ylabel = ylab,
            margin=2mm,
            args...; kwargs...);
end

function points!(L::AbstractVector,
                 lab::String, args...; kwargs...)
    scatter!(1:length(L), L, label = lab, args...; kwargs...)
end

function points!(n::Integer,
                 L::AbstractVector,
                 lab::String,
                 args...; kwargs...)
    sel = sample(1:length(L), n)
    scatter!(sel, L[sel], label = lab, args...; kwargs...)
end

function points!(L)
    scatter!(1:length(L[1]), L[1], label = L[2],
             markershape = :cross, markersize = 0.5)
end

function plotloss(L::AbstractVector,
                  labs::AbstractArray{String},
                  ylab::String,
                  args...; markershape = :cross, markersize = 1,
                  path::String = "", file::String = "", kwargs...)
    p = newplot(ylab, args...; kwargs...);
    map(zip(L, labs)) do x
        points!(x[1], x[2];
                markershape = markershape, markersize = markersize);
    end
    if length(file) > 0
        savepath(p, path, file);
    else
        p;
    end
end

function plotloss(L::Dict,
                  labs::Dict,
                  ylab::String,
                  path::String,
                  file::String,
                  args...; markershape = :cross, markersize = 1,
                  kwargs...)
    p = newplot(ylab, args...; kwargs...);
    mapkey(L,labs) do y, lab
        points!(y, lab;
                markershape = markershape, markersize = markersize);
    end
    savepath(p, path, file);
end

function plotloss(n::Integer,
                  L::AbstractVector,
                  labs::AbstractArray{String},
                  ylab::String,
                  args...; markershape = :cross, markersize = 0.5, 
                  kwargs...)
    p = newplot(ylab,args...; kwargs...);
    map(zip(L, labs)) do x
        points!(n, x[1], x[2], markershape = markershape,
                markersize = markersize);
    end
end
    
function plotstat(Ls::AbstractVector, param::String, stat::String, 
                  labs::AbstractVector, args...; 
                  markershape = :cross, markersize = 5, kwargs...)
    xs = map(Ls, L->L[:, param], Ls)
    ys = map(Ls, L->L[:, stat], Ls)
    p = newplot(stat, param, args...; kwargs...);
    map(zip(xs, ys, labs)) do x, y, lab
        scatter!(x, y; label = lab, markersize = markersize, 
                 markershape = markershape)
    end
