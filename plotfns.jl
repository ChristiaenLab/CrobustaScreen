using Plots
#using Plots.PlotMeasures

function newplot(ylab::String, xlab::String = "batch", args...; kwargs...)
    scatter(xlabel = xlab, ylabel = ylab,
            #margin=10mm,
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
                  path::String,
                  file::String,
                  args...; kwargs...)
    p = newplot(ylab, args...; kwargs...);
    map(zip(L, labs)) do x
        points!(x[1], x[2], markershape = :cross, markersize = 1);
    end
    savepath(p, path, file);
end

function plotloss(L::Dict,
                  labs::Dict,
                  ylab::String,
                  path::String,
                  file::String,
                  args...;kwargs...)
    p = newplot(ylab, args...; kwargs...);
    mapkey(L,labs) do y, lab
        points!(y, lab);
    end
    savepath(p, path, file);
end

function plotloss(n::Integer,
                  L::AbstractVector,
                  labs::AbstractArray{String},
                  ylab::String,
                  path::String,
                  file::String,
                  args...; kwargs...)
    p = newplot(ylab,args...; kwargs...);
    map(zip(L, labs)) do x
        points!(n, x[1], x[2], markershape = :cross, markersize = 0.5);
    end
    savepath(p, path, file);
end
    
    
