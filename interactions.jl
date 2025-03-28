using JLD2, ProgressMeter, Statistics
using SparseArrays, OneHotArrays
using TrainingIO

interactions = readcsv("data/interactions.csv")
groups = readcsv("data/groups.csv")

labels = groups.Condition
conds = unique(labels)

n = size(groups,1)
c = length(conds)

rconds = map(x->replace(x,"/" => ".",
                        "FZD" => "Fzd",
                        "GNA-12" => "Gna12",
                        "Gna-L" => "GNA.L"),
             conds)
labelkey = DataFrame(hcat(1:c,conds,rconds),[:index,:label,:rlabel])
writecsv(labelkey,"data/conditionkey.csv")

label_int = mapreduce(hcat,1:c) do i
    mapreduce(vcat,1:c) do j
        rconds[i]*"->"*rconds[j]
    end
end
writecsv(label_int,"data/interactionlabels.csv")

from = mapreduce(i->findall(i .== rconds),vcat,interactions."name.x")
to = mapreduce(i->findall(i .== rconds),vcat,interactions."name.y")
interactions[:,:id_from] .= from
interactions[:,:id_to] = to
writecsv(interactions,"data/interactions.csv")

label_stringdb = map((i,j)->label_int[i,j],to,from)
writecsv(label_stringdb,"data/stringdblabels.csv")

#C = mapreduce(c->labels .== c,hcat,conds);
C = onehotbatch(labels,conds)
df = DataFrame(C',rconds)
writecsv(df,"data/conds.csv")

fill_i = I->begin
    n = size(I,1)
    M = zeros((n,n))
    M[I,:] .= 1
    M
end

v_G = map((i,j)->fill_i(C[i,:]) .* fill_i(C[j,:])',from,to)
G_STRINGdb = Bool.(sum(v_G))
writecsv(G_STRINGdb,"data/G_STRINGdb.csv")

sel = map(i->BitArray(i * ones(Bool,1,n)),eachslice(C,dims=2));
edges = map(i->map(j->i .* j',sel),sel);

mask_GSEA = sparse(from,to,true,c,c)
sel_GSEA = reshape(mask_GSEA,c^2)

mask_int = map((i,j)->edges[i][j],from,to)
G_int = foldl(.|,tmp)
writecsv(G_int,"data/G_interactions.csv")

