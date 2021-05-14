module BNrules

using JSON

# get variable names
function getnames(model,ng)
    GeneNames=Vector{String}(undef,ng)
    for i = 1:ng
        GeneNames[i]= model["Model"]["Variables"][i]["Name"];
    end
    GeneNames
end
# getnv takes a model dictionary in and outputs the in-degree for each variable
function getnv(model,nv,ng)
    #add initializations
    for i=1:ng
        funS= model["Model"]["Variables"][i]["Formula"];
        stIdx=maximum.(findall("r",funS)).+2;
        endIdx=maximum.(findnext.(")",funS,stIdx)).-1;
        k=String[]
        for j = 1:length(stIdx)
            push!(k,funS[stIdx[j]:endIdx[j]])
        end
        nv[i]=length(unique(k))
    end
end #function

function getvarF(model,varF,nv,ng)
    #add initializations
    for i=1:ng
        funS= model["Model"]["Variables"][i]["Formula"];
        stIdx=maximum.(findall("r",funS)).+2;
        endIdx=maximum.(findnext.(")",funS,stIdx)).-1;
        k=String[]
        for j = 1:length(stIdx)
            push!(k,funS[stIdx[j]:endIdx[j]])
        end
        varF[1:nv[i],i]=sort(Meta.parse.(unique(k))).+1
    end
end #function

function getF(model,F,varF,nv,ng)
    global var=Float64[]
    funS=String[]
    funE=Expr[]
    for i=1:ng # for each gene
        var=zeros(ng)
        funS= model["Model"]["Variables"][i]["Formula"]; # update function as a string
        for j = varF[1:nv[i],i] # for each input 1:nv[i]
            funS = replace(funS, "var($(j-1))"=>"var[$(j)]") # replace () with [] indexing and 0-based with 1-based numbering for julia
        end
        funE=Meta.parse(funS); # evaluable expression form of the ith update function
        for j=0:2^(nv[i])-1
            var[varF[1:nv[i],i]] = reverse(digits(j,base=2,pad=nv[i])) #
            F[j+1,i]=eval(funE)
        end #for j
    end #for i
end #function
# generate takes a model file in and outputs the BN matrices corresponding to it
# nv - a ng length vector with the in-degree of each variable
# varF - a kmax x ng matrix in which each column gives the inputs of each variable(sorted low to high)
# F - a 2^kmax x ng matrix in which each column gives the outputs of each variables to all of its inputs, [k1 k2...kmax] (sorted lexicographically)
function getNetwork(model_path)
    model=JSON.parsefile(model_path); #dictionary of model
    ng=size(model["Model"]["Variables"],1); # number of genes in model
    nv=fill(0,ng); getnv(model,nv,ng)# in-degree of each gene
    kmax=maximum(nv); # max in-degree
    varF=fill(-1,kmax,ng); getvarF(model,varF,nv,ng);#initialized varF matrix
    F=fill(-1,2^kmax,ng); getF(model,F,varF,nv,ng); # initialized F matrix
    GeneNames=getnames(model,ng);
    GeneNames,nv,varF,F
end #function

end #module
