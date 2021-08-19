#DataSubset provides the indices for a given subset of genes in a given data matrix
module DataSubset

using DelimitedFiles

include("MatrixImport.jl") #for MatrixImport.genes

# subset takes two directories, one for the full data and one for the set, and outputs the indices for that gene set
function subset(data_dir,set_dir,set_type)
    if set_type == "gene" # if taking a subset of genes
        SetNames=MatrixImport.genes(data_dir); #Matrix of all names in the full data set
    elseif set_type == "cell" # if taking a subset of genes
        SetNames=MatrixImport.barcodes(data_dir); #Matrix of all names in the full data set
    end
    if isfile(set_dir) #if the subset list is in file form
        SubNames = readdlm(set_dir); #Matrix of all names in the desired subset
        Matches = SetNames[:,end] .== reshape(SubNames,1,length(SubNames)); #Compare the two sets to find the indices of the subset in the full set
        CartIdx=findall(Matches); #The Cartesian indices of the subset
        SetIdx = repeat([1:size(SetNames,1);],1,length(SubNames))[CartIdx] #The Linear Indices of the subset
    else
        Matches = occursin.(set_dir, SetNames)[:]; #BitVector of all elements in SetNames that contain set_dir
        SetIdx = (1:size(SetNames,1))[Matches]; #The Linear Indices of the subset
    end
    SetIdx
end#function

function uniquestates(state_mtx,dim)
    b = 2 .^[size(state_mtx)[dim]-1:-1:0;]
    if dim == 1
        integer_states = permutedims(state_mtx) * b
        unique_states = unique(integer_states)
        unique_idx = indexin(unique_states,integer_states)
    elseif dim == 2
        integer_states = state_mtx * b
        unique_states = unique(integer_states)
        unique_idx = indexin(unique_states,integer_states)
    end
    unique_idx
end #function

end#module
