module Adjacency


# This is a quick and dirty way of doing this.
# Creates the table with geneID,cell1ID, and cell2ID for all hamming neighbors in a cellxgene matrix
function edgelist(state_list)
    state_mtx = state_list[2:end,2:end] #remove cell and gene labels
    AM = adjmatrix(state_mtx) # Get adjacency matrix
    num_neighbors = sum(AM) #Get number of hamming neighbors
    neighbor_idx = findall(AM.==1) #Get coordiates for each hamming neighbor, includes duplicates
    edge_list = Array{String,2}(undef,num_neighbors+1,3) #initialize list for edges
    edge_list[1,:]=["Gene","StateA","StateB"] # first row labels
    for i = 1:num_neighbors #for each hamming neighbor
        pair_id = neighbor_idx[i] # get index of each cell
        edge_list[i+1,2:3]= state_list[[pair_id[2] ; pair_id[1]].+1,1] # fill edge_list with names of each cell in pair
        gene_idx = findall(xor.(state_mtx[pair_id[1],:],state_mtx[pair_id[2],:]))[1] #get index of gene that differs
        edge_list[i+1,1]=state_list[1,gene_idx+1] # fill edge_list with name of gene that differs
    end #for
    edge_list #output edge_list
end #function






# Creates the hamming neighbor adjacency matrix from a Cell x Gene state matrix
function adjmatrix(state_mtx)
    C = size(state_mtx,1); #Number of cells, number of genes
    AM=fill(0,C,C) #Initialized adjacency matrix
    for i=1:C
        hamming_mtx = xor.(state_mtx[i,:]',state_mtx) # Get hamming difference between cell i and all other cells
        AM[i,:] = (sum(hamming_mtx,dims=2).==1) # Record whether cell i is hamming neighbors with all other cells
    end #for
    AM # output adjacency matrix
end # function



end #module
