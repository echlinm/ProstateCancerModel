# Generate_SCNS_Input takes sequence data, cytotrace data, and cluster data to generate the inputs for the SCNS command line program
# required scripts: MatrixImport.jl, Binarize.jl,DataSubset.jl,Adjacency.jl
# required files: CytoTRACE_plot_table, ClusterIDs, genelist
# output: _States.csv, _initial_states.txt, target_states.txt, _Edges.csv, Parameters.csv

module Generate_SCNS_Input

using DelimitedFiles

include("MatrixImport.jl")
include("Binarize.jl")
include("DataSubset.jl")
include("Adjacency.jl")

## Data specifiers
data_directory = "C:/Users/hcmoec/github/ProstateCancerModel"
data_source = "GSE117403_D17_FACS"
data_filepath = "$(data_directory)/$(data_source)/$(data_source)"
gene_filepath = "$(data_directory)/prostate_tfs.txt"

## States.csv
# a Cell x Genes Boolean matrix with header row and column for cell and gene names
# Import Genes x Cells SeqValue matrix
states = MatrixImport.seqvals(data_filepath,"mtx")
# Convert to a Cells x Genes SeqValue matrix
states = permutedims(states)
# Convert to Boolean Genes x Cells
states = Binarize.convert(states)
# Reduce Genes to desired subset
set_index = DataSubset.subset(data_filepath,gene_filepath,"gene") # row indices of desired subset
states = states[:,set_index]
# Reduce Cells to unique cells
unique_index = DataSubset.uniquestates(states,2)
states = states[unique_index,:]
# Create dense matrix with column and row names
SCNS_states = Array{Any,2}(undef,size(states)[1]+1,size(states)[2]+1) #initialize matrix
SCNS_states[1,1]=""
SCNS_states[2:end,2:end] = states # assing state values
SCNS_states[2:end,1] = MatrixImport.barcodes(data_filepath)[unique_index]
SCNS_states[1,2:end] = readdlm(gene_filepath)
#write file
writedlm("$(data_directory)/$(data_source)/$(data_source)_States.csv",SCNS_states,",")

## Initial_states.txt
# Use cytotrace values and clustering data to choose which cells are initial cells
# Clustering parameters used are highly specific to the particular data set
# Load cytotrace data matrix and keep barcodes and cytotrace values
cytotrace_values = readdlm("$(data_directory)/$(data_source)/$(data_source)_filtered_CytoTRACE_plot_table.txt")[2:end,1:2]
# Load cluster data
numdims= 10 # number of dimensions used in cluster analysis
rez = 0.5 # resolution used in cluster analysis
cluster_values = readdlm("$(data_directory)/$(data_source)/$(data_source)_filtered_ClusterIDs_d$(numdims)r$(rez).csv",',')[2:end,2]
# Select cell barcodes for input cells
cluster_id = 2 # Cluster Id chosen to select cells
cytotrace_threshold = 0.8 # Cytotrace value chosen to select cells
selected_barcodes =  cytotrace_values[(cytotrace_values[:,2].>= cytotrace_threshold) .& (cluster_values.==cluster_id),1]
# Load States.csv
states_barcodes = readdlm("$(data_directory)/$(data_source)/$(data_source)_States.csv",',')[2:end,1]
# Find selected barcodes in states barcodes
initial_barcodes = intersect(states_barcodes,selected_barcodes)
# Save input barcodes to file
writedlm("$(data_directory)/$(data_source)/$(data_source)_initial_states.txt",initial_barcodes)

## Target_states.txt
# use initial states to identify all other states in the state matrix as target states
# Load initial states
initial_barcodes = readdlm("$(data_directory)/$(data_source)/$(data_source)_initial_states.txt")
# Load States.csv
states_barcodes = readdlm("$(data_directory)/$(data_source)/$(data_source)_States.csv",',')[2:end,1]
# Identify all non-initial states
target_barcodes = setdiff(states_barcodes,initial_barcodes)
# Save target barcodes to file
writedlm("$(data_directory)/$(data_source)/$(data_source)_target_states.txt",target_barcodes)

## Edges.csv
# A list of all edges between states in States.csv. Each edge is a row with three columns, Gene, StateA, StateB
# Load States.csv
states = readdlm("$(data_directory)/$(data_source)/$(data_source)_States.csv",',')
# Create edge list with Gene, StateA, StateB names for each edge, in both directions
edge_list = Adjacency.edgelist(states)
# Save edge list to file
writedlm("$(data_directory)/$(data_source)/$(data_source)_Edges.csv",edge_list,",")

## Parameters.csv
# A list of MaxActivators,MaxRepressors, and Threshold values for each gene with gene names
# MaxActivators = MaxRepressors=3, Threshold=80
# Load Gene list
gene_list = readdlm(gene_filepath,',')
# Initialize parameter_values
parameter_values = Array{Any,2}(undef,length(gene_list)+1,4)
# Add headers
parameter_values[1,:] = ["Gene", "MaxActivators","MaxRepressors","Threshold"]
# Add gene names
parameter_values[2:end,1] = gene_list
# Add parameter values
parameter_values[2:end,2:end] = repeat([3 3 80],length(gene_list))
# Save to file
writedlm("$(data_directory)/$(data_source)/$(data_source)_Parameters.csv",parameter_values,",")

end #module
