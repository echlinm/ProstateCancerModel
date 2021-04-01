module SummaryAnalysis

include("MatrixImport.jl") #function for importing data
include("ExpressionCount.jl") # function for generating a histogram of expression counts per cell
include("RNA_HD.jl") # function for generating a histogram of Hamming distances between pairs of cells
using SparseArrays
#export

#function run()
## Load Data
#Load D17 data
data_dir = "GSE117403_D17_FACS\\GSE117403_D17_FACS"; #Specify a data directory
SeqVals_17 = MatrixImport.mtxImport(data_dir,"mtx"); #load data

#Load D27 data
data_dir = "GSE117403_D27_FACS\\GSE117403_D27_FACS"; #Specify a data directory
SeqVals_27 = MatrixImport.mtxImport(data_dir,"mtx"); #load data

#Load Pd data
data_dir = "GSE117403_Pd\\GSE117403_Pd"; #Specify a data directory
SeqVals_Pd = MatrixImport.mtxImport(data_dir,"mtx"); #load data

## Compare Data
# How many nonzero expression values does each cell contain?
ExpressionCount.expHist(SeqVals_17,SeqVals_27,SeqVals_Pd)
# What is the average and deviation of the Hamming distance between cells?
# Is this different in different patients?
RNA_HD.analysis([SeqVals_17 SeqVals_27],10000,"D17D27")
#end #function
end #module
