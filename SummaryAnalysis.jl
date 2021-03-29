module SummaryAnalysis

using MatrixMarket,SparseArrays,DelimitedFiles,Plots,StatsBase
#export

#Define all functions here
#Function for importing data of different formats (MatrixMarket-mtx,...)
function mtxImport(filepath,type)
    if type =="mtx"  # MatrixMarket data format
        MatrixMarket.mmread("$(filepath)_filtered_matrix.mtx")
    end
end

#Specify a data directory
data_dir = "GSE117403_D27_FACS\\GSE117403_D27_FACS"

#Import RNASeq data, Gene Names, and Barcodes
SeqVals = mtxImport(data_dir,"mtx")
GeneNames = readdlm("$(data_dir)_filtered_genes.tsv")
BarCodes = readdlm("$(data_dir)_filtered_barcodes.tsv")

# How many nonzero expression values does each cell contain?
CountPerCell = sum(SeqVals.!=0,dims=1) #number of nonzeros per cell barcode
Count_hist=fit(Histogram,CountPerCell[:],[0:250:4500;]) # histogram of the number of nonzeros per cell barcode
bar(Count_hist.weights)
end
