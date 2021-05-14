module SummaryAnalysis

include("MatrixImport.jl") #function for importing data
include("ExpressionCount.jl") # function for generating a histogram of expression counts per cell
include("RNA_HD.jl") # function for generating a histogram of Hamming distances between pairs of cells
include("DataSubset.jl") # function for generating the indices for a desired gene subset in the full data set
using SparseArrays
#export

#function run()
## Load Full Data
#Load D17 data
data_dir = "GSE117403_D17_FACS\\GSE117403_D17_FACS"; #Specify a data directory
SeqVals_17 = MatrixImport.seqvals(data_dir,"mtx"); #load data

#Load D27 data
data_dir = "GSE117403_D27_FACS\\GSE117403_D27_FACS"; #Specify a data directory
SeqVals_27 = MatrixImport.seqvals(data_dir,"mtx"); #load data

#Load Pd data
data_dir = "GSE117403_Pd\\GSE117403_Pd"; #Specify a data directory
SeqVals_Pd = MatrixImport.seqvals(data_dir,"mtx"); #load data

#Data subset
set_dir = "prostate_tfs.txt"
SetIdx = DataSubset.subset(data_dir,set_dir);
SeqVals_17_tfs = SeqVals_17[SetIdx,:]; #take subset
SeqVals_27_tfs = SeqVals_27[SetIdx,:]; #take subset
SeqVals_Pd_tfs = SeqVals_Pd[SetIdx,:]; #take subset

#Data subset
SetIdx = DataSubset.subset("GSE117403_D17_FACS\\GSE117403_D17_FACS","1","cell");
D17_BE_tfs = SeqVals_17_tfs[:,SetIdx]; #take subset
SetIdx = DataSubset.subset("GSE117403_D17_FACS\\GSE117403_D17_FACS","3","cell");
D17_LE_tfs = SeqVals_17_tfs[:,SetIdx]; #take subset
SetIdx = DataSubset.subset("GSE117403_D17_FACS\\GSE117403_D17_FACS","2","cell");
D17_OE_tfs = SeqVals_17_tfs[:,SetIdx]; #take subset
SetIdx = DataSubset.subset("GSE117403_D17_FACS\\GSE117403_D17_FACS","4","cell");
D17_FMST_tfs = SeqVals_17_tfs[:,SetIdx]; #take subset
SetIdx = DataSubset.subset("GSE117403_D27_FACS\\GSE117403_D27_FACS","1","cell");
D27_BE_tfs = SeqVals_27_tfs[:,SetIdx]; #take subset
SetIdx = DataSubset.subset("GSE117403_D27_FACS\\GSE117403_D27_FACS","2","cell");
D27_LE_tfs = SeqVals_27_tfs[:,SetIdx]; #take subset
SetIdx = DataSubset.subset("GSE117403_D27_FACS\\GSE117403_D27_FACS","3","cell");
D27_OE_tfs = SeqVals_27_tfs[:,SetIdx]; #take subset
SetIdx = DataSubset.subset("GSE117403_D27_FACS\\GSE117403_D27_FACS","4","cell");
D27_Edn_tfs = SeqVals_27_tfs[:,SetIdx]; #take subset
SetIdx = DataSubset.subset("GSE117403_D27_FACS\\GSE117403_D27_FACS","5","cell");
D27_StPDPNn_tfs = SeqVals_27_tfs[:,SetIdx]; #take subset
SetIdx = DataSubset.subset("GSE117403_D27_FACS\\GSE117403_D27_FACS","6","cell");
D27_StPDPNp_tfs = SeqVals_27_tfs[:,SetIdx]; #take subset


## Compare Data
# How many nonzero expression values does each cell contain?
ExpressionCount.expHist(SeqVals_17,SeqVals_27,SeqVals_Pd)

# What is the average and deviation of the Hamming distance between cells?
# Is this different in different patients?
set_dir = "prostate_tfs.txt"
#Data = SeqVals_17_tfs ; Source = "D17tfs";
#Data = SeqVals_27_tfs ; Source = "D27tfs";
#Data = SeqVals_Pd_tfs ; Source = "Pdtfs";
#Data = [SeqVals_17_tfs SeqVals_27_tfs SeqVals_Pd_tfs] ; Source = "AllCellstfs";
#Data = [D17_BE_tfs D27_BE_tfs] ; Source = "BEtfs";
#Data = [D17_LE_tfs D27_LE_tfs] ; Source = "LEtfs";
#Data = [D17_OE_tfs D27_OE_tfs] ; Source = "OEtfs";
#Data = D17_FMST_tfs ; Source = "FMSTtfs";
Data = D27_Edn_tfs ; Source = "Edntfs";
#Data = D27_StPDPNn_tfs ; Source = "StPDPNntfs";
#Data = D27_StPDPNp_tfs ; Source = "StPDPNptfs";
RNA_HD.HDhist(Data,10^6,Source)
RNA_HD.HDhistgene(Data,10^6,Source)
RNA_HD.HDbargene(Data,10^6,set_dir,Source)
RNA_HD.HNhist(Data,Source)
RNA_HD.HNbargene(Data,set_dir,Source)

# How many Hamming neighbors are there in the data?

#end #function
end #module
