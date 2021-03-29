module SummaryAnalysis

using MatrixMarket,SparseArrays
#export

#Define all functions here
#Function for importing data of different formats (MatrixMarket-mtx,...)
function mtxImport(filename,type)
    if type =="mtx"  # MatrixMarket data format
        MatrixMarket.mmread(filename)
    end
end

#Specify a filename to be read
filename = "C:\\Users\\hcmoec\\Documents\\Data_Files\\GSE117403_D17_FACS_filtered_matrix.mtx\\GSE117403_D17_FACS_filtered_matrix.mtx"

#import RNASeq data
SeqVals = mtxImport(filename,"mtx")
end
