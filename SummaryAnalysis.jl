module SummaryAnalysis

using MatrixMarket
#export

#Define all functions here
#Function for importing data of different formats
function mtxImport(filename,type)
    if type =="mtx"  # MatrixMarket data format
        M=MatrixMarket.mmread(filename)
    end
end

filename = "C:\\Users\\hcmoec\\Documents\\Data_Files\\GSE117403_D17_FACS_filtered_matrix.mtx\\GSE117403_D17_FACS_filtered_matrix.mtx"

@btime M=mtxImport(filename,"mtx")

end
