# MatrixImport is for defining functions for importing data in different file formats

module MatrixImport

using MatrixMarket, DelimitedFiles

# mtxImport imports data of different formats (MatrixMarket-mtx,...)
function mtxImport(filepath,type)
    if type =="mtx"  # MatrixMarket data format
        MatrixMarket.mmread("$(filepath)_filtered_matrix.mtx")
    end #if
end #function

#mtxLoad imports RNASeq data, Gene Name, and Barcode files
function mtxLoad(filepath,type)
    if type == "mtx"
        mtxImport(filepath,"mtx"), readdlm("$(filepath)_filtered_genes.tsv"), readdlm("$(filepath)_filtered_barcodes.tsv")
    end #if
end #function

end #module
