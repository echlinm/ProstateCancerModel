# MatrixImport is for defining functions for importing data in different file formats

module MatrixImport

using MatrixMarket, DelimitedFiles

#seqvals imports sequence data for different format data structures (MatrixMarket-mtx,...)
function seqvals(filepath,type)
    if type =="mtx"  # MatrixMarket data format
        MatrixMarket.mmread("$(filepath)_filtered_matrix.mtx")
    end #if
end #function

#genes imports Gene Names for different format data structures
function genes(filepath)
    readdlm("$(filepath)_filtered_genes.tsv")
end #function

#genes imports Gene Names for different format data structures
function barcodes(filepath)
    readdlm("$(filepath)_filtered_barcodes.tsv")
end #function

#loadall imports Seq values, Gene Names, and Barcodes for different format data structures
function loadall(filepath,type)
    if type == "mtx"
        MatrixMarket.mmread("$(filepath)_filtered_matrix.mtx"), readdlm("$(filepath)_filtered_genes.tsv"), readdlm("$(filepath)_filtered_barcodes.tsv")
    end #if
end #function

end #module
