module CellMarkers

using DelimitedFiles


function select(Seqvals,markers,filepath,type,thresh::Int64=0)
    SeqVals=SeqVals_17.>thresh;
    #if type == "mtx"
    Genes = readdlm("$(filepath)_filtered_genes.tsv")[:,2]
    #end
        for i = 1:size(M,1)
            geneState = Int64(M[i,2] .== "pos")
            geneIndx=findall(Genes.==M[i,1])
            cellIndx=findnz(SeqVals[geneIndx,:] .== geneState)[3]
            SeqVals = SeqVals[:,cellIndx]
        end
    SeqVals
end#function

#PTPRC=CD45, PECAM1 = CD31, EPCAM = CD326, DPP4=CD26, NGFR = CD271
M = ["PTPRC" "neg";
    "PECAM1" "neg";
    "EPCAM" "neg"];
T=0
SelectVals= select(SeqVals_17,M,data_dir,"mtx", T)
end#module
