module SCNS

G,C = size(SeqVals)
OutputMat = Array{Any,2}(undef,C+1,G+2);
OutputMat[2:end,3:end] .= SeqVals';

2 functions
one for full datat
-    transpose matrix
    add text column for cells
    add text column for gene names
    add text column for initial and non-initial cells
        this will probably come from another file or script
    export as csv
one for subset data



end #module
