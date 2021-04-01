module RNA_HD

using StatsBase,StatsPlots,Statistics
# HD will create a Hamming distance vector of N randomly sampled pairs of cells, without replacement.
function HD(SeqVals,N)
    SeqVals .= SeqVals.>0 #Binarize SeqVals
    G,C = size(SeqVals) #Number of genes and cells
    Seq_HD = fill(0.0,N) #Initialized vector of Hamming Distance between cells
    PairIJ = reshape(sample([1:C;],N*2,replace=true),N,2) # Cell pair indices
    for i = 1:N
        Seq_HD[i] = sum(xor.(SeqVals[:,PairIJ[i,1]],SeqVals[:,PairIJ[i,2]])) / G # Compare Cell 1 to Cell 2,
    end #for
    Seq_HD
end #function

function analysis(SeqVals,N,source)
    SeqHD = HD(SeqVals,N) # get HD values
    fig = histogram(SeqHD,  #generate histogram of values
        bins = [0:0.01:1;],
        xlim = (0,maximum(SeqHD)+0.02),
        ylim = (0,N/2),
        xlabel = "Hamming Distance Between $(N) Sampled Cell Pairs",
        #ylabel = "",
        #xticks = ([0:0:0;],string.([0:0:0;])),
        #yticks = ([0:0:0;],string.([0:0:0;])),
        xtickfont=font(10),
        ytickfont=font(10),
        guidefont=font(12),
        legendfont = font(10),
        legend = false,
        #label = [""],
        title = "Henry et al. Single Cell Data",
        framestyle=:box)
    savefig(fig,"Henryetal_CellHD_$(source)N$(N).png"); #savefig
end #function

end#module
