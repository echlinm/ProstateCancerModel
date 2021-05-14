module ExpressionCount

using StatsBase,StatsPlots

function expHist(SeqVals_17,SeqVals_27,SeqVals_Pd)
    CountPerCell = sum.((SeqVals_17.!=0,SeqVals_27.!=0,SeqVals_Pd.!=0),dims=1); #number of nonzeros per cell barcode
    Count_hist_17=fit(Histogram,CountPerCell[1][:],[0:500:7500;]); # histogram of the number of nonzeros per cell barcode
    Count_hist_27=fit(Histogram,CountPerCell[2][:],[0:500:7500;]); # histogram of the number of nonzeros per cell barcode
    Count_hist_Pd=fit(Histogram,CountPerCell[3][:],[0:500:7500;]); # histogram of the number of nonzeros per cell barcode
    #Grouped (by donor) histogram of cells by the number of expressed genes
    fig=groupedbar([Count_hist_17.weights Count_hist_27.weights Count_hist_Pd.weights],
        xlim = (0,17),
        ylim = (0,2*10^4),
        xlabel = "Number of Reported Genes Expressed (1x10³)",
        ylabel = "Cell Counts (1x10³)",
        xticks = ([0:2:17;],string.([0:1:17/2;])),
        yticks = ([0:5*10^3:2*10^4;],string.([0:5:20;])),
        xtickfont=font(10),
        ytickfont=font(10),
        guidefont=font(12),
        legendfont = font(10),
        label = ["Donor 17" "Donor 27" "D17/27/35"],
        title = "Henry et al. Single Cell Data",
        framestyle=:box)
    savefig(fig,"Henryetal_ExpressedGenesHist.png");
end #function
end #module
