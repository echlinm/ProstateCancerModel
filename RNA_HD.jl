# RNA_HD calculates Hamming distances between cells of a given set and outputs plots
# HD - creates a HD vector of randomly sampled pairs of cells

module RNA_HD

using StatsBase,StatsPlots,Statistics,SparseArrays,LaTeXStrings
using LinearAlgebra, DelimitedFiles, Plots.PlotMeasures

# HD will create a Hamming distance vector of N randomly sampled pairs of cells, with replacement. Gene = false will output plain HD. Gene = true will output HD broken into the contibuting genes.
function HD(SeqVals,N,Gene)
    SeqVals .= SeqVals.>0 #Binarize SeqVals
    G,C = size(SeqVals) #Number of genes and cells
    if Gene == false # SeqHD output
        Seq_HD = fill(0.0,N) #Initialized vector of Hamming Distance between cells
    else # GeneHD output
        Gene_HD = fill(0.0,G,G+1) # Initialized matrix for counting the occurrence of each gene in each Hamming Distance, Gene x Distance
    end
    PairIJ = reshape(sample([1:C;],N*2,replace=true),N,2) # Sampled cell pair indices
    for i = 1:N
        Seq_xor = xor.(SeqVals[:,PairIJ[i,1]],SeqVals[:,PairIJ[i,2]]) # Compare Cell i to Cell j
        if Gene == false #SeqHD output
            Seq_HD[i] = sum(Seq_xor) / G # Add all differences and normalize
        else #GeneHD output
            Gene_HD[findnz(Seq_xor)[1],sum(Seq_xor)+1] .+= 1/sum(Seq_xor) # Normalize each genes contribution and record
        end
    end #for
    if Gene == false #SeqHD output
        Seq_HD
    else # GeneHD output
        Gene_HD
    end
end #function


# hist generates a histogram of HDs from N cell pairs in given data
function HDhist(SeqVals,N,source)
    SeqHD = HD(SeqVals,N,false) # get HD values
    fig = histogram(SeqHD,  #generate histogram of values
        bins = [0:0.025:1;],
        xlim = (0,maximum(SeqHD)+0.02),
        ylim = (0,N/5),
        xlabel = "Hamming Distance",
        #ylabel = "",
        xticks = ([0:0.1:maximum(SeqHD)+0.02;],string.([0:0.1:maximum(SeqHD)+0.02;])),
        #yticks = ([0:0:0;],string.([0:0:0;])),
        xtickfont=font(10),
        ytickfont=font(10),
        guidefont=font(12),
        legendfont = font(10),
        legend = false,
        #label = [""],
        title = "Hamming Distance Between 10^($(Int64(log10(N)))) Sampled Cell Pairs \\n $(source)",
        framestyle=:box,
        top_margin = 3mm,
        size = (600,400))
    savefig(fig,"Henryetal_CellHD_$(source).png"); #savefig
end #function

# histgene generates a histogram of HDs from N cell pairs in given data, broken into the frequency of contributing genes
function HDhistgene(SeqVals,N,source)
    GeneHD = HD(SeqVals,N,true) # get HD values
    G = size(GeneHD,1) # number of genes
    lastcol = findlast(GeneHD.!=0)[2] # last nonzero HD
    cmap = distinguishable_colors(40) # colormap for stacked bar plot
    annot = string.(sum(GeneHD[:,2:lastcol].!=0,dims=1)[:]); # number of contributing genes in each HD
    a_dist = N*(1/50) # height above histogram column for annotation
    fig = groupedbar(GeneHD',
        bar_position=:stack,
        color= reshape(cmap,1,G),
        xlim=(0.5,lastcol+0.5),
        ylim = (0,N/5),
        xlabel = "HD",
        #ylabel = "",
        xticks = ([2:3:lastcol+0.5;],string.([1/G:3/G:lastcol/G;])),
        #yticks = ([0:0:0;],string.([0:0:0;])),
        xtickfont=font(10),
        ytickfont=font(10),
        guidefont=font(12),
        legendfont = font(10),
        legend = false,
        #label = [""],
        title = "Hamming Distance with Gene Contributions, 10^($(Int64(log10(N)))) Cell Pairs \\n $(source)",
        framestyle=:box,
        top_margin = 3mm,
        size=(600,400))
        annotate!([2:lastcol;],sum(GeneHD[:,2:lastcol],dims=1).+a_dist,Plots.text.(annot, 10, :black, :center)) # annotate the number of genes contributing to each HD value
        savefig(fig,"Henryetal_CellHDwithGenes_$(source).png"); #savefig
end #function

# HDbargene generates a bar graph of the contribution of each gene to cell-cell distance
function HDbargene(SeqVals,N,GeneList,source)
    GeneHD = HD(SeqVals,N,true) # get gene by HD matrix
    GeneNames=readdlm(GeneList) # get list of gene names
    data = sum(GeneHD,dims=2)./sum(GeneHD) # get fraction of cells for each gene
    fig = bar(GeneNames[:],data, #generate barplot
        xrotation = 90,
        bar_width=1,
        xlim = (-.75,length(data)+1),
        ylim = (0,maximum(data)*1.1),
        xlabel = "Genes",
        ylabel = "Fraction of All Cell Pairs",
        xticks = ([0:1:40;],GeneNames),
        #yticks = ([0:0:0;],string.([0:0:0;])),
        xtickfont=font(8),
        ytickfont=font(10),
        guidefont=font(12),
        legendfont = font(10),
        legend = false,
        #label = [""],
        title = "Genes Different Between 10^($(Int64(log10(N)))) Cell Pairs \\n $(source)",
        framestyle=:box,
        grid=false,
        bottom_margin = 12mm,
        top_margin = 3mm,
        size = (600,400))
    savefig(fig,"Henryetal_CellHDGenes_$(source).png"); #savefig
end

# neighbors finds all hamming neighbors in a given data and records the number of neighbors and location of difference for each cell
function neighbors(SeqVals)
    SeqVals .= SeqVals.>0; #Binarize SeqVals
    G,C = size(SeqVals); #Number of genes and cells
    b = 2 .^[G-1:-1:0;]'; #vector of binary place values
    SeqInt = b*SeqVals; #Integer value for each cell
    HN = fill(0,G,C); # Initialize Hamming Neighbor gene X cell matrix
    for i = 1:C #for each cell
        NeighBin = Matrix(I,G,G);
        NeighBin = xor.(SeqVals[:,i],NeighBin);
        NeighInt = b*NeighBin;
        HN[:,i]=sum(NeighInt[:].==SeqInt,dims=2);
    end
    HN
end #function

#histogram of number of hamming neighbors a cell has
function HNhist(SeqVals,source)
    HN = neighbors(SeqVals) # get HN values
    NumHN = round(div(sum(HN),2)/(size(SeqVals,2))^2,digits=3) # fraction of Hamming Neighbor Pairs
    NumHCells = round(sum(sum(HN,dims=1).>=1)/size(SeqVals,2),digits=2) # fraction of cells that have Hamming Neighbors
    annot = "Fraction of HN ≈ $(NumHN)\\nFraction of Cells ≈ $(NumHCells)"; #
    data = log10.(sum(HN,dims=1)[:])
    fig = histogram(data, #generate histogram
        nbins = [0:0.1:maximum(data)+0.1;],
        xlim = (0,maximum(data)+0.1),
        ylim = (0,length(data)/10),
        xlabel = "Number of HN Per Cell, logscale",
        ylabel = "Number of Cells",
        #xticks = ([0:0:0;],string.([0:0:0;])),
        #yticks = ([0:0:0;],string.([0:0:0;])),
        xtickfont=font(10),
        ytickfont=font(10),
        guidefont=font(12),
        legendfont = font(10),
        legend = false,
        title = "Distribution of Hamming Neighbors Across Cells \\n $(source)",
        framestyle=:box,
        top_margin = 3mm,
        size = (600,400))
        annotate!(0.7*maximum(data),0.9*length(data)/10,Plots.text.(annot, 10, :black, :left)) # annotate
    savefig(fig,"Henryetal_CellHN_$(source).png"); #savefig
end

# histogram of contibuting genes
function HNbargene(SeqVals,GeneList,source)
    GeneNames=readdlm(GeneList)
    HN = neighbors(SeqVals) # get HN values
    data = sum(HN,dims=2)./sum(HN)
    fig = bar(GeneNames[:],data, #generate barplot
        xrotation = 90,
        bar_width=1,
        xlim = (-0.05,length(data)+1),
        ylim = (0,maximum(data)*1.1),
        xlabel = "Genes",
        ylabel = "Fraction of all HN",
        xticks = ([0:1:40;],GeneNames),
        #yticks = ([0:0:0;],string.([0:0:0;])),
        xtickfont=font(8),
        ytickfont=font(10),
        guidefont=font(12),
        legendfont = font(10),
        legend = false,
        #label = [""],
        title = "Genes Different Between Hamming Neighbors \\n $(source)",
        framestyle=:box,
        grid=false,
        bottom_margin = 12mm,
        top_margin = 3mm,
        size = (600,400))
    savefig(fig,"Henryetal_CellHNGenes_$(source).png"); #savefig
end

end#module

## it might make sense to expand this so that it takes in the full data set with an option for a subset
