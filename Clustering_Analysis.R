## Clustering_Analysis.R takes scRNA expression data and assigns expression-based clusters to cells.
## Requires: Seurat, Matrix
## Inputs: _matrix.mtx, _barcodes.tsv, _genes.tsv, _phenotypes.tsv

### load libraries
library(Seurat)
library(Matrix)


### Input parameters
### Command inputs
arg = commandArgs(trailingOnly=TRUE)
#arg[1] = filename stem



### load RNAseq data
# load expression matrix
expression_matrix <- readMM(paste0(arg[1],"_matrix.mtx"))
# convert from dgTMatrix to dgCMatrix
expression_matrix <- as(expression_matrix,"dgCMatrix")
## Add cell names 
colnames(expression_matrix) <- read.delim(paste0(arg[1],"_barcodes.tsv"),header=FALSE)$V1
## Add gene names
row.names(expression_matrix) <- read.delim(paste0(arg[1],"_genes.tsv"),header=FALSE)$V2

### Seurat Analysis
#create Seurat object to store analysis
expression_obj <- CreateSeuratObject(counts = expression_matrix, project = "proj1", min.cells = 0, min.features = 0)
# Normalize data for dimensionality analysis
expression_obj <- NormalizeData(expression_obj, normalization.method = "LogNormalize", scale.factor = 10000)
# Find features for dimensionality analysis
expression_obj <- FindVariableFeatures(expression_obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
#top10 <- head(VariableFeatures(expression_obj),10)

# Scale data based on gene expression
all.genes <- rownames(expression_obj)
expression_obj <- ScaleData(expression_obj, features = all.genes)

# Run PCA analysis
expression_obj <- RunPCA(expression_obj, features = VariableFeatures(object = expression_obj))
# Analyze pca to look for relevant number of dimensions
png(file=paste0(arg[1],"_elbow_plot.png"))
ElbowPlot(expression_obj)
dev.off()


# Cluster parameters
ndims<- c(10,15,20)
rez <- c(0.5,0.75,1)

for(i in 1:length(ndims)) {

  for(j in 1:length(rez)) {
    
    #Cluster Cells
    expression_obj <- FindNeighbors(expression_obj, dims = 1:ndims[i])
    expression_obj <- FindClusters(expression_obj, resolution = rez[j])
    cluster_ids <- Idents(expression_obj)
    
    
    ## Visualize Cluster data
    # Note: Use the same viz(UMAP/tSNE) as CytoTRACE
    
    # Visualize cell clusters with UMAP
    #expression_obj <- RunUMAP(expression_obj, dims = 1:ndims[i])
    #DimPlot(expression_obj, reduction = "umap")
    
    # Visualize cell clusters with tsne
    expression_obj <- RunTSNE(expression_obj, dims = 1:ndims[i])
    png(file=paste0(arg[1],"_tSNE_d",ndims[i],"r",rez[j],".png"))
    DimPlot(expression_obj, reduction = "tsne")
    dev.off()
    
    
    ## Compare CytoTRACE results with cluster results
    #load cytoTRACE data
    cyto_data <- read.delim(paste0(arg[1],"_CytoTRACE_plot_table.txt"), header=TRUE)
    
    # add cluster data
    cyto_data$Cluster <-cluster_ids
    
    #plot Cluster ID versus Phenotype
    #par(mfcol=c(5,1)) 
    #hist(as.numeric(cluster_ids[cyto_data$Phenotype=="BE"]))
    #hist(as.numeric(cluster_ids[cyto_data$Phenotype=="Club"]))
    #hist(as.numeric(cluster_ids[cyto_data$Phenotype=="Hillock"]))
    #hist(as.numeric(cluster_ids[cyto_data$Phenotype=="LE"]))
    #hist(as.numeric(cluster_ids[cyto_data$Phenotype=="NE"]))
    
    #plot cytoTRACE verus cluster ID
    png(file=paste0(arg[1],"_boxplot_d",ndims[i],"r",rez[j],".png"))
    boxplot(CytoTRACE~Cluster,data=cyto_data)
    dev.off()

    
    ## Save cluster data
    write.csv(cluster_ids,paste0(arg[1],"_ClusterIDs_d",ndims[i],"r",rez[j],".csv"))

  }
}

















