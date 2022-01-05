## Pre-processing of single cell RNAseq data from Henry et al(2018)
## Pre-processing steps include:
##    Aggregating data from all single cell assays (PD, D17, D27)
##    Attaching metadata (Phenotype, Donor,Zone,Barcoding Well, Multiplet Status, Percent Stress Signature)
##    Filter cells with <200 unique genes
##    Filter genes expressed in <3 cells
##    Annotate cells with cell cycle prediction
##    Remove cells with <500 unique genes
##    Remove cells with <10% mitochondrial gene expression
##    Annotate cells with percent stress signature gene expression
##    Remove mitochondrial genes
##    Annotate Doublets with DoubletFinder
##    Remove cells with <3000 unique genes expressed



gc()

### load libraries
library(Matrix)
library(Seurat)
library(SeuratDisk)
library(DoubletFinder)
library(ggplot2)

### Define data sources
data_source = c("GSE117403_PD","GSE117403_D17","GSE117403_D27")


### load RNAseq data
Expression_matrix <- Read10X(data.dir=paste0("./Data/10x/",data_source[1],"/"))
for (i in 2:length(data_source)) {
  Expression_matrix <- cbind(Expression_matrix, Read10X(data.dir= paste0("./Data/10x/",data_source[i],"/")) )
}


### Create Seurat Object
Exp_mtx <- CreateSeuratObject(counts = Expression_matrix)


### Attach metadata
meta_data <- read.delim(paste0("./Data/10x/",data_source[1],"/metadata.tsv"),header=TRUE)
for (i in 2:length(data_source)) {
  meta_data <- rbind(meta_data, read.delim(paste0("./Data/10x/",data_source[i],"/metadata.tsv"),header=TRUE) )
}
#Phenotype
Exp_mtx[["Phenotype"]] <- meta_data$Phenotype
#Donor
Exp_mtx[["Donor"]] <- meta_data$Donor
#Zone
Exp_mtx[["Zone"]] <- meta_data$Zone
#10x Well
Exp_mtx[["SeqWell"]] <- meta_data$SeqWell
#Multiplet Status
Exp_mtx[["Multiplet"]] <- "N/A"
## Stress status
Exp_mtx[["percent.stress"]] <- "N/A"


### Filter low quality cells
# Remove cells with fewer than 200 unique genes expressed
Exp_mtx <- subset(x = Exp_mtx, subset = nFeature_RNA >= 200)
# Remove genes expressed in 3 or fewer cells
Exp_mtx <- Exp_mtx[rowSums(Exp_mtx@assays$RNA@counts>0) >= 3,]


### Predict cell cycle stage
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
# NormalizeData and Scale for annotation
Exp_mtx <- NormalizeData(Exp_mtx)
Exp_mtx <- FindVariableFeatures(Exp_mtx, selection.method = "vst")
Exp_mtx <- ScaleData(Exp_mtx, features = rownames(Exp_mtx))
#Assign cell-cycle score
Exp_mtx <- CellCycleScoring(Exp_mtx, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


### Remove cells with < 500
Exp_mtx <- subset(x=Exp_mtx,subset= nFeature_RNA > 500)
### Remove cells with >10% of their transcriptome as mitochondrial genes
Exp_mtx[["percent.mt"]] <- PercentageFeatureSet(Exp_mtx, pattern = "^MT-")
Exp_mtx <- subset(x=Exp_mtx,subset= percent.mt < 10)


### Annotate stressed cells
Stress_Genes <- read.csv("./Analysis/PreProcessing/StressGenes.csv",header=FALSE)
Exp_mtx$percent.stress <- PercentageFeatureSet(Exp_mtx, features = Stress_Genes$V1)

  
### Remove Mitochondrial Genes
mito.genes <- grep(pattern="^MT-",x=rownames(Exp_mtx),value=TRUE)
Exp_mtx <- Exp_mtx[!rownames(Exp_mtx) %in% mito.genes,]


### Do not perform this step unless there is clear evidence that it should be regressed out.
### Scale UMI counts and remove variation due to differences in UMI/cell and percent mitochondrial genes. 
### Leave cell cycle phase related variation
#Exp_mtx <- ScaleData(object=Exp_mtx,vars.to.regress=c("nCount_RNA","percent.mt"),display.progress=FALSE,do.par=TRUE,num.cores=45)


### Annotate potential doublets using DoubletFinder
Doublet_batch <- unique(Exp_mtx@meta.data$SeqWell) # label of unique single cell labeling batches
for(i in 1:length(Doublet_batch)){
  ## Isolate single 10x sequencing well
  Doublet_mtx <- Exp_mtx[,Exp_mtx$SeqWell==Doublet_batch[i]]
  ## Dimensional Analysis
  # LogNormalize to ln(10k transcripts/cell +1)
  Doublet_mtx <- NormalizeData(object=Doublet_mtx,normalization.method="LogNormalize",scale.factor=10000)
  # Find Variable Features
  Doublet_mtx <- FindVariableFeatures(Doublet_mtx, selection.method = "vst")
  # Scale data based on gene expression
  all.genes <- rownames(Doublet_mtx)
  Doublet_mtx <- ScaleData(Doublet_mtx, features = all.genes)
  # Run PCA analysis
  Doublet_mtx <- RunPCA(Doublet_mtx, features = VariableFeatures(object = Doublet_mtx))
  # Analyze pca to look for relevant number of dimensions
  fig <- ElbowPlot(Doublet_mtx)
  ggsave(file=paste0("./Analysis/PreProcessing/DoubletBatch",Doublet_batch[i],"_elbow_plot.png"),fig)
  # Cluster parameters
  ndims<- 10
  rez <- 0.5
  #Cluster Cells
  Doublet_mtx <- FindNeighbors(Doublet_mtx, dims = 1:ndims)
  Doublet_mtx <- FindClusters(Doublet_mtx, resolution = rez)
  ## Visualize Cluster data
  # Visualize cell clusters with UMAP/tSNE
  Doublet_mtx <- RunUMAP(Doublet_mtx, dims = 1:ndims)
  # Plot UMAP of Phenotypes
  fig <- DimPlot(Doublet_mtx, reduction = "umap",group.by=c("Phenotype","seurat_clusters"),label = TRUE)
  ggsave(file=paste0("./Analysis/PreProcessing/DoubletBatch",Doublet_batch[i],"_UmapPhenotypes.png"),fig)
  ## Estimate pK
  sweep.res.list_mtx <- paramSweep_v3(Doublet_mtx, PCs = 1:10, sct = FALSE)
  sweep.stats_mtx <- summarizeSweep(sweep.res.list_mtx, GT = FALSE)
  bcmvn_mtx <- find.pK(sweep.stats_mtx)
  fig<- ggplot(data=bcmvn_mtx,aes(x=pK,y=BCmetric))
  fig <- fig+geom_point() + scale_x_discrete(breaks = c(0, 0.1, 0.2, 0.3))
  ggsave(file=paste0("./Analysis/PreProcessing/DoubletBatch",Doublet_batch[i],"_bcmvn_plot.png"),fig)
  pK <- as.double(as.character(bcmvn_mtx$pK[which.max(bcmvn_mtx$BCmetric)]))
  # Estimate homotypic Doublet Proportion
  annotations <- Doublet_mtx@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(0.075*nrow(Doublet_mtx@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  # Identify Doublets
  Doublet_mtx <- doubletFinder_v3(Doublet_mtx, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  Doublet_mtx <- doubletFinder_v3(Doublet_mtx, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = paste0("pANN_0.25_",pK,"_",nExp_poi), sct = FALSE)
  # Record status
  Exp_mtx@meta.data$Multiplet[Exp_mtx$SeqWell==Doublet_batch[i]] <- Doublet_mtx@meta.data[[paste0("DF.classifications_0.25_",pK,"_",nExp_poi.adj)]]

  # Plot UMAP with Doublets
  fig <- DimPlot(Doublet_mtx, reduction = "umap",group.by=paste0("DF.classifications_0.25_",pK,"_",nExp_poi.adj),pt.size=0.5)
  ggsave(file=paste0("./Analysis/PreProcessing/DoubletBatch",Doublet_batch[i],"_UMAPDoublet_plot.png"),fig)
}


### Remove cells with > 3000 genes expressed
Exp_mtx <- subset(x=Exp_mtx,subset= nFeature_RNA < 3000)

SaveH5Seurat(Exp_mtx,"Henryetal_PreprocessedData",overwrite=TRUE)













