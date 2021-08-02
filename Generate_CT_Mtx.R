## Generate_CT_Mtx.R takes in separated expression data, combines it, subsets it, and creates a phenotype string
## Takes in a MatrixMarket RNA expression matrix, a cell barcode .tsv, and a gene name .tsv and combines them
## (WIP) Takes in subset arguments and subsets matrix  
## Uses a phenos.tsv to record phenotypes for the subset data 


### Command inputs
arg = commandArgs(trailingOnly=TRUE)
#arg[1] = file name stem

## Import expression matrix
library(Matrix)
# load expression matrix
expression_matrix <- readMM(paste0(arg[1],"_matrix.mtx"))

## Add cell names 
colnames(expression_matrix) <- read.delim(paste0(arg[1],"_barcodes.tsv"),header=FALSE)$V1

## Add gene names
rownames(expression_matrix) <- read.delim(paste0(arg[1],"_genes.tsv"),header=FALSE)$V2

## Subset expression matrix (WIP)
#QC_IDX = read.delim("Pd_NonStressedCell_Idx.tsv", header=FALSE) #indices of QC cells
#rowIDX = 1:dim(expression_matrix)[1] #row indices used to subset matrix
#expression_matrix <- expression_matrix[rowIDX,QC_IDX$V1]


### Subset expression matrix for cells of interest
# Phenotypes of Interest
#LOI = "Epi"
#DOI = "D35"
#POI = "BE"

#Lineage
#lineage = read.delim("PdF_Barcode_id.Lineage.aggr.csv", ',', header = FALSE)$V2 
#Donor
#donor = read.delim("PdF_Barcode_id.Patient.aggr.csv",',',header=FALSE)$V2 
#Phenotype
#pheno = read.delim("PdF_Barcode_id.Populations.aggr.csv",',',header =FALSE)$V2
#Subset index
#SS_IDX = lineage==LOI & donor==DOI & pheno==POI

## Generate phenoString
phenoString <- read.delim(paste0(arg[1],"_phenotypes.tsv"),header=TRUE)
phenoString <- phenoString[colnames(expression_matrix),]

## Save files
write.table(phenoString,past0(arg[1],"phenotypes_input"))
writeMM(expression_matrix,paste0(arg[1],"matrix_input",".mtx"))