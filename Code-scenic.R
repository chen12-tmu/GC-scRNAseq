rm(list=ls())
library(Seurat)
library(dplyr)
library(future)
library(future.apply)
plan("multiprocess", workers = 6) ###set the compute core
options(future.globals.maxSize = 60000 * 1024^2)
getwd()



#rm(list=ls())

## Required
BiocManager::install(c("AUCell", "RcisTarget"))
BiocManager::install(c("GENIE3")) # Optional. Can be replaced by GRNBoost

## Optional (but highly recommended):
# To score the network on cells (i.e. run AUCell):
BiocManager::install(c("zoo", "mixtools", "rbokeh"))
# For various visualizations and perform t-SNEs:
BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"))
# To support paralell execution (not available in Windows):
BiocManager::install(c("doMC", "doRNG"))
# To export/visualize in http://scope.aertslab.org
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)

devtools::install_github("aertslab/SCENIC") 



library(RcisTarget)
library(SCENIC)
library(doMC)
library(doRNG)
library(AUCell)
library(GENIE3)
library(SCENIC)
library(dplyr)
library(Seurat)
library(future)
library(future.apply)
plan("multiprocess", workers = 6) 
###step_1 setup the options###
library(SCENIC)

setwd("~/Documents_PC/scRNA-seq/Data/SCENIC")

org="hgnc"
dbDir="./cisTarget_databases"
myDatasetTitle="SCENIC_myeloid"
#data(defaultDbNames)
dbs <- defaultDbNames[["hgnc"]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=6)
getwd()


####step_1 load the data###
dc<-readRDS("int/Myeloid.Rds")

exprMat<-as.matrix(dc@assays$RNA@counts)
dim(exprMat)
####step_3 filter the genes#####
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept,]
dim(exprMat_filtered)
exprMat_filtered <- log2(exprMat_filtered+1)

###step_4 run correlation
runCorrelation(exprMat_filtered, scenicOptions)
runGenie3(exprMat_filtered, scenicOptions)
##we could save the data to release the ###
save(exprMat_filtered,scenicOptions,file = "input_GENIE3_dc.Rdata")



load("./input_GENIE3_dc.Rdata")
#load("input_GENIE3_endo.Rdata")
scenicOptions@settings$nCores<-10
#scenicOptions@settings$nCores
#scenicOptions@settings$dbs
#setwd("../SCENIC/")

runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions) 

save.image(file="Myeloid_SCENIC.Rdata")

scenicOptions@settings$nCores<-1 ###change into 1 core or wise there will be error
runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)


load("./Myeloid_SCENIC.Rdata")
####next we analysis  
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
#table(tam$seurat)
library(Seurat)
#dc<-readRDS("./int/dc.Rds")

dim(regulonAUC@assays@data$AUC)
table(rownames(dc@meta.data)==colnames(regulonAUC@assays@data$AUC))


data<-getAUC(regulonAUC)

dc$CellType


group<-factor(dc$CellType,levels=c("Monocytes", "TAM", "CD1C_DC", "CLEC9A_DC"))
table(group)
FeaturePlot(dc,reduction = "tsne",features=c("TOP2A","S100A8","CD1C","AIF1","APOE"))
levels(group)
table(names(group)==colnames(data))

data


library(singleseqgset)

logfc.data <- logFC(cluster.ids=group,expr.mat=as.matrix(data))
names(logfc.data)

logfc<-do.call(cbind,logfc.data[[2]])
logfc<-as.data.frame(logfc)

gene1<-row.names(logfc[order(logfc$Monocytes,decreasing = T),])[1:5]
gene2<-row.names(logfc[order(logfc$TAM,decreasing = T),])[1:5]
gene3<-row.names(logfc[order(logfc$CD1C_DC,decreasing = T),])[1:5]
gene4<-row.names(logfc[order(logfc$CLEC9A_DC,decreasing = T),])[1:5]

genes<-unique(c(gene1,gene2,gene3,gene4))
genes
#gene3


library(ComplexHeatmap)
names(group)
dc<-dc[,names(group)]

Idents(dc)<-group
dc$seurat_final<-Idents(dc)

cluster_info<-sort(dc$seurat_final)
cluster_info[1:5]

library(circlize)
library(RColorBrewer)
mat <- as.matrix(data[genes, names(cluster_info)])

dc@meta.data[1:5,]
table(dc$tissue_type)
annotation_col<-data.frame(Cell_group=Idents(dc)[names(cluster_info)],
                           Site=dc$tissue_type[names(cluster_info)])

cor<-readRDS("./int/colVars.Rds")
cor
ann_colors = list(
  Cell_group = c("Monocytes"="#FED439FF", "TAM"="#709AE1FF", "CD1C_DC"="#8A9197FF","CLEC9A_DC"="#D2AF81FF"),
  #Patient=c("SC120"="#FF6F00FF","SC153"="#C71000FF","SC235"="#008EA0FF","SC243"="#8A4198FF","SC265"="#5A9599FF","SC266"="#FF6348FF"),
  Site=c("Normal"="navyblue","Tumor"="firebrick")
)
?ComplexHeatmap::pheatmap
ComplexHeatmap::pheatmap(mat,color =colorRamp2(c(-2.5, 0, 2.5), c("#003366", "white", "#993333")),annotation_col = annotation_col,annotation_colors = ann_colors,cluster_rows = F,cluster_cols = F,border_color = NA,
                         show_colnames = F,fontsize = 8,scale= "row",use_raster=F)


save.image(file="Myeloid_50_cell.Rdata")



