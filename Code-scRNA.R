
#install.packages("Seurat")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("GSVA")
#BiocManager::install("GSEABase")
#BiocManager::install("limma")
#BiocManager::install("SingleR")
#BiocManager::install("celldex")
#BiocManager::install("qlcMatrix’")



library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(SingleR)
library(monocle)



logFCfilter=1               
adjPvalFilter=0.05          
#inputFile="27N.csv"       

setwd("C:\\Users\\Lenovo\\Desktop\\gastric cancer")    

data<-Read10X(data.dir = "C:\\Users\\Lenovo\\Desktop\\---/")


#质控
pbmc=CreateSeuratObject(counts = data,project = "seurat", min.cells=3, min.features=200, names.delim = "_")


pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^mt-")


pdf(file="01.featureViolin.pdf", width=10, height=6)
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pbmc=subset(x = pbmc, subset = nFeature_RNA >500 & nFeature_RNA <6000 & percent.mt < 20)   


#可视化
pdf(file="01.featureCor.pdf",width=10,height=6)
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",,pt.size=1.5)
CombinePlots(plots = list(plot1, plot2))
dev.off()

#数据标准化
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
pbmc <- CellCycleScoring(pbmc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
pbmc<-ScaleData(pbmc,vars.to.regress = c("percent.mt","S.Score","G2M.Score"),verbose = FALSE)

#寻找高变基因
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 3000)

#高变基因可视化
top10 <- head(x = VariableFeatures(object = pbmc), 10)
pdf(file="01.featureVar.pdf",width=10,height=6)
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()


###################################02.PCA主成分分析###################################
##PCA
pbmc=ScaleData(pbmc)          #数据标准化

library(harmony)

###perform linear dimensional reduction###


pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))   

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
#dev.off()
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
pbmc<-RunHarmony(pbmc,group.by.vars = c("orig.ident"), plot_convergence = TRUE) #去除批次效应


harmony_embeddings <- Embeddings(pbmc, 'harmony')
dim(harmony_embeddings)
p3 <- DimPlot(object = pbmc, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p4 <- VlnPlot(object = pbmc, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p3,p4))


###未进行去除批次效应的图
p1 <- DimPlot(object = pbmc, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = pbmc, features = "PC_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p1,p2))

#对比去掉批次
p2+ p4

##其实差距不大（去掉批次效应）

ElbowPlot(pbmc,ndims = 50)###查看最合适的PC

###################################03.TSNE###################################
##TSNE
pcSelect=50    ###需要选择最合适的PC
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)       
pbmc <- FindClusters(object = pbmc, resolution =0.1)         
pbmc <- RunUMAP(object = pbmc, dims = 1:pcSelect)  

#可视化           
pdf(file="03.TSNE.pdf",width=7,height=7)
UMAPPlot(object = pbmc, pt.size = 0.2, label = F)   
dev.off()



##查找每个聚类的差异基因
pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.markers,file="03.clusterMarkers.txt",sep="\t",row.names=F,quote=F)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#绘制marker在各个cluster的热图
pdf(file="03.tsneHeatmap.pdf",width=12,height=9)
DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()
dev.off()






VlnPlot(pbmc, features = c("Map2k4"))
FeaturePlot(pbmc, features = c("MS4A1", "CD79A"))
RidgePlot(pbmc, features = c("CD163", "C1QA", "C1QB", "C1QC", "TREM2"), ncol = 1)
DotPlot(pbmc, features = c("Hnf1a", "Apoa1", "Pecam1", "Cd163", "Csf1r","Dcn", "Ebf1","Cd19", "Cd3d", "Cd3e", "Epcam","Hnf1b"),,cols=c("gray","red")) + RotatedAxis()
DoHeatmap(subset(pbmc, downsample = 100), features = c("CD163", "C1QA", "C1QB", "C1QC", "TREM2"), size = 3)




######biomarkers###

FeaturePlot(pbmc,features = c("ALB","TF","APOB","CYP3A4","CYP2E1","ASGR1","PCK1","ASS1","APOE","AFP","CYP2C8"),reduction = "tsne",label = T)
FeaturePlot(pbmc,features = c("EPCAM","KRT19"),reduction = "tsne",label = T)

##免疫
FeaturePlot(pbmc,features = c("PTPRC"),reduction = "tsne",label = T,cols = c("grey","red"))

##TNK cell
FeaturePlot(pbmc,features = c("CD3D","CD3E"),reduction = "tsne",label = T,cols = c("grey","red"))

##treg
FeaturePlot(pbmc,features = c("CD4","IL2RA","FOXP3"),reduction = "tsne",label = T,cols = c("grey","red"))

##nk
FeaturePlot(pbmc,features = c("KLRF1","FGFBP2","KLRD1"),reduction = "tsne",label = T,cols = c("grey","red"))

##endothelial cell
FeaturePlot(pbmc,features = c("PECAM1","VWF"),reduction = "tsne",label = T,cols = c("grey","red"))


####fibroblasts
FeaturePlot(pbmc,features = c("ACTA2","COL1A1","DCN"),reduction = "tsne",label = T,cols = c("grey","red"))

###B cell
FeaturePlot(pbmc,features = c("CD19","MS4A1","CD79A"),reduction = "tsne",label = T,cols = c("grey","red"))


FeaturePlot(pbmc,features = c("IGHA1","MZB1"),reduction = "tsne",label = T,cols = c("grey","red"))


###myeloid cell
FeaturePlot(pbmc_mac,features = c("CD14","CD68"),reduction = "umap",label = T,cols = c("grey","red"))


###mast cell
FeaturePlot(pbmc,features = c("MS4A2"，"MS4A1"),reduction = "umap",label = T,cols = c("grey","red"),pt.size=0.8)

###tumor cell
FeaturePlot(pbmc,features = c("EPCAM"),reduction = "tsne",label = T,cols = c("grey","red"))


###pDC cell
FeaturePlot(pbmc,features = c("IL3RA","CLEC4C"),reduction = "tsne",label = T,cols = c("grey","red"))


pbmc<-RenameIdents(pbmc,"0"="POSTN fibroblast","1"="POSTN fibroblast","2"="LOX+ fibroblast","3"="POSTN fibroblast","4"="HP fibroblast","5"="unspecific fibroblast","6"="SFRP2 fibroblast",
                              "7"="CCL5 iCAF","8"="CCL11 iCAF","9"="TFF1 fibroblast","10"="CCL3 iCAF","11"="MZB1 iCAF","12"="12","13"="13","14"="14","15"="15","16"="16","17"="17","18"="18",
                              "19"="19","20"="20")



VlnPlot(pbmc,features = c("IL6"),reduction = "tsne",label = T)


###提取各个细胞群

c1<- Idents(pbmc)%in%c("LOX+ fibroblast")
c2<-pbmc[,c1]
table(c1)
saveRDS(c2,file = "LOX+fibroblast.Rds")





###################################05.monocle R模拟时序分析###################################
monocle.matrix=as.matrix(pbmc@assays$RNA@data)
monocle.sample=pbmc@meta.data
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
monocle.clusterAnn=clusterAnn
monocle.markers=sig.markers

data <- as(as.matrix(monocle.matrix), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = monocle.sample)
fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd)
names(pData(cds))[names(pData(cds))=="seurat_clusters"]="Cluster"
pData(cds)[,"Cluster"]=paste0("cluster",pData(cds)[,"Cluster"])

clusterAnn=as.character(monocle.clusterAnn[,2])
names(clusterAnn)=paste0("cluster",monocle.clusterAnn[,1])
pData(cds)$cell_type2 <- plyr::revalue(as.character(pData(cds)$Cluster),clusterAnn)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- setOrderingFilter(cds, as.vector(sig.markers$gene))
#plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2, reduction_method = 'DDRTree')
cds <- orderCells(cds)


pdf(file="05.trajectory.State.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "State")
dev.off()

pdf(file="05.trajectory.Pseudotime.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "Pseudotime")
dev.off()

pdf(file="05.trajectory.cellType.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "cell_type2")
dev.off()

pdf(file="05.trajectory.cluster.pdf",width=6.5,height=6)
plot_cell_trajectory(cds, cols= cell_type_cols,color_by = "Cluster")
dev.off()


groups=subset(pData(cds),select='State')
pbmc=AddMetaData(object=pbmc, metadata=groups, col.name="group")
geneList=list()
for(i in levels(factor(groups$State))){
	pbmc.markers=FindMarkers(pbmc, ident.1 = i, group.by = 'group')
	sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
	sig.markers=cbind(Gene=row.names(sig.markers), sig.markers)
	write.table(sig.markers,file=paste0("05.monocleDiff.", i, ".txt"),sep="\t",row.names=F,quote=F)
	geneList[[i]]=row.names(sig.markers)
}

unionGenes=Reduce(union,geneList)
write.table(file="05.monocleDiff.union.txt",unionGenes,sep="\t",quote=F,col.names=F,row.names=F)



#可视化
library(RColorBrewer)
cell_type_cols <- c(brewer.pal(9, "Set1"), "#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F","#FF6A6A","#7FFFD4", "#AB82FF","#90EE90","#00CD00","#008B8B","#6495ED","#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00","#000000","#708090")  
DimPlot(pbmc, label = T, cols= cell_type_cols, pt.size = 0.3, repel = T)+ NoLegend()


cell_type_cols = c("#FF6A6A","#FF34B3", "#AB82FF","#CD5C5C")

DimPlot(pbmc,group.by = "orig.ident",label = T,cols= cell_type_cols, pt.size = 0.3, repel = T)


library(RColorBrewer)
cell_type_cols <- c(brewer.pal(9, "Set1"), "#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F","#FF6A6A","#7FFFD4", "#AB82FF","#90EE90","#00CD00","#008B8B","#6495ED","#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00","#000000","#708090")  
 library(MySeuratWrappers)  
 markers <- c("SPP1")  
 VlnPlot(pbmc, features = markers,stacked=T,pt.size=0,cols = cell_type_cols, x.lab = " ", y.lab = " ")+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())#不显示坐标刻度 
 

 FeaturePlot(pbmc,features = c("MMP19"),reduction = "tsne",label = T)
 
 
 
 
#相关性
 ggscatterstats(dat, 
               y =RNASE1, 
               x =CD163,
               type = "pearson",
               centrality.para = "mean",                              
               margins = "both",                                         
               xfill = "#009E73",
               yfill = "#D55E00", 
               marginal.type = "density",    #类型可以换成density,boxplot,violin,densigram
               title = "Relationship between LDHA and PKM")
               
               
               
               
                        
               
#需要展示的基因，可以修改
showGenes=c("SPP1") 

#绘制marker在各个cluster的散点图
pdf(file="03.markerScatter.pdf",width=10,height=12)
FeaturePlot(object = pbmc, features = showGenes, cols = c("blue", "red"))
dev.off()



#AddModuleScore计算
library(tidyverse)
library(Matrix)
library(cowplot)

##计算评分AddModuleScore：
gene <- read.table("C:\\Users\\chen\\Desktop\\单细胞亚群分析\\联合分析/gene.txt", quote="\"", comment.char="")
gene = gene$V1

gene = list(gene)

pbmc = AddModuleScore(object = pbmc , features =gene,ctrl = 100 ,name = 'CD_Features')

colnames(pbmc@meta.data)[8] = 'Score'
VlnPlot(pbmc,features = 'Score')

pbmc$cell_type = pbmc$orig.ident
data = FetchData(pbmc,vars = c("cell_type","Score"))
p = ggplot(data,aes(cell_type,Score))
p+geom_boxplot()+theme_bw()+RotatedAxis()


write.table(data,file="data.txt",sep="\t",row.names=F,quote=F)

#可视化
library(ggplot2)

mydata<- FetchData(pbmc,vars = c("umap_1","umap_2","Score"))
a <- ggplot(mydata,aes(x = umap_1,y =umap_2,colour = Score))+
  geom_point(size = 0.3)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colours = c('#333366',"#6666FF",'#FFCC33','red'))
                                             
                                             
                                             
a+ theme_bw() + theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
                      
                      
                      
                      

#组间差异分析
logFCfilter=1
adjPvalFilter=0.05
groups=Idents(pbmc)
names(groups)=colnames(pbmc)


pbmc=AddMetaData(object=pbmc, metadata=groups, col.name="group")
pbmc.markers=FindMarkers(pbmc, ident.1 = "high", ident.2 = "low", group.by = 'group')
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
sig.markers=cbind(Gene=row.names(sig.markers), sig.markers)
write.table(sig.markers,file="diffGene.txt",sep="\t",row.names=F,quote=F)

##NicheNet--cell-cell communications


library(nichenetr)
library(Seurat)
library(tidyverse)



ligand_target_matrix <- readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
# weighted_networks列表包含两个数据框，lr_sig是配体-受体权重信号网络，gr是配体-靶基因权重调控网络



##读入单细胞数据
scRNA <- pbmc
scRNA <- UpdateSeuratObject(scRNA)


nichenet_output = nichenet_seuratobj_aggregate(seurat_obj = scRNA, 
                                               top_n_ligands = 20,
                                               receiver = "Monocyte", 
                                               sender = "LOX+ fibroblast",
                                               condition_colname = "tissue_type", 
                                               condition_oi = "T", 
                                               condition_reference = "N", 
                                               ligand_target_matrix = ligand_target_matrix, 
                                               lr_network = lr_network, 
                                               weighted_networks = weighted_networks, 
                                               organism = "human")
 
## 查看配体活性分析结果
# bona_fide_ligand=False代表PPI预测未经实验证实的配体-受体。
x <- nichenet_output$ligand_activities
write.csv(x, "ligand_activities.csv", row.names = F)
# 查看top20 ligands
nichenet_output$top_ligands
# 查看top20 ligands在各个细胞亚群中表达情况
p = DotPlot(scRNA, features = nichenet_output$top_ligands, cols = "RdYlBu") + RotatedAxis()
ggsave("top20_ligands.png", p, width = 12, height = 6)
# 按"aggregate"的分类对比配体的表达情况
p = DotPlot(scRNA, features = nichenet_output$top_ligands, split.by = "aggregate") + RotatedAxis()
ggsave("top20_ligands_compare.png", p, width = 12, height = 8)
# 用小提琴图对比配体的表达情况
p = VlnPlot(scRNA, features = nichenet_output$top_ligands, 
            split.by = "aggregate",  pt.size = 0, combine = T)
ggsave("VlnPlot_ligands_compare.png", p, width = 12, height = 8)




## 查看配体调控靶基因
p = nichenet_output$ligand_target_heatmap
ggsave("Heatmap_ligand-target.png", p, width = 12, height = 6)
# 更改热图的风格
p = nichenet_output$ligand_target_heatmap + 
     scale_fill_gradient2(low = "whitesmoke",  high = "royalblue", breaks = c(0,0.0045,0.009)) + 
     xlab("anti-LCMV response genes in CD8 T cells") + 
     ylab("Prioritized immmune cell ligands")
ggsave("Heatmap_ligand-target2.png", p, width = 12, height = 6)
# 查看top配体调控的靶基因及其评分
x <- nichenet_output$ligand_target_matrix
#x2 <- nichenet_output$ligand_target_df
write.csv(x, "ligand_target.csv", row.names = F)




## 查看受体情况
# 查看配体-受体互作
p = nichenet_output$ligand_receptor_heatmap
ggsave("Heatmap_ligand-receptor.png", p, width = 12, height = 6)
x <- nichenet_output$ligand_receptor_matrix
#x <- nichenet_output$ligand_receptor_df
write.csv(x, "ligand_receptor.csv", row.names = F)
# 查看受体表达情况
p = DotPlot(scRNA %>% subset(idents = "CD8 T"), 
             features = nichenet_output$top_receptors, 
             split.by = "aggregate") + RotatedAxis()
ggsave("Receptors_Expression_dotplot.png", p, width = 12, height = 6)
p = VlnPlot(scRNA %>% subset(idents = "CD8 T"), features = nichenet_output$top_receptors, 
             split.by = "aggregate", pt.size = 0, combine = T, ncol = 8)
ggsave("Receptors_Expression_vlnplot.png", p, width = 12, height = 8)
# 有文献报道的配体-受体
# Show ‘bona fide’ ligand-receptor links that are described in the literature and not predicted based on PPI
p = nichenet_output$ligand_receptor_heatmap_bonafide
ggsave("Heatmap_ligand-receptor_bonafide.png", p, width = 8, height = 4)
x <- nichenet_output$ligand_receptor_matrix_bonafide
#x <- nichenet_output$ligand_receptor_df_bonafide
write.csv(x, "ligand_receptor_bonafide.csv", row.names = F)