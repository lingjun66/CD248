#title: single cell analysis-GSE145173
#date: 2024.2.16

'Human and Mouse PDGFRab+ Kidney Fibrosis 10X genomics'

#Load software
library(Seurat)
library(SingleR)
library(celldex)
library(tidyverse)
library(cowplot)
library(patchwork)
library(R.utils)
library(data.table)
library(Matrix)
library(pacman)
library(future) 
options(future.globals.maxSize = 10* 1024^3)
rm(list = ls())

setwd("~/GSE145173/Result")
dir.create('RDS')


#Set data path
dir='～/GSE145173/Data/Human_PDGFRb/' 


#Load data
expression_matrix <- ReadMtx(mtx = paste(dir,"pdgfrbMap_UMI_counts.mtx", sep="/"),
                             cells = paste(dir,"pdgfrbMap_UMI_counts_colData.txt", sep="/"),
                             features = paste(dir,"pdgfrbMap_UMI_counts_rowData.txt", sep="/"),
                             skip.cell = 1,
                             skip.feature = 1,
                             feature.sep = ",",
                             feature.column = 1)

#Create seurat
scRNA <- CreateSeuratObject(counts = expression_matrix,
                                 min.features = 200,
                                 min.cells = 3)


#Data standardization and normalization
set.seed(12345)

scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 1e4) %>% 
         FindVariableFeatures(selection.method = "vst", Features = 2000, verbose = FALSE) %>% 
         ScaleData()


scRNA <- RunPCA(scRNA, npcs = 50,
                features = VariableFeatures(object = scRNA), 
                )

plot1 <- DimPlot(scRNA, shuffle =T,
                 raster=FALSE,
                 reduction = "pca", 
                 group.by="orig.ident") 
plot2 <- ElbowPlot(scRNA, 
                   ndims=50, 
                   reduction="pca") 
plotc <- plot1+plot2
ggsave("QC/pca.pdf", plot = plotc, width = 10, height = 4) 
ggsave("QC/pca.png", plot = plotc, width = 10, height = 4) 


#Pick PCA number
pc.num=1:40

#Reduction
dir.create("Cluster") 

scRNA <- RunUMAP(scRNA, 
                 dims = pc.num, 
                 reduction = "pca")
scRNA <- RunTSNE(scRNA, 
                  dims = pc.num, 
                  reduction = "pca")
#Clustering
scRNA <- FindNeighbors(scRNA, 
                       dims = pc.num, 
                       reduction = "pca")

scRNA <- FindClusters(scRNA, graph.name = "RNA_snn",
                      algorithm = 1,
                      resolution =c(seq(0,1,.1)), 
                      verbose = FALSE)
#Check result
library(clustree)
p1 <- clustree(scRNA, prefix = "RNA_snn_res.")
ggsave("Cluster/clustree.pdf", p1, width = 12, height = 16)
ggsave("Cluster/clustree.png", p1, width = 12, height = 16)

p2 <- DimPlot(scRNA,
              reduction = 'umap', 
              group.by = "RNA_snn_res.0.1",  
              repel=T, label =T)
ggsave("Cluster/UMAP_cluster_res.0.1.pdf", p2, width = 8, height = 8)
ggsave("Cluster/UMAP_cluster_res.0.1.png", p2, width = 8, height = 8)


p3 <- DimPlot(scRNA,
              reduction = 'tsne', 
              group.by = "RNA_snn_res.0.1",  
              repel=T, label =T)
ggsave("Cluster/tSNE_cluster_res.0.1.pdf", p3, width = 8, height = 8)
ggsave("Cluster/tSNE_cluster_res.0.1.png", p3, width = 8, height = 8)



#Marker gene mapping
DefaultAssay(scRNA) <- "RNA"


markerGenes <- c("COX4I2",        #locate Pericyte (Pe)
                 "DCN",           #locate Fibroblast (Fib)
                 "POSTN",         #locate Myofibroblast (MF)
                 "RERGL",         #locate Vascular smooth muscle cells (VSMCs)
                 "GATA3",         #locate Mesangial cells (Mesa)
                 "PECAM1",        #locate Endothelial cells (Endo)
                 "EPCAM",         #locate Epithelial cells (Epi)
                 "PTPRC"          #locate Immune cells
)

p4.1 <- FeaturePlot(scRNA,
                  features = markerGenes, 
                  reduction ="umap", ncol=4)
ggsave("Cluster/UMAP_markergenes_res.0.1.pdf", p4.1, bg="#ffffff", width=18, height=8)
ggsave("Cluster/UMAP_markergenes_res.0.1.png", p4.1, bg="#ffffff", width=18, height=8)


p4.2 <- FeaturePlot(scRNA,
                  features = markerGenes, 
                  reduction ="tsne", ncol=4)
ggsave("Cluster/tSNE_markergenes_res.0.1.pdf", p4.2, bg="#ffffff", width=18, height=8)
ggsave("Cluster/tSNE_markergenes_res.0.1.png", p4.2, bg="#ffffff", width=18, height=8)


p5 <- VlnPlot(scRNA, pt.size=0, 
              features = markerGenes, 
              group.by = "RNA_snn_res.0.1", ncol=4) 
ggsave("Cluster/VlnPlot_markergenes_res.0.1.pdf", p5, dpi=300, width=18, height=8)
ggsave("Cluster/VlnPlot_markergenes_res.0.1.png", p5,dpi=300, bg="#ffffff", width=18, height=8)

p6 <- DotPlot(scRNA,  
              features = markerGenes, 
              group.by = "RNA_snn_res.0.1") 
ggsave("Cluster/DotPlot_markergenes_res.0.1.pdf", p6, dpi=300, width=10, height=8)
ggsave("Cluster/DotPlot_markergenes_res.0.1.png", p6, dpi=300, bg="#ffffff", width=10, height=8)


##Save data
saveRDS(scRNA, file="RDS/scRNA.rds")



##Define cell type
dir.create("CellType")

rm(list=ls())
library(SingleR)
library(celldex)

scRNA2 <- readRDS("RDS/scRNA.rds")
DefaultAssay(scRNA2) <- "RNA"

#Reference
refdata <- get(load("~/Reference/celldex/HumanPrimaryCellAtlasData.Rdata"))
testdata <- GetAssayData(scRNA2, slot="data")
clusters <- scRNA2@meta.data$RNA_snn_res.0.1

#Cell annotation
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")


#Check result
cellpred$labels
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)

#Redefine cell types
celltype[celltype$ClusterID %in% c(0,1,2,7),2]='Fib/MF'
celltype[celltype$ClusterID %in% c(3),2]='Pe/SMCs'
celltype[celltype$ClusterID %in% c(4),2]='Immune'
celltype[celltype$ClusterID %in% c(5),2]='Endo'
celltype[celltype$ClusterID %in% c(6),2]='Mesa'
celltype[celltype$ClusterID %in% c(8),2]='Epi'

#Save result
write.csv(celltype,"CellType/Celltype.csv",row.names = F)

#Mapping cell types
scRNA2@meta.data$Celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA2@meta.data[which(scRNA2@meta.data$RNA_snn_res.0.1 == celltype$ClusterID[i]),'Celltype'] <- celltype$celltype[i]
}

p1 = DimPlot(scRNA2, group.by="Celltype", repel=T, 
             label=T, label.size=3, reduction='tsne')
p2 = DimPlot(scRNA2, group.by="Celltype", repel=T, label=T, 
             label.size=3, reduction='umap')
p3 = plotc <- p1+p2+ plot_layout(guides = "collect")
ggsave("CellType/tSNE_celltype.pdf", p1, width=7, height=6)
ggsave("CellType/tSNE_celltype.png", p1, width=7, height=6)
ggsave("CellType/UMAP_celltype.pdf", p2, width=7, height=6)
ggsave("CellType/UMAP_celltype.png", p2, width=7, height=6)
ggsave("CellType/Celltype.pdf", p3, width=14, height=6)
ggsave("CellType/Celltype.png", p3, width=14, height=6)


#Save data
saveRDS(scRNA2, file="RDS/scRNA2.rds")



##Fibroblast reclustering
rm(list=ls())
dir.create("Subcluster")
scRNA3 <- readRDS("RDS/scRNA2.rds")


#Extract cell subset
Cells.sub <- subset(scRNA3@meta.data, Celltype=="Fib/MF") 
scRNAsub_Fb <- subset(scRNA3, cells=row.names(Cells.sub))

#Build a new data set
cellinfo <- subset(scRNAsub_Fb@meta.data, select = c("orig.ident", "nCount_RNA", "nFeature_RNA")) 
scRNAsub <- CreateSeuratObject(scRNAsub_Fb@assays$RNA@counts, meta.data = cellinfo)

rm(scRNA3,Cells.sub,scRNAsub_Fb,cellinfo)


#Data standardization and normalization
scRNAsub <- NormalizeData(scRNAsub, normalization.method = "LogNormalize", scale.factor = 1e4) %>% 
            FindVariableFeatures(selection.method = "vst", Features = 3000, verbose = FALSE) %>% 
            ScaleData()

scRNAsub <- RunPCA(scRNAsub, 
                   features = VariableFeatures(object = scRNAsub), 
                   assay = "RNA")

plot1 <- DimPlot(scRNAsub, shuffle =T, 
                 reduction = "pca", 
                 group.by="orig.ident") 
plot2 <- ElbowPlot(scRNAsub, 
                   ndims=50, 
                   reduction="pca") 
plotc <- plot1+plot2
ggsave("Subcluster/pca.pdf", plot = plotc, width = 8, height = 4) 
ggsave("Subcluster/pca.png", plot = plotc, width = 8, height = 4) 


#Pick PCA number
pc.num=1:40

#Reduction
scRNAsub <- RunUMAP(scRNAsub, 
                    dims = pc.num, 
                    reduction = "pca")

scRNAsub <- RunTSNE(scRNAsub, 
                    dims = pc.num, 
                    reduction = "pca")

#Clustering
scRNAsub <- FindNeighbors(scRNAsub, 
                          dims = pc.num, 
                          reduction = "pca")

scRNAsub <- FindClusters(scRNAsub, graph.name = "RNA_snn",
                         algorithm = 1,
                         resolution =c(seq(0,1,.1)), 
                         verbose = FALSE)
#Check result
library(clustree)
p1 <- clustree(scRNAsub, prefix = "RNA_snn_res.")
ggsave("Subcluster/clustree.pdf", p1, width = 12, height = 16)
ggsave("Subcluster/clustree.png", p1, width = 12, height = 16)


p2 <- DimPlot(scRNAsub,
              reduction = 'umap', 
              group.by = "RNA_snn_res.0.1",  
              repel=T, label =T)
ggsave("Subcluster/UMAP_cluster_res.0.1.pdf", p2, width = 8, height = 8)
ggsave("Subcluster/UMAP_cluster_res.0.1.png", p2, width = 8, height = 8)


p3 <- DimPlot(scRNAsub,
              reduction = 'tsne', 
              group.by = "RNA_snn_res.0.1",  
              repel=T, label =T)
ggsave("Subcluster/tSNE_cluster_res.0.1.pdf", p3, width = 8, height = 8)
ggsave("Subcluster/tSNE_cluster_res.0.1.png", p3, width = 8, height = 8)


#Check cluster quality
theme.set = theme(axis.title.x = element_blank())
plot.features = c("nFeature_RNA", "nCount_RNA")

plots = list()
for(i in seq_along(plot.features)){
  plots[[i]] =  VlnPlot(scRNAsub, pt.size = 0, 
                        group.by = "RNA_snn_res.0.1", 
                        features = plot.features[i]) + 
                        theme.set + NoLegend()
}
p4 <- wrap_plots(plots = plots, nrow=1)
ggsave("Subcluster/vlnplot_cluster_res.0.1.pdf", p4, width = 14, height = 8)
ggsave("Subcluster/vlnplot_cluster_res.0.1.png", p4, width = 14, height = 8)


#Significantly overexpressed genes
DefaultAssay(scRNAsub) <- "RNA" 
Idents(scRNAsub) <- "RNA_snn_res.0.1"
plan("multisession", workers = 8)

diff.wilcox = FindAllMarkers(scRNAsub,
                             only.pos = FALSE, 
                             min.pct = 0, 
                             thresh.use = 0)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(all.markers, "Subcluster/Diff_genes_wilcox_res.0.1_RNA.csv", row.names = F)
write.csv(top10, "Subcluster/Top10_diff_genes_wilcox_res.0.1_RNA.csv", row.names = F)
plan("sequential") 

#Top10 genes
top10_genes = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10_genes = CaseMatch(search = as.vector(top10_genes$gene), match = rownames(scRNAsub)) 
plot1 = DoHeatmap(scRNAsub, 
                  slot = "data",
                  features = top10_genes, 
                  group.by = "RNA_snn_res.0.1", 
                  group.bar = T, size = 4)
ggsave("Subcluster/Top10_markers_res.0.1_RNA.pdf", plot=plot1, width=30, height=40) 
ggsave("Subcluster/Top10_markers_res.0.1_RNA.png", plot=plot1, width=30, height=40)


#Marker gene
select_genes <- c("DCN", "POSTN",  #Fib/MF
                  "CD248", "PTPRC" 
                   )

#Vlnplot
p1 <- VlnPlot(scRNAsub, 
              features = select_genes, pt.size=0, 
              group.by="RNA_snn_res.0.1", ncol=2)
ggsave("Subcluster/Selectgenes_VlnPlot.png", p1, width=12, height=10)

#Featureplot
p2 <- FeaturePlot(scRNAsub, 
                  features = select_genes, 
                  reduction = "umap", ncol=2)
ggsave("Subcluster/Selectgenes_umapplot.png", p2, width=10, height=8)

#DotPlot
p3 <- DotPlot(scRNAsub, 
              features = select_genes, 
              group.by="RNA_snn_res.0.1")
ggsave("Subcluster/Selectgenes_DotPlot.png", p3, bg="#ffffff", width=10, height=8)

#Heatmap
p4 <- DoHeatmap(scRNAsub, 
                features = select_genes, 
                group.by = "RNA_snn_res.0.1", 
                group.bar = T, size = 4, slot="data") 
ggsave("Subcluster/Selectgenes_Heatmap.png", p4, width=10, height=8)



#Redefine subgroups
Idents(scRNAsub) <- "RNA_snn_res.0.1"
new.cluster.ids <- c("F1", "F2", "F3","F4", "F5", "F6", "F7")

names(new.cluster.ids) <- levels(scRNAsub)

#Add information in metadata
scRNAsub <- RenameIdents(scRNAsub, new.cluster.ids)
scRNAsub$original_Cluster <- Idents(scRNAsub)


#Load software
library(SCP)
library(paletteer)
pal <- paletteer_d( "ggsci::nrc_npg")[c(1:20)];pal


CellDimPlot(
  srt = scRNAsub, group.by = c("original_Cluster"),
  reduction = "umap", theme_use = "theme_blank",label_insitu = T,palette = "simpsons",
  label = F,raster = F
)
ggsave("Subcluster/Kidney-Fibrosis_CellDimPlot_label.pdf",width = 6,height = 4.8)
ggsave("Subcluster/Kidney-Fibrosis_CellDimPlot_label.tiff",width = 6,height = 4.8)
ggsave("Subcluster/Kidney-Fibrosis_CellDimPlot_label.svg",width = 6,height = 4.8)


FeatureStatPlot(scRNAsub, stat.by = c("CD248"), 
                group.by = "original_Cluster",palette = "simpsons")+ #palette = "npg"
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1), geom = "pointrange",show.legend = F)+
  stat_summary(fun = "mean", geom = "point",color = "white",show.legend = F)

ggsave("Subcluster/Kidney-Fibrosis_FeatureDimPlot_Cd248.pdf",width = 6,height = 4)
ggsave("Subcluster/Kidney-Fibrosis_FeatureDimPlot_Cd248.tiff",width = 6,height = 4)
ggsave("Subcluster/Kidney-Fibrosis_FeatureDimPlot_Cd248.svg",width = 6,height = 4)


#Save data
saveRDS(scRNAsub, file="RDS/scRNAsub.rds")



##Differential genes expression analysis
library(Seurat)
library(tidyverse)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
rm(list=ls())

dir.create("Subcluster/Enrich")
scRNAsub <- readRDS("RDS/scRNAsub.rds")


#cluster5 
Idents(scRNAsub) <- "RNA_snn_res.0.1"

dge.cluster <- FindMarkers(scRNAsub,
                           ident.1 = 5,
                           only.pos = TRUE)

sig_dge.cluster <- subset(dge.cluster, p_val_adj<0.05&abs(avg_log2FC)>0.5)
write.csv(sig_dge.cluster,'Subcluster/Enrich/sig_dge.cluster.csv')


##GO enrichment analysis
ego_ALL <- enrichGO(gene          = row.names(sig_dge.cluster),
                    OrgDb         = 'org.Hs.eg.db',
                    keyType       = 'SYMBOL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)

ego_all <- data.frame(ego_ALL)
write.csv(ego_all,'Subcluster/Enrich/enrichGO.csv')  

ego_CC <- enrichGO(gene          = row.names(sig_dge.cluster),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_MF <- enrichGO(gene          = row.names(sig_dge.cluster),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_BP <- enrichGO(gene          = row.names(sig_dge.cluster),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)           


#Visual result
ego_CC@result$Description <- substring(ego_CC@result$Description,1,70)
ego_MF@result$Description <- substring(ego_MF@result$Description,1,70)
ego_BP@result$Description <- substring(ego_BP@result$Description,1,70)
p_BP <- barplot(ego_BP,showCategory = 10) + ggtitle("barplot for Biological process")
p_CC <- barplot(ego_CC,showCategory = 10) + ggtitle("barplot for Cellular component")
p_MF <- barplot(ego_MF,showCategory = 10) + ggtitle("barplot for Molecular function")
plota <- p_BP/p_CC/p_MF
ggsave('Subcluster/Enrich/enrichGO_barplot.pdf', plota, width = 10,height = 18)
ggsave('Subcluster/Enrich/enrichGO_barplot.png', plota, width = 10,height = 18)


p_BP2 <- dotplot(ego_BP,showCategory = 10) + ggtitle("dotplot for Biological process")
p_CC2 <- dotplot(ego_CC,showCategory = 10) + ggtitle("dotplot for Cellular component")
p_MF2 <- dotplot(ego_MF,showCategory = 10) + ggtitle("dotplot for Molecular function")
plotb <- p_BP2/p_CC2/p_MF2
ggsave('Subcluster/Enrich/enrichGO_dotplot.pdf', plotb, width = 10,height = 18)
ggsave('Subcluster/Enrich/enrichGO_dotplot.png', plotb, width = 10,height = 18)


#Beautify the picture
GO_BP_up = as.data.frame(ego_BP)

df = GO_BP_up[1:10,]
df$LogP <- -log10(df$pvalue)

df$labelx=rep(0,nrow(df))
df$labely=seq(nrow(df),1)

ggplot(data = df, 
       aes(LogP, reorder(Description,LogP))) +
  geom_bar(stat="identity",
           alpha=1.0,
           fill="#AA2123",color = "black",
           width = 0.8) +
  theme_classic()+
  theme(
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_line(colour = 'black', linewidth = 1),
    axis.text.y = element_text(colour = 'black', size = 15),
    axis.line.x = element_line(colour = 'black', linewidth = 1,lineend = "butt"),
    axis.text.x = element_text(colour = 'black', size = 15),
    axis.ticks.x = element_line(colour = 'black', linewidth = 1),
    axis.title.x = element_text(colour = 'black', size = 12))+
  scale_x_discrete(labels = function() str_wrap(x,width = 30))+
  xlab("-log10(pvalue)")+
  ggtitle("Cd248+FIB GO term")+
  scale_x_continuous(expand = c(0,0))
ggsave(filename = "Subcluster/Enrich/Kidney-Fibrosis_Cd248+FIB_GO_barplot.pdf",width = 12,height = 6)



#KEGG enrichment analysis
genelist <- bitr(row.names(sig_dge.cluster), fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db') ##注意物种更换

genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa')

KEGG <- data.frame(ekegg)
write.csv(KEGG,'Subcluster/Enrich/enrichKEGG.csv') 

p1 <- barplot(ekegg, showCategory=20)
p2 <- dotplot(ekegg, showCategory=20)


ggsave("Subcluster/Enrich/enrichKEGG_barplot.pdf", p1, width = 8, height = 8)
ggsave("Subcluster/Enrich/enrichKEGG_barplot.png", p1, width = 8, height = 8)

ggsave("Subcluster/Enrich/enrichKEGG_dotplot.pdf", p2, width = 8, height = 8)
ggsave("Subcluster/Enrich/enrichKEGG_dotplot.png", p2, width = 8, height = 8)


#Visual result
df = KEGG[1:10,]
df$LogP <- -log10(df$pvalue)

df$labelx=rep(0,nrow(df))
df$labely=seq(nrow(df),1)

ggplot(data = df, 
       aes(LogP, reorder(Description,LogP))) +
  geom_bar(stat="identity",
           alpha=1.0,
           fill="#AA2123",color = "black",
           width = 0.8) +
  theme_classic()+
  theme(
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_line(colour = 'black', linewidth = 1),
    axis.text.y = element_text(colour = 'black', size = 15),
    axis.line.x = element_line(colour = 'black', linewidth = 1,lineend = "butt"),
    axis.text.x = element_text(colour = 'black', size = 15),
    axis.ticks.x = element_line(colour = 'black', linewidth = 1),
    axis.title.x = element_text(colour = 'black', size = 12))+
  scale_x_discrete(labels = function() str_wrap(x,width = 30))+
  xlab("-log10(pvalue)")+
  ggtitle("Cd248+FIB KEGG term")+
  scale_x_continuous(expand = c(0,0))
ggsave(filename = "Subcluster/Enrich/Kidney-Fibrosis_Cd248+FIB_KEGG_barplot.pdf",width = 12,height = 6)



