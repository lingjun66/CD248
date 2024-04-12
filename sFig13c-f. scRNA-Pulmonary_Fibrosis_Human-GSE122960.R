#title: single cell analysis-GSE122960
#author: Li Mo
#date: 2024.2.16

'2019 Single-Cell Analysis of Human Pulmonary Fibrosis 10X genomics'


#Load software
library(Seurat)
library(SingleR)
library(celldex)
library(tidyverse)
library(cowplot)
library(patchwork)
library(R.utils)
library(data.table)
library(pacman)
library(future) 
options(future.globals.maxSize = 100* 1024^3)
rm(list = ls())

setwd("～/GSE122960/Result")
dir.create('RDS')

#Set data path
dir='～/GSE122960/Rawdata' 

#Load data
h5_files <- list.files(dir, pattern = "\\_filtered_gene_bc_matrices_h5.h5$")

#Create seurat
seurat_list <- list()

for (h5_file in h5_files) {
  data.path <- paste0(dir, h5_file)
  seurat_data <- Read10X_h5(filename = data.path)
  sample_name <- sub("_.*", "", h5_file)
  seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   project = sample_name,
                                   min.features = 200,
                                   min.cells = 3)
  
  seurat_list <- append(seurat_list, seurat_obj)
}

sample_names <- sub("_.*", "", h5_files)

scRNA <- merge(seurat_list[[1]],
               y = seurat_list[-1],
               add.cell.ids = sample_names)

rm(seurat_data,seurat_obj,seurat_list)


#Calculate the proportion of fetures and counts
scRNA$log10GenesPerUMI <- log10(scRNA$nFeature_RNA) / log10(scRNA$nCount_RNA)


#Data quality control
dir.create("QC")

theme.set = theme(axis.title.x = element_blank())
plot.features = c("nFeature_RNA", "nCount_RNA", "log10GenesPerUMI")

#Data quality before quality control
plots = list()
for(i in seq_along(plot.features)){
   plots[[i]] =  VlnPlot(scRNA, group.by = "orig.ident", 
                  pt.size = 0,
                  features = plot.features[i]) + theme.set + NoLegend()
}
violin <- wrap_plots(plots = plots, nrow=1)
ggsave("QC/vlnplot_before_qc.pdf", plot = violin, width = 20, height = 4)
ggsave("QC/vlnplot_before_qc.png", plot = violin, width = 20, height = 4)



#Set quality control standards
minGene=200
minUMI=500


#Data quality after quality control
scRNA <- subset(scRNA, 
                subset = nFeature_RNA > minGene & 
                nCount_RNA > minUMI &
                log10GenesPerUMI > 0.8)

plots = list()
for(i in seq_along(plot.features)){
  plots[[i]] =  VlnPlot(scRNA, group.by = "orig.ident", 
                        pt.size = 0, #需要显示点，可以设置pt.size = 0.01
                        features = plot.features[i]) + theme.set + NoLegend()
}
violin <- wrap_plots(plots = plots, nrow=1)
ggsave("QC/vlnplot_after_qc.pdf", plot = violin, width = 20, height = 4)
ggsave("QC/vlnplot_after_qc.png", plot = violin, width = 20, height = 4)


#Data standardization and normalization
set.seed(1000)

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

##Save data
saveRDS(scRNA, file="RDS/scRNA.rds")



#Correction of batch variance
library(harmony)
rm(list = ls())

scRNA <- readRDS("RDS/scRNA.rds")


scRNA <- RunHarmony(scRNA,
                    max.iter.harmony = 20, 
                    group.by.vars = "orig.ident")


plot3 <- DimPlot(scRNA,shuffle =T,
                 raster=FALSE,
                 reduction = "harmony", 
                 group.by="orig.ident") 
plot4 <- ElbowPlot(scRNA, 
                   ndims=50, 
                   reduction="harmony") 
plotb <- plot3+plot4
ggsave("QC/Harmony_int.pdf", plot = plotb, width = 10, height = 4) 
ggsave("QC/Harmony_int.png", plot = plotb, width = 10, height = 4)


#Pick PCA number
pc.num=1:20

dir.create("Cluster") 

##Reduction
scRNA <- RunUMAP(scRNA, 
                 dims = pc.num, 
                 reduction = "harmony")
scRNA <- RunTSNE(scRNA, 
                  dims = pc.num, 
                  reduction = "harmony")

p1 = DimPlot(scRNA, shuffle =T, 
             reduction = 'umap',
             raster=FALSE,
             group.by = "orig.ident")
ggsave("Cluster/UMAP_orig.ident_int.pdf", p1, width = 10, height = 8)
ggsave("Cluster/UMAP_orig.ident_int.png", p1, width = 10, height = 8)

p2 = DimPlot(scRNA, shuffle =T, 
             reduction = 'tsne',
             raster=FALSE,
             group.by = "orig.ident")
ggsave("Cluster/tSNE_orig.ident_int.pdf", p2, width = 10, height = 8)
ggsave("Cluster/tSNE_orig.ident_int.png", p2, width = 10, height = 8)


#Clustering
scRNA <- FindNeighbors(scRNA, 
                       dims = pc.num, 
                       reduction = "harmony")

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
              group.by = "RNA_snn_res.0.3",  
              repel=T, label =T)
ggsave("Cluster/UMAP_cluster_res.0.3.pdf", p2, width = 8, height = 8)
ggsave("Cluster/UMAP_cluster_res.0.3.png", p2, width = 8, height = 8)


p3 <- DimPlot(scRNA,
              reduction = 'tsne', 
              group.by = "RNA_snn_res.0.3",  
              repel=T, label =T)
ggsave("Cluster/tSNE_cluster_res.0.3.pdf", p3, width = 8, height = 8)
ggsave("Cluster/tSNE_cluster_res.0.3.png", p3, width = 8, height = 8)



#Marker gene mapping
DefaultAssay(scRNA) <- "RNA"


markerGenes <- c("AGER", "SFTPC", #locate Alveolar Type  Cells
                 "COL1A1","DCN",  #locate Fibroblasts
                 "VWF",           #locate Endothelial/Lymphatic Cells
                 "TPPP3",         #locate Ciliated Cells
                 "SCGB3A2",       #locate Club Cells
                 "KRT5",          #locate Basal Cells
                 "CD68",          #locate Macrophages
                 "FCN1",          #locate Monocytes
                 "CLEC10A",       #locate DCs
                 "TPSB2",         #locate Mast Cells
                 "IGHG4",         #locate Plasma Cells
                 "MS4A1",         #locate B Cells
                 "CD3D"           #locate T Cells
)

p4.1 <- FeaturePlot(scRNA,
                  features = markerGenes, 
                  reduction ="umap", ncol=4)
ggsave("Cluster/UMAP_markergenes_res.0.3.pdf", p4.1, bg="#ffffff", width=18, height=16)
ggsave("Cluster/UMAP_markergenes_res.0.3.png", p4.1, bg="#ffffff", width=18, height=16)


p4.2 <- FeaturePlot(scRNA,
                  features = markerGenes, 
                  reduction ="tsne", ncol=4)
ggsave("Cluster/tSNE_markergenes_res.0.3.pdf", p4.2, bg="#ffffff", width=18, height=16)
ggsave("Cluster/tSNE_markergenes_res.0.3.png", p4.2, bg="#ffffff", width=18, height=16)


p5 <- VlnPlot(scRNA, pt.size=0, 
              features = markerGenes, 
              group.by = "RNA_snn_res.0.3", ncol=4) 
ggsave("Cluster/VlnPlot_markergenes_res.0.3.pdf", p5, dpi=300, width=18, height=16)
ggsave("Cluster/VlnPlot_markergenes_res.0.3.png", p5,dpi=300, bg="#ffffff", width=18, height=16)

p6 <- DotPlot(scRNA,  
              features = markerGenes, 
              group.by = "RNA_snn_res.0.3") 
ggsave("Cluster/DotPlot_markergenes_res.0.3.pdf", p6, dpi=300, width=18, height=16)
ggsave("Cluster/DotPlot_markergenes_res.0.3.png", p6, dpi=300, bg="#ffffff", width=18, height=16)


##Save data
saveRDS(scRNA, file="RDS/scRNA_harmony_int.rds")



##Define cell type
dir.create("CellType")

rm(list=ls())
library(SingleR)
library(celldex)

scRNA2 <- readRDS("RDS/scRNA_harmony_int.rds")
DefaultAssay(scRNA2) <- "RNA"

#Reference
refdata <- get(load("/Users/limo/Documents/Tutorials/Single_cell/Reference/celldex/HumanPrimaryCellAtlasData.Rdata"))
testdata <- GetAssayData(scRNA2, slot="data")
clusters <- scRNA2@meta.data$RNA_snn_res.0.3

#Cell annotation
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")


#Check result
cellpred$labels
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)

#Redefine cell types
celltype[celltype$ClusterID %in% c(0),2]='AT2 Cells'
celltype[celltype$ClusterID %in% c(1,2,8,9,12),2]='Macrophages'
celltype[celltype$ClusterID %in% c(3),2]='Monocytes' 
celltype[celltype$ClusterID %in% c(4),2]='Basal Cells'
celltype[celltype$ClusterID %in% c(5),2]='T Cells'
celltype[celltype$ClusterID %in% c(6),2]='AT1 Cells'
celltype[celltype$ClusterID %in% c(7),2]='Ciliated Cells'
celltype[celltype$ClusterID %in% c(10,16),2]='Endothelial Cells'
celltype[celltype$ClusterID %in% c(11),2]='B Cells'
celltype[celltype$ClusterID %in% c(13),2]='Fibroblasts'
celltype[celltype$ClusterID %in% c(14),2]='Mast cells'
celltype[celltype$ClusterID %in% c(15),2]='Club Cells'


#Save result
write.csv(celltype,"CellType/Celltype.csv",row.names = F)

#Mapping cell types
scRNA2@meta.data$Celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA2@meta.data[which(scRNA2@meta.data$RNA_snn_res.0.3 == celltype$ClusterID[i]),'Celltype'] <- celltype$celltype[i]
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
Cells.sub <- subset(scRNA3@meta.data, Celltype=="Fibroblasts") 
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


#Correction of batch variance
DefaultAssay(scRNAsub) = 'RNA'

scRNAsub <- RunHarmony(scRNAsub,
                      max.iter.harmony = 20, 
                      group.by.vars = "orig.ident", 
                      assay.use = "RNA")

plot3 <- DimPlot(scRNAsub,shuffle =T,
                 reduction = "harmony", 
                 group.by="orig.ident") 
plot4 <- ElbowPlot(scRNAsub, 
                   ndims=50, 
                   reduction="harmony") 
plotb <- plot3+plot4
ggsave("Subcluster/Harmony_int.pdf", plot = plotb, width = 8, height = 4) 
ggsave("Subcluster/Harmony_int.png", plot = plotb, width = 8, height = 4)


#Pick PCA number
pc.num=1:20

#Reduction
scRNAsub <- RunUMAP(scRNAsub, 
                    dims = pc.num, 
                    reduction = "harmony")
scRNAsub <- RunTSNE(scRNAsub, 
                    dims = pc.num, 
                    reduction = "harmony")

#Chart showing
p1 = DimPlot(scRNAsub, shuffle =T, 
             reduction = 'umap', 
             group.by = "orig.ident")
ggsave("Subcluster/UMAP_orig.ident_int.pdf", p1, width = 10, height = 8)
ggsave("Subcluster/UMAP_orig.ident_int.png", p1, width = 10, height = 8)

p2 = DimPlot(scRNAsub, shuffle =T, 
             reduction = 'tsne', 
             group.by = "orig.ident")
ggsave("Subcluster/tSNE_orig.ident_int.pdf", p2, width = 10, height = 8)
ggsave("Subcluster/tSNE_orig.ident_int.png", p2, width = 10, height = 8)


#Clustering
scRNAsub <- FindNeighbors(scRNAsub, 
                          dims = pc.num, 
                          reduction = "harmony")

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
plan("multisession", workers = 10)

diff.wilcox = FindAllMarkers(scRNAsub,
                             only.pos = FALSE, 
                             min.pct = 0, 
                             thresh.use = 0)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(all.markers, "Subcluster/Diff_genes_wilcox_res.0.1_RNA.csv", row.names = F)
write.csv(top10, "Subcluster/Top10_diff_genes_wilcox_res.0.1_RNA.csv", row.names = F)
plan("sequential") 

##Top10 genes
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
select_genes <- c("COL1A1", "DCN", #Fibroblast
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
ggsave("Subcluster/Pulmonary-Fibrosis_CellDimPlot_label.pdf",width = 6,height = 4.8)
ggsave("Subcluster/Pulmonary-Fibrosis_CellDimPlot_label.tiff",width = 6,height = 4.8)
ggsave("Subcluster/Pulmonary-Fibrosis_CellDimPlot_label.svg",width = 6,height = 4.8)

FeatureStatPlot(scRNAsub, stat.by = c("CD248"), 
                group.by = "original_Cluster",palette = "simpsons")+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1), geom = "pointrange",show.legend = F)+
  stat_summary(fun = "mean", geom = "point",color = "white",show.legend = F)

ggsave("Subcluster/Pulmonary-Fibrosis_FeatureDimPlot_Cd248.pdf",width = 6,height = 4)
ggsave("Subcluster/Pulmonary-Fibrosis_FeatureDimPlot_Cd248.tiff",width = 6,height = 4)
ggsave("Subcluster/Pulmonary-Fibrosis_FeatureDimPlot_Cd248.svg",width = 6,height = 4)


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


#cluster4
Idents(scRNAsub) <- "RNA_snn_res.0.1"

dge.cluster <- FindMarkers(scRNAsub,
                           ident.1 = 4,
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
plotb <- p_BP2|p_CC2|p_MF2
ggsave('Subcluster/Enrich/enrichGO_dotplot.pdf', plotb, width = 18,height = 6)
ggsave('Subcluster/Enrich/enrichGO_dotplot.png', plotb, width = 18,height = 6)


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
ggsave(filename = "Subcluster/Enrich/Pulmonary-Fibrosis_Cd248+FIB_GO_barplot.pdf",width = 12,height = 6)



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
df = KEGG[1:6,]
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
ggsave(filename = "Subcluster/Enrich/Pulmonary-Fibrosis_Cd248+FIB_KEGG_barplot.pdf",width = 12,height = 4)
