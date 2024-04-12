#title: single cell analysis
#author: Li Mo
#date: 2024.2.22

'20240222 human heart scRNA seq MI infarct and remote 10X genomics'

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
library(future) # Seurat并行运算的一个包
options(future.globals.maxSize = 100* 1024^3)
rm(list = ls())

setwd("~/MI_data/Human/Result")
dir.create('RDS')


#Set data path
dir = c('~/Rawdata/Infarct/filtered_feature_bc_matrix',
        '~/Rawdata/Remote/filtered_feature_bc_matrix'
       )

names = c('Infarct',"Remote")


#Create seurat
scRNAlist <- list()

for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=names[i],min.cells=3,min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = names[i])
  if(T){
    scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
  }
  if(T){
    scRNAlist[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^RP[SL]")
  }
  if(T){
    HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
    HB.genes <- CaseMatch(HB.genes, rownames(scRNAlist[[i]])) 
    scRNAlist[[i]][["percent.hb"]] <- PercentageFeatureSet(scRNAlist[[i]], features=HB.genes) 
  }
}

names(scRNAlist) <- names

scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])


#Calculate the proportion of fetures and counts
scRNA$log10GenesPerUMI <- log10(scRNA$nFeature_RNA) / log10(scRNA$nCount_RNA)


##Data quality control
dir.create("QC")


theme.set = theme(axis.title.x = element_blank())
plot.features = c("nFeature_RNA", "nCount_RNA","log10GenesPerUMI","percent.mt", "percent.rb", "percent.hb")

#Data quality before quality control
plots = list()
for(i in seq_along(plot.features)){
   plots[[i]] =  VlnPlot(scRNA, 
                         group.by = "orig.ident", 
                         pt.size = 0, 
                         features = plot.features[i]) + theme.set + NoLegend()
}
violin <- wrap_plots(plots = plots, nrow=2)
ggsave("QC/vlnplot_before_qc.pdf", plot = violin, width = 12, height = 8)
ggsave("QC/vlnplot_before_qc.png", plot = violin, width = 12, height = 8)


#Remove the double cells
library(dplyr)
library(scDblFinder)

scRNA <- as.SingleCellExperiment(scRNA)
scRNA <- scDblFinder(scRNA, samples = "orig.ident")

Dbcells=table(scRNA$scDblFinder.sample, scRNA$scDblFinder.class)
write.csv(Dbcells,"QC/scDblFinder_sample_cells.csv",row.names = T)

scRNA <- as.Seurat(scRNA)
scRNA <- subset(scRNA, scDblFinder.class=="singlet")


#Set quality control standards
minGene=300
minUMI=500
minRatio=0.8
maxpctMT=25
maxpctRB=25
maxpctHB=5


#Confirm the maxGene value according to top1% genes
rank_nGene = scRNA$nFeature_RNA[order(scRNA$nFeature_RNA,decreasing = T)]
maxGene = as.numeric(rank_nGene[ceiling(length(rank_nGene)/100)])

#Confirm the maxUMI value according to top2% UMI
rank_nUMI= scRNA$nCount_RNA[order(scRNA$nCount_RNA,decreasing = T)]
maxUMI = as.numeric(rank_nUMI[ceiling(length(rank_nUMI)/50)])

#Data quality after quality control
scRNA <- subset(scRNA, 
                subset = nFeature_RNA > minGene & nFeature_RNA < maxGene &
                  nCount_RNA > minUMI & nCount_RNA < maxUMI &
                  log10GenesPerUMI > minRatio & 
                  percent.mt < maxpctMT &
                  percent.rb < maxpctRB &
                  percent.hb < maxpctHB)

plots = list()
for(i in seq_along(plot.features)){
  plots[[i]] =  VlnPlot(scRNA, group.by = "orig.ident", 
                        pt.size = 0, #需要显示点，可以设置pt.size = 0.01
                        features = plot.features[i]) + theme.set + NoLegend()
}
violin <- wrap_plots(plots = plots, nrow=2)
ggsave("QC/vlnplot_after_qc.pdf", plot = violin, width = 12, height = 8)
ggsave("QC/vlnplot_after_qc.png", plot = violin, width = 12, height = 8)


#Data standardization and normalization
scRNA <-  NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 1e4) %>%
          FindVariableFeatures(selection.method = "vst", nfeatures = 3000, verbose = FALSE) 


##Calculated cell cycle
CaseMatch(c(cc.genes$s.genes,cc.genes$g2m.genes),VariableFeatures(scRNA))
g2m_genes = cc.genes$g2m.genes
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(scRNA))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(scRNA))

scRNA <- CellCycleScoring(object=scRNA,  
                          g2m.features=g2m_genes,  
                          s.features=s_genes,
                          set.ident =T)

scRNA$CC.Difference <- scRNA$S.Score - scRNA$G2M.Score

p1 <- VlnPlot(scRNA, features = c("S.Score", "G2M.Score","CC.Difference"), 
              pt.size=0, group.by="orig.ident", ncol=3)
ggsave("QC/vlnplot_CellCycle.pdf", p1, width=12, height=4)
ggsave("QC/vlnplot_CellCycle.png", p1, width=12, height=4)


##Deletion of mitochondrial and ribosomal genes
scRNA <- scRNA[!grepl("^MT-", rownames(scRNA),ignore.case = T), ]
scRNA <- scRNA[!grepl("^RP[SL]", rownames(scRNA),ignore.case = T), ]


##Corrected cell cycle effects
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 3000, verbose = FALSE) %>% 
         ScaleData(vars.to.regress  = c("CC.Difference"))


##Save data
saveRDS(scRNA, file="RDS/scRNA.rds")



#Reload data
rm(list = ls())
set.seed(12345)
scRNA <- readRDS("RDS/scRNA.rds")


#SCTransform standardization and normalization
scRNA <- SCTransform(scRNA, verbose = F, 
                     variable.features.n = 3000,
                     vars.to.regress = c("CC.Difference"))


scRNA <- RunPCA(scRNA, 
                npcs = 50,
                features = VariableFeatures(object = scRNA), 
                assay = "SCT")

plot1 <- DimPlot(scRNA, shuffle =T, 
                 reduction = "pca", 
                 group.by="orig.ident") 
plot2 <- ElbowPlot(scRNA, 
                   ndims=50, 
                   reduction="pca") 
plota <- plot1+plot2
ggsave("QC/pca.pdf", plot = plota, width = 8, height = 4) 
ggsave("QC/pca.png", plot = plota, width = 8, height = 4) 


##Correction of batch variance
library(harmony)

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
ggsave("QC/Harmony.pdf", plot = plotb, width = 8, height = 4) 
ggsave("QC/Harmony.png", plot = plotb, width = 8, height = 4)


#Pick PCA number
pc.num=1:40
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
ggsave("Cluster/UMAP_orig.ident_int.pdf", p1, width = 8, height = 8)
ggsave("Cluster/UMAP_orig.ident_int.png", p1, width = 8, height = 8)


p2 = DimPlot(scRNA, shuffle =T, 
             reduction = 'tsne',
             raster=FALSE,
             group.by = "orig.ident")
ggsave("Cluster/tSNE_orig.ident_int.pdf", p2, width = 8, height = 8)
ggsave("Cluster/tSNE_orig.ident_int.png", p2, width = 8, height = 8)


#Clustering
scRNA <- FindNeighbors(scRNA, 
                       dims = pc.num, 
                       reduction = "harmony")

scRNA <- FindClusters(scRNA, graph.name = "SCT_snn",
                      algorithm = 1,
                      resolution =c(seq(0,1,.1)), 
                      verbose = FALSE)
#Check result
library(clustree)
p1 <- clustree(scRNA, prefix = "SCT_snn_res.")
ggsave("Cluster/clustree.pdf", p1, width = 12, height = 16)
ggsave("Cluster/clustree.png", p1, width = 12, height = 16)

p2 <- DimPlot(scRNA,
              reduction = 'umap', 
              group.by = "SCT_snn_res.0.4",  
              repel=T, label =T)
ggsave("Cluster/UMAP_cluster_res.0.4.pdf", p2, width = 8, height = 8)
ggsave("Cluster/UMAP_cluster_res.0.4.png", p2, width = 8, height = 8)


p3 <- DimPlot(scRNA,
              reduction = 'tsne', 
              group.by = "SCT_snn_res.0.4",  
              repel=T, label =T)
ggsave("Cluster/tSNE_cluster_res.0.4.pdf", p3, width = 8, height = 8)
ggsave("Cluster/tSNE_cluster_res.0.4.png", p3, width = 8, height = 8)


#Marker gene mapping
DefaultAssay(scRNA) <- "RNA"

markerGenes <- c("COX4I2","PDGFRB", #locate Pericyte (Pe)
                 "DCN","COL1A1",    #locate Fibroblast (Fib)
                 "POSTN",           #locate Myofibroblast (MF)
                 "RERGL","MYH11",   #locate Vascular smooth muscle cells (VSMCs)
                 "KRT5",            #locate Basal Cells
                 "PECAM1","VWF",    #locate Endothelial cells (Endo)
                 "EPCAM","WT1",     #locate Epithelial cells (Epi)
                 "PTPRC",           #locate Immune cells
                 "CD68","C1QA",     #locate Macrophages
                 "FCN1","PLAC8",    #locate Monocytes
                 "TPSAB1",          #locate Mast cells
                 "CLEC10A",         #locate DCs
                 "S100A8",          #locate Granulocyte 
                 "MS4A1","CD79A",   #locate B Cells
                 "CD3D","NKG7"      #locate T/NK Cells                 
)

p4.1 <- FeaturePlot(scRNA,
                  features = markerGenes, 
                  reduction ="umap", ncol=4)
ggsave("Cluster/UMAP_markergenes_res.0.4.pdf", p4.1, bg="#ffffff", width=18, height=16)
ggsave("Cluster/UMAP_markergenes_res.0.4.png", p4.1, bg="#ffffff", width=18, height=16)


p4.2 <- FeaturePlot(scRNA,
                  features = markerGenes, 
                  reduction ="tsne", ncol=4)
ggsave("Cluster/tSNE_markergenes_res.0.4.pdf", p4.2, bg="#ffffff", width=18, height=16)
ggsave("Cluster/tSNE_markergenes_res.0.4.png", p4.2, bg="#ffffff", width=18, height=16)


p5 <- VlnPlot(scRNA, pt.size=0, 
              features = markerGenes, 
              group.by = "SCT_snn_res.0.4", ncol=4) 
ggsave("Cluster/VlnPlot_markergenes_res.0.4.pdf", p5, dpi=300, width=18, height=16)
ggsave("Cluster/VlnPlot_markergenes_res.0.4.png", p5,dpi=300, bg="#ffffff", width=18, height=16)

p6 <- DotPlot(scRNA,  
              features = markerGenes, 
              group.by = "SCT_snn_res.0.4") 
ggsave("Cluster/DotPlot_markergenes_res.0.4.pdf", p6, dpi=300, width=18, height=16)
ggsave("Cluster/DotPlot_markergenes_res.0.4.png", p6, dpi=300, bg="#ffffff", width=18, height=16)


##Save data
saveRDS(scRNA, file="RDS/scRNA_harmony_int.rds")



##Define cell type
dir.create("CellType")

rm(list=ls())
set.seed(12345)
library(SingleR)
library(celldex)

scRNA2 <- readRDS("RDS/scRNA_harmony_int.rds")
DefaultAssay(scRNA2) <- "RNA"

#Reference
refdata <- get(load("/PUBLIC/limo/reference/singlecell/celldex/HumanPrimaryCellAtlasData.Rdata"))
testdata <- GetAssayData(scRNA2, slot="data")
clusters <- scRNA2@meta.data$SCT_snn_res.0.4

#Cell annotation
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")


#Check result
cellpred$labels
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)


#Redefine cell types
celltype[celltype$ClusterID %in% c(0),2]='Pericytes/SMCs'
celltype[celltype$ClusterID %in% c(1,4),2]='Fibroblasts'
celltype[celltype$ClusterID %in% c(2,5,6,10),2]='Endothelial cells'
celltype[celltype$ClusterID %in% c(3),2]='T cells'
celltype[celltype$ClusterID %in% c(7),2]='NK cells'
celltype[celltype$ClusterID %in% c(8),2]='Macrophages'
celltype[celltype$ClusterID %in% c(9),2]='Granulocytes'
celltype[celltype$ClusterID %in% c(11),2]='Undefined'


#Save result
write.csv(celltype,"CellType/Main_Cell_type.csv",row.names = F)

#Mapping cell types
scRNA2@meta.data$Main_Cell_type = "NA"
for(i in 1:nrow(celltype)){
  scRNA2@meta.data[which(scRNA2@meta.data$SCT_snn_res.0.4 == celltype$ClusterID[i]),'Main_Cell_type'] <- celltype$celltype[i]
}

p1 = DimPlot(scRNA2, group.by="Main_Cell_type", repel=T, 
             label=T, label.size=5, reduction='tsne')
p2 = DimPlot(scRNA2, group.by="Main_Cell_type", repel=T, label=T, 
             label.size=5, reduction='umap')
p3 = plotc <- p1+p2+ plot_layout(guides = "collect")
ggsave("CellType/tSNE_Main_Cell_type.pdf", p1, width=8, height=7)
ggsave("CellType/tSNE_Main_Cell_type.png", p1, width=8, height=7)
ggsave("CellType/UMAP_Main_Cell_type.pdf", p2, width=8, height=7)
ggsave("CellType/UMAP_Main_Cell_type.png", p2, width=8, height=7)
ggsave("CellType/Main_Cell_type.pdf", p3, width=16, height=8)
ggsave("CellType/Main_Cell_type.png", p3, width=16, height=8)


#Save data
saveRDS(scRNA2, file="RDS/scRNA_Main_Cell_type.rds")



##Fibroblast reclustering
rm(list=ls())
dir.create("Subcluster")
scRNA3 <- readRDS("RDS/scRNA_Main_Cell_type.rds")


#Extract cell subset
Cells.sub <- subset(scRNA3@meta.data, Main_Cell_type=="Fibroblasts") 
scRNAsub_Fb <- subset(scRNA3, cells=row.names(Cells.sub))

#Build a new data set
cellinfo <- subset(scRNAsub_Fb@meta.data, select = c("orig.ident", "nCount_RNA", "nFeature_RNA","orig.ident")) 
scRNAsub <- CreateSeuratObject(scRNAsub_Fb@assays$RNA@counts, meta.data = cellinfo)


#Data standardization and normalization
scRNAsub <- NormalizeData(scRNAsub, normalization.method = "LogNormalize", scale.factor = 1e4) %>% 
            FindVariableFeatures(selection.method = "vst", Features = 3000, verbose = FALSE) 


##Calculated cell cycle
CaseMatch(c(cc.genes$s.genes,cc.genes$g2m.genes),VariableFeatures(scRNAsub))

g2m_genes = cc.genes$g2m.genes
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(scRNAsub))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(scRNAsub))

scRNAsub <- CellCycleScoring(object=scRNAsub,  
                          g2m.features=g2m_genes,  
                          s.features=s_genes,
                          set.ident =T)

scRNAsub$CC.Difference <- scRNAsub$S.Score - scRNAsub$G2M.Score

p1 <- VlnPlot(scRNAsub, features = c("S.Score", "G2M.Score","CC.Difference"), 
              pt.size=0, group.by="orig.ident", ncol=3)
ggsave("Subcluster/vlnplot_CellCycle.pdf", p1, width=12, height=4)
ggsave("Subcluster/vlnplot_CellCycle.png", p1, width=12, height=4)


##Corrected cell cycle effects
scRNAsub <- FindVariableFeatures(scRNAsub, selection.method = "vst", nfeatures = 3000, verbose = FALSE) %>% 
            ScaleData(vars.to.regress  = c("CC.Difference"))


set.seed(12345)

#SCTransform standardization and normalization
scRNAsub <- SCTransform(scRNAsub, verbose = F, 
                     variable.features.n = 3000,
                     vars.to.regress = c("CC.Difference"))


scRNAsub <- RunPCA(scRNAsub, 
                npcs = 50,
                features = VariableFeatures(object = scRNAsub), 
                assay = "SCT")

plot1 <- DimPlot(scRNAsub, shuffle =T, 
                 reduction = "pca", 
                 group.by="orig.ident") 
plot2 <- ElbowPlot(scRNAsub, 
                   ndims=50, 
                   reduction="pca") 
plota <- plot1+plot2
ggsave("Subcluster/pca.pdf", plot = plota, width = 8, height = 4) 
ggsave("Subcluster/pca.png", plot = plota, width = 8, height = 4) 


#Correction of batch variance
scRNAsub <- RunHarmony(scRNAsub,
                      max.iter.harmony = 20, 
                      group.by.vars = "orig.ident", 
                      assay.use = "RNA")
                  
names(scRNAsub@reductions)

plot3 <- DimPlot(scRNAsub,shuffle =T,
                 reduction = "harmony", 
                 group.by="orig.ident") 
plot4 <- ElbowPlot(scRNAsub, 
                   ndims=50, 
                   reduction="harmony") 
plotb <- plot3+plot4
ggsave("Subcluster/Harmony.pdf", plot = plotb, width = 8, height = 4) 
ggsave("Subcluster/Harmony.png", plot = plotb, width = 8, height = 4)


#Pick PCA number
pc.num=1:40

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
ggsave("Subcluster/UMAP_orig.ident_int.pdf", p1, width = 8, height = 8)
ggsave("Subcluster/UMAP_orig.ident_int.png", p1, width = 8, height = 8)

p2 = DimPlot(scRNAsub, shuffle =T, 
             reduction = 'tsne', 
             group.by = "orig.ident")
ggsave("Subcluster/tSNE_orig.ident_int.pdf", p2, width = 8, height = 8)
ggsave("Subcluster/tSNE_orig.ident_int.png", p2, width = 8, height = 8)


#Clustering
scRNAsub <- FindNeighbors(scRNAsub, 
                          dims = pc.num, 
                          reduction = "harmony")

scRNAsub <- FindClusters(scRNAsub, graph.name = "SCT_snn",
                         algorithm = 1,
                         resolution =c(seq(0,1,.1)), 
                         verbose = FALSE)
#Check result
library(clustree)
p1 <- clustree(scRNAsub, prefix = "SCT_snn_res.")
ggsave("Subcluster/clustree.pdf", p1, width = 12, height = 16)
ggsave("Subcluster/clustree.png", p1, width = 12, height = 16)


p2 <- DimPlot(scRNAsub,
              reduction = 'umap', 
              group.by = "SCT_snn_res.0.3",  
              repel=T, label =T)
ggsave("Subcluster/UMAP_cluster_res.0.3.pdf", p2, width = 8, height = 8)
ggsave("Subcluster/UMAP_cluster_res.0.3.png", p2, width = 8, height = 8)


p3 <- DimPlot(scRNAsub,
              reduction = 'tsne', 
              group.by = "SCT_snn_res.0.3",  
              repel=T, label =T)
ggsave("Subcluster/tSNE_cluster_res.0.3.pdf", p3, width = 8, height = 8)
ggsave("Subcluster/tSNE_cluster_res.0.3.png", p3, width = 8, height = 8)


#Check cluster quality
theme.set = theme(axis.title.x = element_blank())
plot.features = c("nFeature_RNA", "nCount_RNA")

plots = list()
for(i in seq_along(plot.features)){
  plots[[i]] =  VlnPlot(scRNAsub, pt.size = 0, 
                        group.by = "SCT_snn_res.0.3", 
                        features = plot.features[i]) + 
                        theme.set + NoLegend()
}
p4 <- wrap_plots(plots = plots, nrow=1)
ggsave("Subcluster/vlnplot_cluster_res.0.3.pdf", p4, width = 12, height = 6)
ggsave("Subcluster/vlnplot_cluster_res.0.3.png", p4, width = 12, height = 6)


#Significantly overexpressed genes
DefaultAssay(scRNAsub) <- "RNA" 
Idents(scRNAsub) <- "SCT_snn_res.0.3"
plan("multisession", workers = 10)

diff.wilcox = FindAllMarkers(scRNAsub,
                             only.pos = FALSE, 
                             min.pct = 0, 
                             thresh.use = 0)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(all.markers, "Subcluster/Diff_genes_wilcox_res.0.3_RNA.csv", row.names = F)
write.csv(top10, "Subcluster/Top10_diff_genes_wilcox_res.0.3_RNA.csv", row.names = F)
plan("sequential") 

#Top10 genes
top10_genes = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10_genes = CaseMatch(search = as.vector(top10_genes$gene), match = rownames(scRNAsub)) 
plot1 = DoHeatmap(scRNAsub, 
                  slot = "data",
                  features = top10_genes, 
                  group.by = "SCT_snn_res.0.3", 
                  group.bar = T, size = 4)
ggsave("Subcluster/Top10_markers_res.0.3_RNA.pdf", plot=plot1, width=30, height=40) 
ggsave("Subcluster/Top10_markers_res.0.3_RNA.png", plot=plot1, width=30, height=40)


#Marker gene
select_genes <- c("COL1A1", "DCN", #Fibroblast
                  "CD248", "PTPRC" 
                   )

#Vlnplot
p1 <- VlnPlot(scRNAsub, features = select_genes, pt.size=0, 
              group.by="SCT_snn_res.0.3", ncol=2)
ggsave("Subcluster/Selectgenes_VlnPlot.png", p1, width=12, height=12)

#Featureplot
p2 <- FeaturePlot(scRNAsub, 
                  features = select_genes, 
                  reduction = "umap", ncol=2)
ggsave("Subcluster/Selectgenes_UMAPplot.png", p2, width=12, height=12)

#DotPlot
p3 <- DotPlot(scRNAsub, 
              features = select_genes, 
              group.by="SCT_snn_res.0.3")
ggsave("Subcluster/Selectgenes_DotPlot.png", p3, bg="#ffffff", width=8, height=4)

#Heatmap
p4 <- DoHeatmap(scRNAsub, 
                features = select_genes, 
                group.by = "SCT_snn_res.0.3", 
                group.bar = T, size = 4, slot="data") 
ggsave("Subcluster/Selectgenes_Heatmap.png", p4, width=8, height=4)



#Visualization
library(SCP)
library(ggplot2) 
library(tidyverse) 
library(patchwork) 
library(BiocParallel)
dir.create("Subcluster/SCP")

#Redefine subgroups
Idents(scRNAsub) <- "SCT_snn_res.0.3"
new.cluster.ids <- c("hF1", "hF2", "hF3","hF4", "hF5", "hF6", "hF7")


names(new.cluster.ids) <- levels(scRNAsub)

#Add information in metadata
scRNAsub <- RenameIdents(scRNAsub, new.cluster.ids)
scRNAsub$original_Cluster <- Idents(scRNAsub)


##Umap plot of samples
p1 <- CellDimPlot(srt = scRNAsub, 
                    group.by = "orig.ident",
                    reduction = "umap", 
                    theme_use = "theme_blank")
ggsave("Subcluster/SCP/Human_MI_CellDimPlot_sample.pdf", p1, width = 6, height = 4.8) 
ggsave("Subcluster/SCP/Human_MI_CellDimPlot_sample.tiff", p1, width = 6, height = 4.8) 


##Umap plot of cell types
p2 <- CellDimPlot(srt = scRNAsub, 
                    group.by = "original_Cluster",
                    reduction = "umap", 
                    theme_use = "theme_blank",
                    label_insitu = T,palette = "simpsons",
                    label = F,raster = F
)
ggsave("Subcluster/SCP/Human_MI_CellDimPlot_Cluster.pdf", p2, width = 6,height = 4.8)
ggsave("Subcluster/SCP/Human_MI_CellDimPlot_Cluster.tiff", p2, width = 6,height = 4.8)


#Umap plot of cell types in different position
p2.3 <- CellDimPlot(srt = scRNAsub, 
                    group.by = "original_Cluster",
                    split.by = "orig.ident",
                    reduction = "umap", 
                    theme_use = "theme_blank",
                    label_insitu = T,palette = "simpsons",
                    label = F,raster = F
)
ggsave("Subcluster/SCP/Human_MI_CellDimPlot_split.pdf", p2.3, width = 12, height = 4.8)
ggsave("Subcluster/SCP/Human_MI_CellDimPlot_split.tiff", p2.3, width = 12, height = 4.8)


#Marker gene display
p3.3 <- FeatureStatPlot(scRNAsub, 
                        stat.by = c("CD248"), 
                        group.by = "original_Cluster",
                        palette = "simpsons")+
        stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1), geom = "pointrange",show.legend = F)+
        stat_summary(fun = "mean", geom = "point",color = "white",show.legend = F)
 
ggsave("Subcluster/SCP/Human_MI_FeatureStatPlot.pdf",width = 6,height = 4)
ggsave("Subcluster/SCP/Human_MI_FeatureStatPlot.tiff",width = 6,height = 4)


##Differential gene analysis
scRNAsub <- RunDEtest(srt = scRNAsub, 
                      group_by = "original_Cluster", 
                      fc.threshold = 1, 
                      only.pos = FALSE)

DEGs <- scRNAsub@tools$DEtest_original_Cluster$AllMarkers_wilcox
DEGs <- DEGs[with(DEGs, avg_log2FC > 0.5 & p_val_adj < 0.05), ]
write.csv(DEGs,'Subcluster/SCP/DEGs.csv') 

top10_DEGs <- DEGs %>% group_by(group1) %>% top_n(n = 10, wt = avg_log2FC)


##Annotate top 10 genes
markers <- c("SOCS3","GADD45B","MIR23AHG","COL3A1","COL1A2","POSTN",
             "ABCA8","SCN7A","COL15A1","SFRP1","PTGDS","NTRK2",
             "PCOLCE2","IGFBP6","S100A10","APOE","APOC1","GGT5",
             "TNMD","KCNMA1","COMP")


#Top 10 genes Heatmap
ht1 <- GroupHeatmap(srt = scRNAsub,
                    features = top10_DEGs$gene,
                    group.by = "original_Cluster",
                    group_palette = "simpsons",
                    show_column_names = TRUE,
                    features_label = markers,
                    heatmap_palcolor = c("navy", "white", "firebrick3"),
)

p6 <- print(ht1$plot)
ggsave("Subcluster/SCP/DEGs_Top10_GroupHeatmap.pdf", p6, width = 8, height = 8) 
ggsave("Subcluster/SCP/DEGs_Top10_GroupHeatmap.png", p6, width = 8, height = 8) 


#Save data
saveRDS(scRNAsub, file="RDS/scRNAsub.rds")



##Differential genes expression analysis
library(Seurat)
library(tidyverse)
library(patchwork)
library(monocle)
library(clusterProfiler)
library(org.Hs.eg.db)
rm(list=ls())

dir.create("Subcluster/Enrich")
scRNAsub <- readRDS("RDS/scRNAsub.rds")

Idents(scRNAsub) <- "SCT_snn_res.0.3"

#cluster4
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
ggsave(filename = "Subcluster/Enrich/Human_MI_Cd248+FIB_GO_barplot.pdf",width = 12,height = 6)



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
ggsave(filename = "Subcluster/Enrich/Human_MI_Cd248+FIB_KEGG_barplot.pdf",width = 12,height = 8)
