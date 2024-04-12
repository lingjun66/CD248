#---
#  title: "scRNA-seq for singlecell transcriptome in Heart Fibroblast"
#author: "Fu zaiyang"
#date: "2024-02-13"
#---

# load R package ---------------------------------------------------------------
library(SeuratObject)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(data.table)
library(paletteer) 
library(vcd)
library(gplots)
library(harmony)
library(clustree)
library(SingleR)
library(celldex)
library(readr)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(paletteer)
library(dplyr) 
library(tidyverse)
library(presto)
library(SCP)
library(BiocParallel)
library(monocle)
library(ComplexHeatmap)
# initial set ------------------------------------------------------------------

setwd(dir = "/new_data/zyfu/Project_sc_RNA_CF")
DefaultAssay(sce) = "RNA"
set.seed(12345)
register(MulticoreParam(workers = 3, progressbar = TRUE,)) 
pal <- paletteer_d("ggsci::nrc_npg")[c(1:20)];pal
pal <- paletteer_d("ggsci::springfield_simpsons")[c(1:16)];pal

# load data --------------------------------------------------------------------
sce = readRDS("Rdata/heart_3_7_Fibroblasts.rds")

table(sce$original_Cluster)
#F1   F2   F3   F4   F5   F6   F7   F8   F9  F10  F11  F12  F13 
#8163 6508 6460 8620 7467 6279 5869 4856 4835 4007 3326 1688 1627

table(sce$label)
#F1:Ground state Fib 1 F2:Ground state Fib 2     F3:Upk3b/Msln Fib      F4:Ccl2/Ptx3 Fib    F5:Saa3/Cxcl13 Fib 
#8163                  6508                  6460                  8620                  7467 
#F6:Apoe/Fmo2 Fib    F7:Thbs4/Postn Fib F8:Matrifibrocyte Fib          F9:Cd248 Fib   F10:Cthrc1/Ccl9 Fib 
#6279                  5869                  4856                  4835                  4007 
#F11:Lars2/Hspg2 Fib F12:Myofibroblast Fib     F13:Atf3/Fosb Fib 
#3326                  1688                  1627


# Umap-------------------------------------------------------------------
# group by cluster
CellDimPlot(
  srt = sce, group.by = c("label"),
  reduction = "umap", theme_use = "theme_blank",label_insitu = T,palette = "simpsons",
  label = F,raster = F#,split.by = c("Time")
)
ggsave("figures_Cd248/figure1_CellDimPlot_label.pdf",width = 6,height = 4.8)
ggsave("figures_Cd248/figure1_CellDimPlot_label.tiff",width = 6,height = 4.8)
ggsave("figures_Cd248/figure1_CellDimPlot_label.svg",width = 6,height = 4.8)

# group by Time
CellDimPlot(
  srt = sce, group.by = c("Time"),
  reduction = "umap", theme_use = "theme_blank",label_insitu = T,palette = "npg",
  label = F,raster = F#,split.by = c("Time")
)
ggsave("figures_Cd248/figure1_CellDimPlot_Time.pdf",width = 6,height = 4.8)
ggsave("figures_Cd248/figure1_CellDimPlot_Time.tiff",width = 6,height = 4.8)
ggsave("figures_Cd248/figure1_CellDimPlot_Time.svg",width = 6,height = 4.8)

#group by cluster split by Time
CellDimPlot(
  srt = sce, group.by = c("label"),split.by = "Time",
  reduction = "umap", theme_use = "theme_blank",label_insitu = T,palette = "simpsons",
  label = F,raster = F#,split.by = c("Time")
)
ggsave("figures_Cd248/figure1_CellDimPlot_labelsplit.pdf",width = 18,height = 9.6)
ggsave("figures_Cd248/figure1_CellDimPlot_labelsplit.tiff",width = 18,height = 9.6)
ggsave("figures_Cd248/figure1_CellDimPlot_labelsplit.svg",width = 18,height = 9.6)


# Proportion of each Fibroblast cluster -----------------------------------------------------------------
CellStatPlot(sce, stat.by = "label", group.by = "Time",label = F,palette = "simpsons")
ggsave("figures_Cd248/figure1_CellStatPlot.pdf",width = 4.8,height = 4.8)
ggsave("figures_Cd248/figure1_CellStatPlot.tiff",width = 4.8,height = 4.8)
ggsave("figures_Cd248/figure1_CellStatPlot.svg",width = 4.8,height = 4.8)

# Cd248 expression plotted in Umap embedding ----------------------------------------------------------
FeatureDimPlot(sce, features = c("Cd248"),split.by = "Time",reduction = "umap", label = F,
               lower_quantile = 0.1,upper_quantile  = 0.95,raster = T,raster.dpi = c(250,250),show_stat = T)
ggsave("figures_Cd248/figure1_FeatureDimPlot_Cd248_split_by_time.pdf",width = 6*3,height = 4.8*2)
ggsave("figures_Cd248/figure1_FeatureDimPlot_Cd248_split_by_time.tiff",width = 6*3,height = 4.8*2)
ggsave("figures_Cd248/figure1_FeatureDimPlot_Cd248_split_by_time.svg",width = 6*3,height = 4.8*2)

#voilin plots ------------------------------------------------------------------

# time point
FeatureStatPlot(sce, stat.by = c("Cd248"), group.by = "Time"
                ,palette = "npg")+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1), geom = "pointrange",show.legend = F)+
  stat_summary(fun = "mean", geom = "point",color = "white",show.legend = F)
ggsave("figures_Cd248/figure1_FeatureStatPlot_Cd248.pdf",width = 5.4,height = 3)
ggsave("figures_Cd248/figure1_FeatureStatPlot_Cd248.tiff",width = 5.4,height = 3)
ggsave("figures_Cd248/figure1_FeatureStatPlot_Cd248.svg",width = 5.4,height = 3)

# fibroblast cluster
FeatureStatPlot(sce, stat.by = c("Cd248"), group.by = "original_Cluster"
                ,palette = "simpsons")+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1), geom = "pointrange",show.legend = F)+
  stat_summary(fun = "mean", geom = "point",color = "white",show.legend = F)
ggsave("figures_Cd248/figure1_FeatureStatPlot_Cd248_in_cluster.pdf",width = 9.6,height = 3)
ggsave("figures_Cd248/figure1_FeatureStatPlot_Cd248_in_cluster.tiff",width = 9.6,height = 3)
ggsave("figures_Cd248/figure1_FeatureStatPlot_Cd248_in_cluster.svg",width = 9.6,height = 3)

FeatureStatPlot(sce, stat.by = c("Ackr3"), group.by = "original_Cluster"
                ,palette = "simpsons")+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1), geom = "pointrange",show.legend = F)+
  stat_summary(fun = "mean", geom = "point",color = "white",show.legend = F)
ggsave("figures_Cd248/figure1_FeatureStatPlot_Ackr3_in_cluster.pdf",width = 9.6,height = 3)
ggsave("figures_Cd248/figure1_FeatureStatPlot_Ackr3_in_cluster.tiff",width = 9.6,height = 3)
ggsave("figures_Cd248/figure1_FeatureStatPlot_Ackr3_in_cluster.svg",width = 9.6,height = 3)


# Top expressed marker genes of fibroblast cluster-----------------------------------------------------------
features = c("Mt1","Tcf21",#1
             "Smoc2","Gsn",#2
             "Upk3b","Msln",#4
             "Ccl2","Ptx3",#4
             "Saa3","Serpinb2",#5
             "Apoe","Fmo2",#6
             "Thbs4","Postn",#7
             "Comp","Sfrp2",#8
             "Cd248","Ackr3","Pi16",#9
             "Ccl9","Cthrc1",#10
             "Lars2","Hspg2",#11
             "Spp1","Timp1","Acta2",#12
             "Atf3","Fosb")#13
ht <- GroupHeatmap(
  srt = sce,
  slot = "data",
  features = features,
  group.by = c("original_Cluster"),
  legend_title = c("Zscore(Average Expression)"),
  group_palette = c("simpsons"),
  reticle_color = "white",
  #heatmap_palette = "YlOrRd",
  heatmap_palcolor = c("navy", "white", "firebrick3"),
  cell_annotation = c("Cd248"),
  #  cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  show_column_names = TRUE,column_names_rot = 90,
  show_row_names = TRUE, row_names_side = "right",
  add_dot = T, add_reticle = F,
)
ht
ggsave("figures_Cd248/figures1_GroupHeatmap.pdf",width = 8,height = 8)
ggsave("figures_Cd248/figures1_GroupHeatmap.tiff",width = 8,height = 8)
ggsave("figures_Cd248/figures1_GroupHeatmap.svg",width = 8,height = 8)

# Heatmap of top marker genes for each fibroblast cluster ---------------------------------------------------------------
load("Rdata/marker_orig_cluster_heart_sham_3_28_Fib.Rdata")
top10 <- all.markers.pos %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top15 <- all.markers.pos %>% group_by(cluster) %>% top_n(20, avg_log2FC)
writexl::write_xlsx(top15,path = "results/top15gene.xlsx")
top10.F10 = top10$gene[top10$cluster == "F10"]
showfeatures = c("Mt1","Sat1",#1
                            "Smoc2","Gsn",#2
                            "Upk3b","Msln",#4
                            "Ccl2","Ptx3",#4
                            "Saa3","Serpinb2",#5
                            "Apoe","Igfbp3",#6
                            "Thbs4","Postn",#7
                            "Comp","Sfrp2",#8
                            "Cd248","Ackr3","Pi16",#9
                            "Ccl9","Cthrc1",#10
                            "Lars2","Hspg2",#11
                            "Spp1","Timp1","Acta2",#12
                            "Atf3","Fosb")#13
top10 = readxl::read_excel("results/top15gene.xlsx")
top10 = top10[!is.na(top10$select),]

ht <- GroupHeatmap(
  srt = sce,
  #slot = "data",
  features = top10$gene,
  group.by = c("original_Cluster"),
  legend_title = c("Zscore(Average Expression)"),
  group_palette = c("simpsons"),
  reticle_color = "white",
  #heatmap_palette = "YlOrRd",
  heatmap_palcolor = c("navy", "white", "firebrick3"),
  #  cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  show_column_names = TRUE,column_names_rot = 90,
  show_row_names = F, row_names_side = "right",
  add_dot = F, add_reticle = F,features_label = showfeatures,
)
ht
ggsave("figures_Cd248/figures1_GroupHeatmap_top10.pdf",width = 8,height = 8)
ggsave("figures_Cd248/figures1_GroupHeatmap_top10.tiff",width = 8,height = 8)
ggsave("figures_Cd248/figures1_GroupHeatmap_top10.svg",width = 8,height = 8)


#Gene set score for fibroblast cell states plotted in Umap embedding
DefaultAssay(sce) = "SCT"
sce =  AddModuleScore(object = sce,list(F1 = c("Mt1","Sat1","Tcf21")),name = "Fib1")
sce =  AddModuleScore(object = sce,list(F1 = c("Smoc2","Gsn","G0s2","Tcf21")),name = "Fib2")
sce =  AddModuleScore(object = sce,list(F1 = c("Upk3b","Msln","Igfbp5")),name = "Fib3")
sce =  AddModuleScore(object = sce,list(F1 = c("Ccl2","Ptx3","Cxcl2")),name = "Fib4")
sce =  AddModuleScore(object = sce,list(F1 = c("Saa3","Serpinb2","Cxcl13")),name = "Fib5")
sce =  AddModuleScore(object = sce,list(F1 = c("Apoe","Fmo2","Igfbp3")),name = "Fib6")
sce =  AddModuleScore(object = sce,list(F1 = c("Thbs4","Cilp","Postn","Aspn")),name = "Fib7")
sce =  AddModuleScore(object = sce,list(F1 = c("Comp","Sfrp2")),name = "Fib8")
sce =  AddModuleScore(object = sce,list(F1 = c("Ly6c1","Cd248","Pi16","Ackr3","Cd55")),name = "Fib9")
sce =  AddModuleScore(object = sce,list(F1 = c("Ccl9","Cthrc1","Nt5dc2")),name = "Fib10")
sce =  AddModuleScore(object = sce,list(F1 = c("Lars2","Hspg2")),name = "Fib11")
sce =  AddModuleScore(object = sce,list(F1 = c("Spp1","Timp1","Acta2","Tagln")),name = "Fib12")
sce =  AddModuleScore(object = sce,list(F1 = c("Atf3","Fosb","Hspa1a","Hspa1b")),name = "Fib13")
DPI = c(250,250)

F1 = FeatureDimPlot(sce, features = c("Fib11"), reduction = "umap", label = F,title = "F1:Mt1 Sat1 Tcf21",upper_quantile = 0.95,
               raster = T,raster.dpi = DPI,show_stat = F)
F2 = FeatureDimPlot(sce, features = c("Fib21"), reduction = "umap", label = F,title = "F2:Smoc2 Gsn Tcf21",upper_quantile = 0.95,
               raster = T,raster.dpi = DPI,show_stat = F)
F3 = FeatureDimPlot(sce, features = c("Fib31"), reduction = "umap", label = F,title = "F3:Upk3b Msln Igfbp5",upper_quantile = 0.95,
               raster = T,raster.dpi = DPI,show_stat = F)
F4 = FeatureDimPlot(sce, features = c("Fib41"), reduction = "umap", label = F,title = "F4:Ccl2 Ptx3 Cxcl2",upper_quantile = 0.95,
               raster = T,raster.dpi = DPI,show_stat = F)
F5 = FeatureDimPlot(sce, features = c("Fib51"), reduction = "umap", label = F,title = "F5:Saa3 Serpinb2 Cxcl13",upper_quantile = 0.95,
               raster = T,raster.dpi = DPI,show_stat = F)
F6 = FeatureDimPlot(sce, features = c("Fib61"), reduction = "umap", label = F,title = "F6:Apoe Fmo2 Igfbp3",upper_quantile = 0.95,
               raster = T,raster.dpi = DPI,show_stat = F)
F7 = FeatureDimPlot(sce, features = c("Fib71"), reduction = "umap", label = F,title = "F7:Thbs4 Cilp Postn Aspn",upper_quantile = 0.95,
               raster = T,raster.dpi = DPI,show_stat = F)
F8 = FeatureDimPlot(sce, features = c("Fib81"), reduction = "umap", label = F,title = "F8:Comp Sfrp2",upper_quantile = 0.95,
               raster = T,raster.dpi = DPI,show_stat = F)
F9 = FeatureDimPlot(sce, features = c("Fib91"), reduction = "umap", label = F,title = "F9:Cd248 Ackr3 Ly6c1 Pi16 Cd55",upper_quantile = 0.95,
               raster = T,raster.dpi = DPI,show_stat = F)
F10 = FeatureDimPlot(sce, features = c("Fib101"), reduction = "umap", label = F,title = "F10:Ccl9 Cthrc1 Nt5dc2",upper_quantile = 0.90,
               raster = T,raster.dpi = DPI,show_stat = F)
F11 = FeatureDimPlot(sce, features = c("Fib111"), reduction = "umap", label = F,title = "F11:Lars2 Hspg2",upper_quantile = 0.95,
               raster = T,raster.dpi = DPI,show_stat = F)
F12 = FeatureDimPlot(sce, features = c("Fib121"), reduction = "umap", label = F,title = "F12:Acta2 Tagln Spp1 Timp1",upper_quantile = 0.95,
               raster = T,raster.dpi = DPI,show_stat = F)
F13 = FeatureDimPlot(sce, features = c("Fib131"), reduction = "umap", label = F,title = "F13:Atf3 Fosb Hspa1a Hspa1b",upper_quantile = 0.95,
               raster = T,raster.dpi = DPI,show_stat = F)

patchwork::wrap_plots(list(F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12, F13), ncol = 4)
ggsave("figures_Cd248/figure1_FeatureDimPlot_F1-F13.pdf",width = 6*4,height = 4.8*4)

