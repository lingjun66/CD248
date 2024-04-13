###shaohui shi
###shaohuishi98@gmail.com
##Nebulosa for gene density

library(Seurat)
library("Nebulosa")
library(ggplot2)
library(ggnewscale)

sample <- c('sham','m14h','m28h')
lapply(sample,function(day){
stereo <- readRDS(paste0('/data/slurm/shish/data/HuLab/Visium/stereo/after_qc/',day,'_heart_bin50.rds'))
coor = data.frame('UMAP_1'= stereo@images$slice1@coordinates$imagerow , 'UMAP_2' = stereo@images$slice1@coordinates$imagecol)
rownames(coor) = colnames(stereo)
coor = as.matrix(coor)
stereo@reductions$umaps <- stereo@reductions$umap
stereo@reductions$umaps@cell.embeddings = coor  
stereo@reductions$umaps@key <- 'umaps'

p1 = plot_density(stereo, c('Cd3e','Cd248'), joint = TRUE) &
     scale_color_gradient2(limits = c(0, 9e-08), low = "#440154", mid = '#299b87',  
                        high = "#fce726",midpoint = 4.5e-08)


p = p1 + scale_color_gradient2(limits = c(0, 8e-16), low = "#440154", mid = '#299b87',  
                        high = "#fce726",midpoint = 4e-16)


ggsave(paste0('/data/slurm/shish/data/HuLab/Visium/stereo/output/',day,'_cd3_7_density.pdf'),p,width = 30,height = 19)
})