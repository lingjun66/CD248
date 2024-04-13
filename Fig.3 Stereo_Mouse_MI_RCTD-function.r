###shaohui shi
###shaohuishi98@gmail.com
##RCTD for stereo-seq bin50 deconvolution

library(Seurat)
library(spacexr)
anno <- readRDS('~/result/process/anno/heart_Annotation.rds')

Idents(anno) <- 'Main_Cell_type'
sub_anno <- subset(anno,downsample=400)

sc_counts <- sub_anno[['RNA']]@counts
str(sc_counts)

cell_types <- Idents(sub_anno)
str(cell_types)

function(day){
    stereo <- readRDS(paste0('~/rawdata/stereo/after_qc/',day,
                            'h_bin50.rds'))

    head(stereo)
    coords <- stereo@meta.data[c('x','y')]
    rownames(coords) <- rownames(stereo@meta.data)
    sp_counts <- stereo[['Spatial']]@counts
    #str(sp_counts)

    reference <- Reference(sc_counts,cell_types)
    spaceRNA <- SpatialRNA(coords,sp_counts)

    myRCTD <- create.RCTD(spaceRNA,reference,max_cores = 64, test_mode = FALSE)
    myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full') 

    saveRDS(myRCTD,paste0('~/result/stereo/data/',day,
                        '/RCTD_full_cd248.rds'))
}