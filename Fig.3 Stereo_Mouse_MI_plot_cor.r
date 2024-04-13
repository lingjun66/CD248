### shaohui shi
### shaohuishi98@gmail.com
##stereo-seq GEF QC&gaussian smooth&convert format
library(CARD)
library(Seurat)
library(tidyverse)

plot_card_inf_data <- function(day){
    myRCTD <- readRDS(paste0('~/result/stereo/data/',day,
                        '/RCTD_full_cd248.rds'))
    RCTD_results_full <- myRCTD@results
    RCTD_norm_weights = as.matrix(sweep(RCTD_results_full$weights, 1, 
                        rowSums(RCTD_results_full$weights), '/'))
    RCTD_norm_weights_tibb = as_tibble(RCTD_norm_weights, rownames = "barcodes")
    res_CARD <- as.data.frame(RCTD_norm_weights_tibb)
    rownames(res_CARD) <- res_CARD$barcodes
    res_CARD <- res_CARD[,-grep("barcodes",colnames(res_CARD))]
    res_CARD = res_CARD[,mixedsort(colnames(res_CARD))]
    return(res_CARD)
    }

plot_card_inf_loc <- function(day){
    stereo <- readRDS(paste0('~/rawdata/stereo/after_qc/',day,
                            'h_bin50.rds'))
    coords <- stereo@meta.data[c('x','y')]
    rownames(coords) <- rownames(stereo@meta.data)
    location = as.data.frame(coords)
    if(nrow(res_CARD)!=nrow(location)){ 
        if(nrow(res_CARD)>nrow(location)){
            res_CARD <- res_CARD[rownames(location),]
        }else{
            location <- location[rownames(res_CARD),]}}
    return(location)
}


day = 'm1h'
res_CARD <- plot_card_inf_data(day)
location <- plot_card_inf_loc(day)
sub_location <- filter(location,
        location$x < 14000 & location$y < 11000)
#res_CARD <- res_CARD[rownames(location),]
res_CARD2 <- res_CARD[rownames(sub_location),]
p1 <- CARD.visualize.Cor(res_CARD1,colors = NULL)
p2 <- CARD.visualize.Cor(res_CARD2,colors = NULL)
proportion = res_CARD[,order(colnames(res_CARD))]
proportion2 = res_CARD2[,order(colnames(res_CARD2))]
cor_CARD = cor.test(proportion$`Cd248 Fibroblasts`,proportion$`T cells`,conf.level = 0.95)
d95 <- cor_CARD$conf.int[c(1,2)]
d95 <- t(as.data.frame(d95))
d0 <- cor_CARD$estimate[1]
d0 <- t(as.data.frame(d0))
df <- cbind(d0,d95)
colnames(df) <- c('cor','ymax','ymin')
df <- as.data.frame(df)
df$celltype <- 'Cd248 Fibroblasts'
df$inf <- 'all'
df$sample <- paste0(day)

cor_CARD = cor.test(proportion$`Remaining Fibroblasts`,proportion$`T cells`,conf.level = 0.95)
d95 <- cor_CARD$conf.int[c(1,2)]
d95 <- t(as.data.frame(d95))
d0 <- cor_CARD$estimate[1]
d0 <- t(as.data.frame(d0))
df2 <- cbind(d0,d95)
colnames(df2) <- c('cor','ymax','ymin')
df2 <- as.data.frame(df2)
df2$celltype <- 'Ramaining Fibroblasts'
df2$inf <- 'all'
df2$sample <- paste0(day)
df_pic <- rbind(df,df2)

cor_CARD = cor.test(proportion2$`Cd248 Fibroblasts`,proportion2$`T cells`,conf.level = 0.95)
d95 <- cor_CARD$conf.int[c(1,2)]
d95 <- t(as.data.frame(d95))
d0 <- cor_CARD$estimate[1]
d0 <- t(as.data.frame(d0))
df <- cbind(d0,d95)
colnames(df) <- c('cor','ymax','ymin')
df <- as.data.frame(df)
df$celltype <- 'Cd248 Fibroblasts'
df$inf <- 'inf'
df$sample <- paste0(day)

cor_CARD = cor.test(proportion2$`Remaining Fibroblasts`,proportion2$`T cells`,conf.level = 0.95)
d95 <- cor_CARD$conf.int[c(1,2)]
d95 <- t(as.data.frame(d95))
d0 <- cor_CARD$estimate[1]
d0 <- t(as.data.frame(d0))
df2 <- cbind(d0,d95)
colnames(df2) <- c('cor','ymax','ymin')
df2 <- as.data.frame(df2)
df2$celltype <- 'Ramaining Fibroblasts'
df2$inf <- 'inf'
df2$sample <- paste0(day)
df_pic <- rbind(df_pic,df,df2)

day = 'm14h'
res_CARD <- plot_card_inf_data(day)
location <- plot_card_inf_loc(day)
sub_location <- filter(location,
        location$y > 18000 | location$x >15500)
#res_CARD <- res_CARD[rownames(location),]
res_CARD2 <- res_CARD[rownames(sub_location),]
p3 <- CARD.visualize.Cor(res_CARD1,colors = NULL)
p4 <- CARD.visualize.Cor(res_CARD2,colors = NULL)

proportion = res_CARD[,order(colnames(res_CARD))]
proportion2 = res_CARD2[,order(colnames(res_CARD2))]
cor_CARD = cor.test(proportion$`Cd248 Fibroblasts`,proportion$`T cells`,conf.level = 0.95)
d95 <- cor_CARD$conf.int[c(1,2)]
d95 <- t(as.data.frame(d95))
d0 <- cor_CARD$estimate[1]
d0 <- t(as.data.frame(d0))
df <- cbind(d0,d95)
colnames(df) <- c('cor','ymax','ymin')
df <- as.data.frame(df)
df$celltype <- 'Cd248 Fibroblasts'
df$inf <- 'all'
df$sample <- paste0(day)

cor_CARD = cor.test(proportion$`Remaining Fibroblasts`,proportion$`T cells`,conf.level = 0.95)
d95 <- cor_CARD$conf.int[c(1,2)]
d95 <- t(as.data.frame(d95))
d0 <- cor_CARD$estimate[1]
d0 <- t(as.data.frame(d0))
df2 <- cbind(d0,d95)
colnames(df2) <- c('cor','ymax','ymin')
df2 <- as.data.frame(df2)
df2$celltype <- 'Ramaining Fibroblasts'
df2$inf <- 'all'
df2$sample <- paste0(day)
df_pic <- rbind(df_pic,df,df2)

cor_CARD = cor.test(proportion2$`Cd248 Fibroblasts`,proportion2$`T cells`,conf.level = 0.95)
d95 <- cor_CARD$conf.int[c(1,2)]
d95 <- t(as.data.frame(d95))
d0 <- cor_CARD$estimate[1]
d0 <- t(as.data.frame(d0))
df <- cbind(d0,d95)
colnames(df) <- c('cor','ymax','ymin')
df <- as.data.frame(df)
df$celltype <- 'Cd248 Fibroblasts'
df$inf <- 'inf'
df$sample <- paste0(day)

cor_CARD = cor.test(proportion2$`Remaining Fibroblasts`,proportion2$`T cells`,conf.level = 0.95)
d95 <- cor_CARD$conf.int[c(1,2)]
d95 <- t(as.data.frame(d95))
d0 <- cor_CARD$estimate[1]
d0 <- t(as.data.frame(d0))
df2 <- cbind(d0,d95)
colnames(df2) <- c('cor','ymax','ymin')
df2 <- as.data.frame(df2)
df2$celltype <- 'Ramaining Fibroblasts'
df2$inf <- 'inf'
df2$sample <- paste0(day)
df_pic <- rbind(df_pic,df,df2)

day = 'm28h'
res_CARD <- plot_card_inf_data(day)
location <- plot_card_inf_loc(day)
sub_location <- filter(location,
        location$x < 10000)
sub_location <- filter(sub_location,
        sub_location$y > 14000 | sub_location$x < 8800)   
#res_CARD <- res_CARD[rownames(location),]
res_CARD2 <- res_CARD[rownames(sub_location),]

p5 <- CARD.visualize.Cor(res_CARD1,colors = NULL)
p6 <- CARD.visualize.Cor(res_CARD2,colors = NULL)

proportion = res_CARD[,order(colnames(res_CARD))]
proportion2 = res_CARD2[,order(colnames(res_CARD2))]
cor_CARD = cor.test(proportion$`Cd248 Fibroblasts`,proportion$`T cells`,conf.level = 0.95)
d95 <- cor_CARD$conf.int[c(1,2)]
d95 <- t(as.data.frame(d95))
d0 <- cor_CARD$estimate[1]
d0 <- t(as.data.frame(d0))
df <- cbind(d0,d95)
colnames(df) <- c('cor','ymax','ymin')
df <- as.data.frame(df)
df$celltype <- 'Cd248 Fibroblasts'
df$inf <- 'all'
df$sample <- paste0(day)

cor_CARD = cor.test(proportion$`Remaining Fibroblasts`,proportion$`T cells`,conf.level = 0.95)
d95 <- cor_CARD$conf.int[c(1,2)]
d95 <- t(as.data.frame(d95))
d0 <- cor_CARD$estimate[1]
d0 <- t(as.data.frame(d0))
df2 <- cbind(d0,d95)
colnames(df2) <- c('cor','ymax','ymin')
df2 <- as.data.frame(df2)
df2$celltype <- 'Ramaining Fibroblasts'
df2$inf <- 'all'
df2$sample <- paste0(day)
df_pic <- rbind(df_pic,df,df2)

cor_CARD = cor.test(proportion2$`Cd248 Fibroblasts`,proportion2$`T cells`,conf.level = 0.95)
d95 <- cor_CARD$conf.int[c(1,2)]
d95 <- t(as.data.frame(d95))
d0 <- cor_CARD$estimate[1]
d0 <- t(as.data.frame(d0))
df <- cbind(d0,d95)
colnames(df) <- c('cor','ymax','ymin')
df <- as.data.frame(df)
df$celltype <- 'Cd248 Fibroblasts'
df$inf <- 'inf'
df$sample <- paste0(day)

cor_CARD = cor.test(proportion2$`Remaining Fibroblasts`,proportion2$`T cells`,conf.level = 0.95)
d95 <- cor_CARD$conf.int[c(1,2)]
d95 <- t(as.data.frame(d95))
d0 <- cor_CARD$estimate[1]
d0 <- t(as.data.frame(d0))
df2 <- cbind(d0,d95)
colnames(df2) <- c('cor','ymax','ymin')
df2 <- as.data.frame(df2)
df2$celltype <- 'Ramaining Fibroblasts'
df2$inf <- 'inf'
df2$sample <- paste0(day)
df_pic <- rbind(df_pic,df,df2)

day = 'sham'
res_CARD <- plot_card_inf_data(day)
location <- plot_card_inf_loc(day)
sub_location <- filter(location,
        location$y-location$x > 5000)
#res_CARD <- res_CARD[rownames(location),]
res_CARD2 <- res_CARD[rownames(sub_location),]

proportion = res_CARD[,order(colnames(res_CARD))]
proportion2 = res_CARD2[,order(colnames(res_CARD2))]
cor_CARD = cor.test(proportion$`Cd248 Fibroblasts`,proportion$`T cells`,conf.level = 0.95)
d95 <- cor_CARD$conf.int[c(1,2)]
d95 <- t(as.data.frame(d95))
d0 <- cor_CARD$estimate[1]
d0 <- t(as.data.frame(d0))
df <- cbind(d0,d95)
colnames(df) <- c('cor','ymax','ymin')
df <- as.data.frame(df)
df$celltype <- 'Cd248 Fibroblasts'
df$inf <- 'all'
df$sample <- paste0(day)

cor_CARD = cor.test(proportion$`Remaining Fibroblasts`,proportion$`T cells`,conf.level = 0.95)
d95 <- cor_CARD$conf.int[c(1,2)]
d95 <- t(as.data.frame(d95))
d0 <- cor_CARD$estimate[1]
d0 <- t(as.data.frame(d0))
df2 <- cbind(d0,d95)
colnames(df2) <- c('cor','ymax','ymin')
df2 <- as.data.frame(df2)
df2$celltype <- 'Ramaining Fibroblasts'
df2$inf <- 'all'
df2$sample <- paste0(day)
df_pic <- rbind(df_pic,df,df2)

cor_CARD = cor.test(proportion2$`Cd248 Fibroblasts`,proportion2$`T cells`,conf.level = 0.95)
d95 <- cor_CARD$conf.int[c(1,2)]
d95 <- t(as.data.frame(d95))
d0 <- cor_CARD$estimate[1]
d0 <- t(as.data.frame(d0))
df <- cbind(d0,d95)
colnames(df) <- c('cor','ymax','ymin')
df <- as.data.frame(df)
df$celltype <- 'Cd248 Fibroblasts'
df$inf <- 'inf'
df$sample <- paste0(day)

cor_CARD = cor.test(proportion2$`Remaining Fibroblasts`,proportion2$`T cells`,conf.level = 0.95)
d95 <- cor_CARD$conf.int[c(1,2)]
d95 <- t(as.data.frame(d95))
d0 <- cor_CARD$estimate[1]
d0 <- t(as.data.frame(d0))
df2 <- cbind(d0,d95)
colnames(df2) <- c('cor','ymax','ymin')
df2 <- as.data.frame(df2)
df2$celltype <- 'Ramaining Fibroblasts'
df2$inf <- 'inf'
df2$sample <- paste0(day)
df_pic <- rbind(df_pic,df,df2)

table(df_pic$sample)
df_pic$sample <- factor(df_pic$sample,levels = c("sham", "m1h", "m3h", "m7h", "m14h", "m28h"))

p = ggplot(df_pic,aes(x=sample,y=cor,fill=celltype)) + 
    geom_bar(stat="identity",position = 'dodge',
            width = 0.8,
            color='black') +
    geom_errorbar(data = df_pic, aes(x = sample,  ymin = ymin, ymax = ymax,
    group = celltype),position = position_dodge(0.8),
    width= 0.2) +
    scale_fill_manual(values=c('#f2b77c', '#099963')) +
  facet_wrap(~inf)
pa <- p1 / p2 +p3 / p4 + p5 / p6
ggsave('/result/stereo/pic/cor_sum_bar.pdf',p)
ggsave('/result/stereo/pic/cor_sum.pdf',pa)