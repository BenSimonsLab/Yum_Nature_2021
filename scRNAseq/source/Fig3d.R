##### Comparing priming scores of WT epithelial cells for Confetti and Red2Onco models toward enterocyte lineage.
##### Author: Seungmin Han (sh906@cam.ac.uk)
##### Last Update: 31/03/2021

library(scater)
library(scran)
library(Seurat)
library(FateID)
library(dplyr)

result_path <-"../data"
setwd(result_file)
load("R2_sce_epi.RData")
load("Geneset_Sec.RData") # Geneset for Secretory lineages
load("Geneset_Ent.RData") # Geneset for Enterocyte lineage

colData(R2_sce_epi)['celltype_int'] <- vector(mode='integer', length=dim(R2_sce_epi)[2])

tmp_idx <- R2_sce_epi$celltype %in% c('F_ISC')
R2_sce_epi$celltype_int[tmp_idx] <- c(1)
tmp_idx <- R2_sce_epi$celltype %in% c('G_TA')
R2_sce_epi$celltype_int[tmp_idx] <- c(2)
tmp_idx <- R2_sce_epi$celltype %in% c('B_EP')
R2_sce_epi$celltype_int[tmp_idx] <- c(3)
tmp_idx <- R2_sce_epi$celltype %in% c('A_Ent')
R2_sce_epi$celltype_int[tmp_idx] <- c(4)
tmp_idx <- R2_sce_epi$celltype %in% c('E_Paneth')
R2_sce_epi$celltype_int[tmp_idx] <- c(5)
tmp_idx <- R2_sce_epi$celltype %in% c('D_Goblet')
R2_sce_epi$celltype_int[tmp_idx] <- c(6)
tmp_idx <- R2_sce_epi$celltype %in% c('H_Tuft')
R2_sce_epi$celltype_int[tmp_idx] <- c(7)
tmp_idx <- R2_sce_epi$celltype %in% c('C_EEC')
R2_sce_epi$celltype_int[tmp_idx] <- c(8)

tmp_idx <- R2_sce_epi$cond_color %in% c("CONF_R", "CONF_Y", "R2KR_Y",  "R2P3_Y")
R2WT_sce_epi <- R2_sce_epi[,tmp_idx]
v <- as.matrix( assays(R2WT_sce_epi)$logcounts )
tmp_idx <- rownames(v) %in% union(Geneset_Sec, Geneset_Ent)
v <- v[tmp_idx,]
y <- R2WT_sce_epi$celltype_int
tar <- c(4, 5, 6, 7, 8)
rc <- reclassify(v, y, tar, clthr=.75, 
                 nbfactor=5, use.dist=FALSE, seed=12345, nbtree=NULL, q=0.9) 
y <- rc$part
x <- v
fb <- fateBias(x, y, tar, z=NULL, minnr=5, minnrh=10, adapt=TRUE, 
               confidence=0.75, nbfactor=5, use.dist=FALSE, seed=12345, nbtree=NULL)
dr <- compdr(x, z=NULL, m=c("tsne","umap"), k=c(2,3), tsne.perplexity=30, seed=12345)
fb_short <- data.frame(fb_t4=fb$probs$t4, fb_t5=fb$probs$t5, fb_t6=fb$probs$t6, fb_t7=fb$probs$t7,
                       fb_t8=fb$probs$t8, stringsAsFactors=FALSE)

colData(R2_sce_epi)['cond_color2'] <- R2_sce_epi$cond_color
tmp_idx <- R2_sce_epi$cond_color %in% c('CONF_R', 'CONF_Y')
R2_sce_epi$cond_color2[tmp_idx] <- c('CONF')
tmp_idx <- R2_sce_epi$cond_color %in% c('CONF_R', 'CONF_Y', 'R2KR_Y', "R2P3_Y")
R2WT_sce_epi <- R2_sce_epi[,tmp_idx]

colData(R2WT_sce_epi)['FB_Ent'] <- vector(mode='numeric', length=dim(R2WT_sce_epi)[2])
R2WT_sce_epi$FB_Ent <- fb_short$fb_t4
colData(R2WT_sce_epi)['FB_Paneth'] <- vector(mode='numeric', length=dim(R2WT_sce_epi)[2])
R2WT_sce_epi$FB_Paneth <- fb_short$fb_t5
colData(R2WT_sce_epi)['FB_Goblet'] <- vector(mode='numeric', length=dim(R2WT_sce_epi)[2])
R2WT_sce_epi$FB_Goblet <- fb_short$fb_t6
colData(R2WT_sce_epi)['FB_Tuft'] <- vector(mode='numeric', length=dim(R2WT_sce_epi)[2])
R2WT_sce_epi$FB_Tuft <- fb_short$fb_t7
colData(R2WT_sce_epi)['FB_EEC'] <- vector(mode='numeric', length=dim(R2WT_sce_epi)[2])
R2WT_sce_epi$FB_EEC <- fb_short$fb_t8

tmp_idx <- R2WT_sce_epi$celltype %in% c('F_ISC', 'G_TA')
R2WT_ScTA_sce <- R2WT_sce_epi[,tmp_idx]

FB_df <- data.frame(FB_score=R2WT_ScTA_sce$FB_Ent, cond_color=R2WT_ScTA_sce$cond_color2, stringsAsFactors = FALSE)
p <- ggplot(FB_df, aes(x=cond_color, y=FB_score, color=cond_color)) + 
  geom_boxplot() +
  coord_flip() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("Priming score of WT epithelial cells for Confetti and Red2Onco models toward enterocyte lineage")


