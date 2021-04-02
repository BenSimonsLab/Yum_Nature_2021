##### Estimating the significance of transcriptome change for mesenchymal cells based on the separability using AUGUR 
##### Author: Seungmin Han (sh906@cam.ac.uk)
##### Last Update: 31/03/2021

library(scater)
library(scran)
library(Seurat)
library(Augur)
library(combinat)

result_file = '../data'
setwd(result_file)
load("R2_sce_Immune.RData")

colData(R2_sce_Immune)['condition'] <- vector(mode = 'character', length = dim(R2_sce_Immune)[2])
tmp_idx <- R2_sce_Immune$sample_name %in% c('CONF_2', 'CONF_3')
R2_sce_Immune$condition[tmp_idx] <- c('CONF')
tmp_idx <- R2_sce_Immune$sample_name %in% c('R2KR_1', 'R2KR_2', 'R2KR_3')
R2_sce_Immune$condition[tmp_idx] <- c('R2KR')
tmp_idx <- R2_sce_Immune$sample_name %in% c('R2P3_1', 'R2P3_2')
R2_sce_Immune$condition[tmp_idx] <- c('R2P3')

all_diffprior_KR <- list()
N_repeat <- 50
set.seed(500)
tmp_seed <- sample(c(1:10000),100)

for (i in c(1:N_repeat)) {
  tmp_idx <- which( R2_sce_Immune$condition %in% c('CONF') )
  set.seed(tmp_seed[i])
  tmp_idx2 <- sample(tmp_idx, round(length(tmp_idx) / 2))
  tmp_idx3 <- setdiff(tmp_idx, tmp_idx2)
  colData(R2_sce_Immune)$condition_pseudo <- colData(R2_sce_Immune)$condition
  colData(R2_sce_Immune)$condition_pseudo[tmp_idx2] <- 'CONF_C'
  colData(R2_sce_Immune)$condition_pseudo[tmp_idx3] <- 'CONF_T'
  
  tmp_idx <- R2_sce_Immune$condition %in% c('R2P3')
  R2KR_sceIm <- R2_sce_Immune[,!tmp_idx]
  
  tmp_idx <- R2KR_sceIm$condition_pseudo %in% c('CONF_C', 'CONF_T')
  R2KR_CT_sceIm <- R2KR_sceIm[,tmp_idx]
  tmp_idx2 <- R2KR_sceIm$condition_pseudo %in% c('CONF_C', 'R2KR')
  R2KR_KR_sceIm <- R2KR_sceIm[,tmp_idx2]

  # This needs to run Augur.
  colData(R2KR_CT_sceIm)$cell_type <- R2KR_CT_sceIm$celltype
  colData(R2KR_CT_sceIm)$label <- R2KR_CT_sceIm$condition_pseudo
  colData(R2KR_KR_sceIm)$cell_type <- R2KR_KR_sceIm$celltype
  colData(R2KR_KR_sceIm)$label <- R2KR_KR_sceIm$condition_pseudo
  
  R2KR_CT_Im_augur = calculate_auc(R2KR_CT_sceIm, subsample_size=6) 
  R2KR_CT_Im_augur_pm = calculate_auc(R2KR_CT_sceIm, subsample_size=6, augur_mode="permute")
  R2KR_KR_Im_augur = calculate_auc(R2KR_KR_sceIm, subsample_size=6) 
  R2KR_KR_Im_augur_pm = calculate_auc(R2KR_KR_sceIm, subsample_size=6, augur_mode="permute") 

  diff_prior_KR <- calculate_differential_prioritization(R2KR_CT_Im_augur, 
                                                         R2KR_KR_Im_augur, 
                                                         R2KR_CT_Im_augur_pm, 
                                                         R2KR_KR_Im_augur_pm)
  all_diffprior_KR[[i]] <- diff_prior_KR
  print( paste0("N = ",as.character(i)) )
}

all_KR_df <- all_diffprior_KR[[1]]
for (i in c(2:50)) {
  all_KR_df <- rbind(all_KR_df, all_diffprior_KR[[i]])
}

celltypes_KR_Im <- names( table(all_diffprior_KR[[1]]$cell_type) )
pv_augur_Im_KR_df <- data.frame(cell_type=character(),
                                condition=character(),
                                pval=double(),
                                stringsAsFactors=FALSE
)

for (celltype in celltypes_KR_Im) {
  tmp_idx <- all_KR_df$cell_type %in% celltype
  tmp_df <- all_KR_df[tmp_idx,]
  tmp_df2 <- data.frame(cell_type=celltype, condition='R2KR', pval=median(tmp_df$pval), stringsAsFactors = F )
  pv_augur_Im_KR_df <- rbind(pv_augur_Im_KR_df, tmp_df2)
}

all_diffprior_P3 <- list()
N_repeat <- 50
set.seed(500)
tmp_seed <- sample(c(1:10000),100)

for (i in c(1:N_repeat)) {
  tmp_idx <- which( R2_sce_Immune$condition %in% c('CONF') )
  set.seed(tmp_seed[i])
  tmp_idx2 <- sample(tmp_idx, round(length(tmp_idx) / 2))
  tmp_idx3 <- setdiff(tmp_idx, tmp_idx2)
  colData(R2_sce_Immune)$condition_pseudo <- colData(R2_sce_Immune)$condition
  colData(R2_sce_Immune)$condition_pseudo[tmp_idx2] <- 'CONF_C'
  colData(R2_sce_Immune)$condition_pseudo[tmp_idx3] <- 'CONF_T'
  
  tmp_idx <- R2_sce_Immune$condition %in% c('R2KR')
  R2P3_sceIm <- R2_sce_Immune[,!tmp_idx]
  
  tmp_idx <- R2P3_sceIm$condition_pseudo %in% c('CONF_C', 'CONF_T')
  R2P3_CT_sceIm <- R2P3_sceIm[,tmp_idx]
  tmp_idx2 <- R2P3_sceIm$condition_pseudo %in% c('CONF_C', 'R2P3')
  R2P3_P3_sceIm <- R2P3_sceIm[,tmp_idx2]

  # This needs to run Augur.
  colData(R2P3_CT_sceIm)$cell_type <- R2P3_CT_sceIm$celltype
  colData(R2P3_CT_sceIm)$label <- R2P3_CT_sceIm$condition_pseudo
  colData(R2P3_P3_sceIm)$cell_type <- R2P3_P3_sceIm$celltype
  colData(R2P3_P3_sceIm)$label <- R2P3_P3_sceIm$condition_pseudo

  R2P3_CT_Im_augur = calculate_auc(R2P3_CT_sceIm, subsample_size=6) 
  R2P3_CT_Im_augur_pm = calculate_auc(R2P3_CT_sceIm, subsample_size=6, augur_mode="permute")
  R2P3_P3_Im_augur = calculate_auc(R2P3_P3_sceIm, subsample_size=6) 
  R2P3_P3_Im_augur_pm = calculate_auc(R2P3_P3_sceIm, subsample_size=6, augur_mode="permute") 

  diff_prior_P3 <- calculate_differential_prioritization(R2P3_CT_Im_augur, 
                                                         R2P3_P3_Im_augur, 
                                                         R2P3_CT_Im_augur_pm, 
                                                         R2P3_P3_Im_augur_pm)
  all_diffprior_P3[[i]] <- diff_prior_P3
  print( paste0("N = ",as.character(i)) )
}

all_P3_df <- all_diffprior_P3[[1]]
for (i in c(2:50)) {
  all_P3_df <- rbind(all_P3_df, all_diffprior_P3[[i]])
}

celltypes_P3_Im <- names( table(all_diffprior_P3[[1]]$cell_type) )
pv_augur_Im_P3_df <- data.frame(cell_type=character(),
                                condition=character(),
                                pval=double(),
                                stringsAsFactors=FALSE
)

for (celltype in celltypes_P3_Im) {
  tmp_idx <- all_P3_df$cell_type %in% celltype
  tmp_df <- all_P3_df[tmp_idx,]
  tmp_df2 <- data.frame(cell_type=celltype, condition='R2P3', pval=median(tmp_df$pval), stringsAsFactors = F )
  pv_augur_Im_P3_df <- rbind(pv_augur_Im_P3_df, tmp_df2)
}

pv_augur_Im_df <- rbind(pv_augur_Im_KR_df, pv_augur_Im_P3_df)
pv_augur_Im_df # This is p-value for the significance of transcriptomic change of immune cells between Confetti and Red2Onco models based on Augur method.

