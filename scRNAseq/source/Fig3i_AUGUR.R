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
load("R2_sce_Mesen.RData")

colData(R2_sce_Mesen)['condition'] <- vector(mode = 'character', length = dim(R2_sce_Mesen)[2])
tmp_idx <- R2_sce_Mesen$sample_name %in% c('CONF_2', 'CONF_3')
R2_sce_Mesen$condition[tmp_idx] <- c('CONF')
tmp_idx <- R2_sce_Mesen$sample_name %in% c('R2KR_1', 'R2KR_2', 'R2KR_3')
R2_sce_Mesen$condition[tmp_idx] <- c('R2KR')
tmp_idx <- R2_sce_Mesen$sample_name %in% c('R2P3_1', 'R2P3_2')
R2_sce_Mesen$condition[tmp_idx] <- c('R2P3')

all_diffprior_KR <- list()
N_repeat <- 50
set.seed(500)
tmp_seed <- sample(c(1:10000),100)

for (i in c(1:N_repeat)) {
  tmp_idx <- which( R2_sce_Mesen$condition %in% c('CONF') )
  set.seed(tmp_seed[i])
  tmp_idx2 <- sample(tmp_idx, round(length(tmp_idx) / 2))
  tmp_idx3 <- setdiff(tmp_idx, tmp_idx2)
  colData(R2_sce_Mesen)$condition_pseudo <- colData(R2_sce_Mesen)$condition
  colData(R2_sce_Mesen)$condition_pseudo[tmp_idx2] <- 'CONF_C'
  colData(R2_sce_Mesen)$condition_pseudo[tmp_idx3] <- 'CONF_T'
  
  # Filtering out intestitial cells of Cajal because their cell number in Confetti is too small.
  tmp_idx <- R2_sce_Mesen$celltype %in% c('Intestitial_cell')
  R2_sce_Mesen2 <- R2_sce_Mesen[,!tmp_idx]
  tmp_idx <- R2_sce_Mesen2$condition %in% c('R2P3')
  R2KR_sceMs <- R2_sce_Mesen2[,!tmp_idx]
  
  tmp_idx <- R2KR_sceMs$condition_pseudo %in% c('CONF_C', 'CONF_T')
  R2KR_CT_sceMs <- R2KR_sceMs[,tmp_idx]
  tmp_idx2 <- R2KR_sceMs$condition_pseudo %in% c('CONF_C', 'R2KR')
  R2KR_KR_sceMs <- R2KR_sceMs[,tmp_idx2]
  
  # This needs to run Augur.
  colData(R2KR_CT_sceMs)$cell_type <- R2KR_CT_sceMs$celltype
  colData(R2KR_CT_sceMs)$label <- R2KR_CT_sceMs$condition_pseudo
  colData(R2KR_KR_sceMs)$cell_type <- R2KR_KR_sceMs$celltype
  colData(R2KR_KR_sceMs)$label <- R2KR_KR_sceMs$condition_pseudo
  
  R2KR_CT_Ms_augur = calculate_auc(R2KR_CT_sceMs, subsample_size=6) 
  R2KR_CT_Ms_augur_pm = calculate_auc(R2KR_CT_sceMs, subsample_size=6, augur_mode="permute")
  R2KR_KR_Ms_augur = calculate_auc(R2KR_KR_sceMs, subsample_size=6) 
  R2KR_KR_Ms_augur_pm = calculate_auc(R2KR_KR_sceMs, subsample_size=6, augur_mode="permute") 
  
  diff_prior_KR <- calculate_differential_prioritization(R2KR_CT_Ms_augur, 
                                                         R2KR_KR_Ms_augur, 
                                                         R2KR_CT_Ms_augur_pm, 
                                                         R2KR_KR_Ms_augur_pm)
  all_diffprior_KR[[i]] <- diff_prior_KR
  print( paste0("N = ",as.character(i)) )
}

all_KR_df <- all_diffprior_KR[[1]]
for (i in c(2:50)) {
  all_KR_df <- rbind(all_KR_df, all_diffprior_KR[[i]])
}

celltypes_KR_Ms <- names( table(all_diffprior_KR[[1]]$cell_type) )
pv_augur_Ms_KR_df <- data.frame(cell_type=character(),
                                condition=character(),
                                pval=double(),
                                stringsAsFactors=FALSE
)

for (celltype in celltypes_KR_Ms) {
  tmp_idx <- all_KR_df$cell_type %in% celltype
  tmp_df <- all_KR_df[tmp_idx,]
  tmp_df2 <- data.frame(cell_type=celltype, condition='R2KR', pval=median(tmp_df$pval), stringsAsFactors = F )
  pv_augur_Ms_KR_df <- rbind(pv_augur_Ms_KR_df, tmp_df2)
}


all_diffprior_P3 <- list()
N_repeat <- 50
set.seed(500)
tmp_seed <- sample(c(1:10000),100)

for (i in c(1:N_repeat)) {
  tmp_idx <- which( R2_sce_Mesen$condition %in% c('CONF') )
  set.seed(tmp_seed[i])
  tmp_idx2 <- sample(tmp_idx, round(length(tmp_idx) / 2))
  tmp_idx3 <- setdiff(tmp_idx, tmp_idx2)
  colData(R2_sce_Mesen)$condition_pseudo <- colData(R2_sce_Mesen)$condition
  colData(R2_sce_Mesen)$condition_pseudo[tmp_idx2] <- 'CONF_C'
  colData(R2_sce_Mesen)$condition_pseudo[tmp_idx3] <- 'CONF_T'

  # Filtering out intestitial cells and myofibroblast cells because their cell numbers are too small.
  tmp_idx <- R2_sce_Mesen$celltype %in% c('Intestitial_cell', 'Myofibroblast')
  R2_sce_Mesen2 <- R2_sce_Mesen[,!tmp_idx]
  tmp_idx <- R2_sce_Mesen2$condition %in% c('R2KR')
  R2P3_sceMs <- R2_sce_Mesen2[,!tmp_idx]

  tmp_idx <- R2P3_sceMs$condition_pseudo %in% c('CONF_C', 'CONF_T')
  R2P3_CT_sceMs <- R2P3_sceMs[,tmp_idx]
  tmp_idx2 <- R2P3_sceMs$condition_pseudo %in% c('CONF_C', 'R2P3')
  R2P3_P3_sceMs <- R2P3_sceMs[,tmp_idx2]

  # This needs to run Augur.
  colData(R2P3_CT_sceMs)$cell_type <- R2P3_CT_sceMs$celltype
  colData(R2P3_CT_sceMs)$label <- R2P3_CT_sceMs$condition_pseudo
  colData(R2P3_P3_sceMs)$cell_type <- R2P3_P3_sceMs$celltype
  colData(R2P3_P3_sceMs)$label <- R2P3_P3_sceMs$condition_pseudo
  
  R2P3_CT_Ms_augur = calculate_auc(R2P3_CT_sceMs, subsample_size=6) 
  R2P3_CT_Ms_augur_pm = calculate_auc(R2P3_CT_sceMs, subsample_size=6, augur_mode="permute")
  R2P3_P3_Ms_augur = calculate_auc(R2P3_P3_sceMs, subsample_size=6) 
  R2P3_P3_Ms_augur_pm = calculate_auc(R2P3_P3_sceMs, subsample_size=6, augur_mode="permute") 

  diff_prior_P3 <- calculate_differential_prioritization(R2P3_CT_Ms_augur, 
                                                         R2P3_P3_Ms_augur, 
                                                         R2P3_CT_Ms_augur_pm, 
                                                         R2P3_P3_Ms_augur_pm)
  all_diffprior_P3[[i]] <- diff_prior_P3
  print( paste0("N = ",as.character(i)) )

}

all_P3_df <- all_diffprior_P3[[1]]
for (i in c(2:50)) {
  all_P3_df <- rbind(all_P3_df, all_diffprior_P3[[i]])
}

celltypes_P3_Ms <- names( table(all_diffprior_P3[[1]]$cell_type) )
pv_augur_Ms_P3_df <- data.frame(cell_type=character(),
                                condition=character(),
                                pval=double(),
                                stringsAsFactors=FALSE
)

for (celltype in celltypes_P3_Ms) {
  tmp_idx <- all_P3_df$cell_type %in% celltype
  tmp_df <- all_P3_df[tmp_idx,]
  tmp_df2 <- data.frame(cell_type=celltype, condition='R2P3', pval=median(tmp_df$pval), stringsAsFactors = F )
  pv_augur_Ms_P3_df <- rbind(pv_augur_Ms_P3_df, tmp_df2)
}

pv_augur_Ms_df <- rbind(pv_augur_Ms_KR_df, pv_augur_Ms_P3_df)
pv_augur_Ms_df # This is p-value for the significance of transcriptomic change of mesenchymal cells between Confetti and Red2Onco models based on AUGUR method.




