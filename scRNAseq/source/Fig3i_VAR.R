##### Estimating the significance of transcriptome change for mesenchymal cells based on cell-cell variability 
##### Author: Seungmin Han (sh906@cam.ac.uk)
##### Last Update: 31/03/2021

library(scater)
library(scran)
library(cowplot)
library(ggplot2)
library(Seurat)
library(stats)
library(plyr)
library(RColorBrewer)

result_file = "../data"
setwd(result_file)
load("R2_sce_Mesen.RData")

colData(R2_sce_Mesen)["sample_name2"] <- vector(mode = "character", length= dim(R2_sce_Mesen)[2])
tmp_idx <- ( colData(R2_sce_Mesen)$sample_name %in% c("CONF_2", "CONF_3") )
colData(R2_sce_Mesen)$sample_name2[tmp_idx] <- "CONF"
tmp_idx <- ( colData(R2_sce_Mesen)$sample_name %in% c("R2KR_1", "R2KR_2", "R2KR_3") )
colData(R2_sce_Mesen)$sample_name2[tmp_idx] <- "R2KR"
tmp_idx <- ( colData(R2_sce_Mesen)$sample_name %in% c("R2P3_1", "R2P3_2") )
colData(R2_sce_Mesen)$sample_name2[tmp_idx] <- "R2P3"

### Module 1
cmp_tstat_SH <- function (tmp_sce, celltype_ofinterest) {
  
  tmp_idx <- colData(tmp_sce)$celltype == celltype_ofinterest
  tmp_sce <- tmp_sce[,tmp_idx]
  
  stat_df <- data.frame(
    celltype = character(),
    med_inter_var = double(),
    med_intra_var = double(),
    t_stat = double(),
    condition = character()
  )
  
  for (i in c(1:2)) {
    compare_list <- list( c("CONF", "R2KR"), c("CONF", "R2P3") )
    compare_of_interest <- compare_list[[i]]
    
    tmp_idx <- colData(tmp_sce)$sample_name2 %in% compare_of_interest
    tmp_sce_1 <- tmp_sce[,tmp_idx]
    
    tmp_idx <- ( Matrix::rowSums(assays(tmp_sce_1)$counts) > 0 )
    tmp_sce_1 <- tmp_sce_1[tmp_idx,]
    
    tmp_sce_1 <- scater::arrange(tmp_sce_1, sample_name2)
    
    tmp_cor_1 <- cor(as.matrix( assays(tmp_sce_1)$logcounts ))
    
    N_CONF <- sum( colData(tmp_sce_1)$sample_name2 == "CONF" )
    N_Treat <- sum( colData(tmp_sce_1)$sample_name2 == compare_of_interest[2] )
    
    intra_var_CONF <- tmp_cor_1[(1:N_CONF), (1:N_CONF)]
    intra_var_Treat <- tmp_cor_1[(N_CONF + 1):(N_CONF + N_Treat), (N_CONF + 1):(N_CONF + N_Treat)]
    intra_var_CONF <- intra_var_CONF[lower.tri(intra_var_CONF)]
    intra_var_Treat <- intra_var_Treat[lower.tri(intra_var_Treat)]
    intra_var <- (1 - c(intra_var_CONF, intra_var_Treat))
    rm(intra_var_CONF)
    rm(intra_var_Treat)
    inter_var <- tmp_cor_1[(N_CONF + 1):(N_CONF + N_Treat), (1:N_CONF)]
    inter_var <- (1 - c(inter_var))
    
    t_stat <- ( ( median(inter_var) - median(intra_var) ) / sqrt( var(inter_var)/length(inter_var) + var(intra_var)/length(intra_var) ) )
    t_stat
    
    tmp_df = data.frame(celltype=celltype_ofinterest, med_inter_var = median(inter_var), med_intra_var = median(intra_var), t_stat = t_stat, condition=compare_of_interest[2] )
    stat_df <- rbind(stat_df,tmp_df)
  }
  return(stat_df)
}

### Module 2
cmp_tstat_rnd_SH <- function (tmp_sce, celltype_ofinterest) {
  
  tmp_idx <- colData(tmp_sce)$celltype == celltype_ofinterest
  sce_rnd <- tmp_sce[,tmp_idx]
  
  cond_N <- table( colData(sce_rnd)$sample_name2 )
  CONF_N <- unname(cond_N[1])
  R2KR_N <- unname(cond_N[2])
  R2P3_N <- unname(cond_N[3])
  
  smpl_num <- c(1:dim(sce_rnd)[2])
  sce_rnd <- sce_rnd[,sample(smpl_num)]
  colData(sce_rnd)$sample_name2 <- vector( mode = "character", length = dim(sce_rnd)[2] )
  colData(sce_rnd)$sample_name2[1:CONF_N] <- "CONF"
  colData(sce_rnd)$sample_name2[(CONF_N + 1):(CONF_N + R2KR_N)] <- "R2KR"
  colData(sce_rnd)$sample_name2[(CONF_N + R2KR_N + 1):(CONF_N + R2KR_N + R2P3_N)] <- "R2P3"
  
  stat_df <- data.frame(
    med_inter_var = double(),
    med_intra_var = double(),
    t_stat = double(),
    cond = character()
  )
  
  for (i in c(1:2)) {
    compare_list <- list(c("CONF", "R2KR"), c("CONF", "R2P3"))
    compare_of_interest <- compare_list[[i]]
    
    tmp_idx <- colData(sce_rnd)$sample_name2 %in% compare_of_interest
    tmp_sce_1 <- sce_rnd[,tmp_idx]
    
    tmp_idx <- ( Matrix::rowSums(assays(tmp_sce_1)$counts) > 0 )
    tmp_sce_1 <- tmp_sce_1[tmp_idx,]
    
    tmp_sce_1 <- scater::arrange(tmp_sce_1, sample_name2)
    
    tmp_cor_1 <- cor(as.matrix( assays(tmp_sce_1)$logcounts ))
    
    N_CONF <- sum( colData(tmp_sce_1)$sample_name2 == "CONF" )
    N_Treat <- sum( colData(tmp_sce_1)$sample_name2 == compare_of_interest[2] )
    
    intra_var_CONF <- tmp_cor_1[(1:N_CONF), (1:N_CONF)]
    intra_var_Treat <- tmp_cor_1[(N_CONF + 1):(N_CONF + N_Treat), (N_CONF + 1):(N_CONF + N_Treat)]
    intra_var_CONF <- intra_var_CONF[lower.tri(intra_var_CONF)]
    intra_var_Treat <- intra_var_Treat[lower.tri(intra_var_Treat)]
    intra_var <- (1 - c(intra_var_CONF, intra_var_Treat))
    rm(intra_var_CONF)
    rm(intra_var_Treat)
    inter_var <- tmp_cor_1[(N_CONF + 1):(N_CONF + N_Treat), (1:N_CONF)]
    inter_var <- (1 - c(inter_var))
    
    t_stat <- ( ( median(inter_var) - median(intra_var) ) / sqrt( var(inter_var)/length(inter_var) + var(intra_var)/length(intra_var) ) )
    t_stat
    
    tmp_df = data.frame(med_inter_var = median(inter_var), med_intra_var = median(intra_var), t_stat = t_stat, cond=compare_of_interest[2] )
    stat_df <- rbind(stat_df,tmp_df)
  }
  return(stat_df)
}

### Module 3
cmp_nulldist_SH <- function (tmp_sce, celltype_ofinterest, N){
  TstatNull_df <- data.frame(celltype=character(),
                             t_stat=numeric(),         
                             condition=character(),
                             stringsAsFactors=FALSE)
  
  for (i in c(1:N)) {
    print(paste0("N=",i," in ",celltype_ofinterest))
    tmp_df <- cmp_tstat_rnd_SH(tmp_sce, celltype_ofinterest)
    
    tmp_df_KR <- data.frame(celltype = celltype_ofinterest, t_stat = tmp_df[tmp_df$cond == "R2KR",]$t_stat, condition = "R2KR", stringsAsFactors = FALSE)
    tmp_df_P3 <- data.frame(celltype = celltype_ofinterest, t_stat = tmp_df[tmp_df$cond == "R2P3",]$t_stat, condition = "R2P3", stringsAsFactors = FALSE)
    TstatNull_df <- rbind(TstatNull_df, tmp_df_KR)  
    TstatNull_df <- rbind(TstatNull_df, tmp_df_P3)     
  }
  return(TstatNull_df)
}

## Main function
N <- 20000 

tstat_pv_df <- data.frame(celltype=character(),
                          t_stat=numeric(),         
                          condition=character(),
                          pv=numeric(),
                          stringsAsFactors=FALSE)

celltypes <- names( table( colData(R2_sce_Mesen)$celltype ) ) 

for (celltype_ofinterest in celltypes) {
  
  tmp_sce <- R2_sce_Mesen 
  tmp_idx <- colData(tmp_sce)$celltype == celltype_ofinterest
  tmp_sce <- tmp_sce[,tmp_idx]
  
  tstatNull_df <- cmp_nulldist_SH(tmp_sce, celltype_ofinterest, N)
  
  tstat_df <- cmp_tstat_SH(tmp_sce, celltype_ofinterest)
  
  tmp_null_df <- tstatNull_df[tstatNull_df$condition == "R2KR",]
  tmp_df <- tstat_df[tstat_df$condition == "R2KR",]
  pv_KR <- sum(tmp_df$t_stat < tmp_null_df$t_stat) / length(tmp_null_df$t_stat)
  
  tmp_null_df <- tstatNull_df[tstatNull_df$condition == "R2P3",]
  tmp_df <- tstat_df[tstat_df$condition == "R2P3",]
  pv_P3 <- sum(tmp_df$t_stat < tmp_null_df$t_stat) / length(tmp_null_df$t_stat)
  
  tmp_tstat_pv_df <- data.frame(tstat_df, pv = c(pv_KR, pv_P3))
  tstat_pv_df <- rbind(tstat_pv_df, tmp_tstat_pv_df)
}

tstat_pv_Ms_df <- tstat_pv_df
rm(tstat_pv_df)
tstat_pv_Ms_df # This is p-value for the significance of transcriptomic change of mesenchymal cells between Confetti and Red2Onco models.



