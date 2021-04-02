##### Calculating the fractions of epithelial cell types in Red2Onco models and Confetti control.
##### Author: Seungmin Han (sh906@cam.ac.uk)
##### Last Update: 31/03/2021

library(scater)
library(scran)
library(ggplot2)
library(Seurat)
library(plyr)

result_path <-"../data"
setwd(result_file)

load("R2_sce_epi.RData")

colData(R2_sce_epi)['celltype_merged'] <- vector(mode = "character", length = dim(R2_sce_epi)[2])
colData(R2_sce_epi)['celltype_merged'] <- R2_sce_epi$celltype

tmp_idx <- colData(R2_sce_epi)$cluster_k6 %in% c('11', '22')
colData(R2_sce_epi)$celltype_merged[tmp_idx] <- c("stem_cell")
tmp_idx <- colData(R2_sce_epi)$cluster_k6 %in% c('14', '15', '19', '20', '23', '29', '30')
colData(R2_sce_epi)$celltype_merged[tmp_idx] <- c("TA_cell")
tmp_idx <- colData(R2_sce_epi)$cluster_k6 %in% c('5', '18', '21', '32', '33','6', '8', '13', '31', '32', '33')
colData(R2_sce_epi)$celltype_merged[tmp_idx] <- c("EP_Ent")
tmp_idx <- colData(R2_sce_epi)$cluster_k6 %in% c('4', '16', '25', '26')
colData(R2_sce_epi)$celltype_merged[tmp_idx] <- c("Paneth_cell")
tmp_idx <- colData(R2_sce_epi)$cluster_k6 %in% c('1', '2', '3', '9', '10', '27')
colData(R2_sce_epi)$celltype_merged[tmp_idx] <- c("Goblet_cell")
tmp_idx <- colData(R2_sce_epi)$cluster_k6 %in% c('7', '17')
colData(R2_sce_epi)$celltype_merged[tmp_idx] <- c("Tuft_cell")
tmp_idx <- colData(R2_sce_epi)$cluster_k6 %in% c('12', '24', '28')
colData(R2_sce_epi)$celltype_merged[tmp_idx] <- c("Enteroendocrine_cell")

yellow_df <- data.frame(cellcount_type=integer(),
                        cellcount_total=integer(),
                        samplename=character(),
                        celltype=character(),
                        stringsAsFactors=FALSE)

red_df <- data.frame(cellcount_type=integer(),
                     cellcount_total=integer(),
                     samplename=character(),
                     celltype=character(),
                     stringsAsFactors=FALSE)

samples <- c("CONF_2", "CONF_3", "R2KR_1", "R2KR_2", "R2KR_3", "R2P3_1", "R2P3_2")
celltypes <- c("stem_cell", "TA_cell", "Paneth_cell", "Goblet_cell", "Tuft_cell", "Enteroendocrine_cell", "EP_Ent")

for (i in celltypes) {
  for (j in samples) {
    tmp_idx <- (colData(R2_sce_epi)$sample_name == j) & (colData(R2_sce_epi)$EYFP_idx == 1)
    tmp_idx2 <- tmp_idx & (colData(R2_sce_epi)$celltype_merged == i)
    tmp_df <- data.frame(cellcount_type = sum(tmp_idx2), cellcount_total = sum(tmp_idx), samplecondition = strsplit(j,"_")[[1]][1], celltype = i, stringsAsFactors=FALSE)
    yellow_df <- rbind(yellow_df, tmp_df)
  }
}

for (i in celltypes) {
  for (j in samples) {
    tmp_idx <- (colData(R2_sce_epi)$sample_name == j) & (colData(R2_sce_epi)$RFP_idx == 1)
    tmp_idx2 <- tmp_idx & (colData(R2_sce_epi)$celltype_merged == i)
    tmp_df <- data.frame(cellcount_type = sum(tmp_idx2), cellcount_total = sum(tmp_idx), samplecondition = strsplit(j,"_")[[1]][1], celltype = i, stringsAsFactors=FALSE)
    red_df <- rbind(red_df, tmp_df)
  }
}

yellow_df2 <- cbind(yellow_df, ratio = yellow_df[ ,1]/yellow_df[ ,2])
red_df2 <- cbind(red_df, ratio = red_df[ ,1]/red_df[ ,2])

cellcounttable <- yellow_df
cellcounttable$samplecondition <- factor(cellcounttable$samplecondition, levels = c("CONF", "R2KR", "R2P3"))
cellcounttable$celltype <- factor(cellcounttable$celltype, levels = c("stem_cell", "TA_cell", "Paneth_cell", "Goblet_cell", "Tuft_cell", "Enteroendocrine_cell", "EP_Ent"))
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  datac <- rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

celltyperatio = cellcounttable[,1] / cellcounttable[,2]
cellcounttable = cbind(cellcounttable,celltyperatio)
cellcount_summary <- summarySE(cellcounttable, measurevar = "celltyperatio", groupvars = c("celltype","samplecondition"))
Csummary <- cellcount_summary

ggplot(Csummary, aes(x=celltype, y=celltyperatio, fill=samplecondition)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=celltyperatio-se, ymax=celltyperatio+se),
                width=.2,                    
                position=position_dodge(.9)) +
  xlab("Cell type") +
  ylab("Cell type fraction") +
  guides(fill=guide_legend(title="Condition")) +
  coord_fixed(ratio = 10) +
  theme(
    axis.text.x = element_text(angle = 90),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("Cell type fraction in YFP+ epithelial cells")

df <- yellow_df
Pr_df <- matrix(0, nrow = length(celltypes), ncol = 2)
colnames(Pr_df) <- c("CONF_vs_R2KR", "CONF_vs_R2P3")
celltypes <- c("stem_cell", "TA_cell", "Paneth_cell", "Goblet_cell", "Tuft_cell", "Enteroendocrine_cell", "EP_Ent")
rownames(Pr_df) <- celltypes

for (celltype in celltypes) {
  for (of_interest in c("R2KR", "R2P3")) {
    comparison = list("R2KR" = c("CONF", "R2KR"), "R2P3" = c("CONF", "R2P3"))
    tmp_df <- df[( df$celltype == celltype & df$samplecondition %in% comparison[[of_interest]] ), ]
    tmp_df <- cbind(tmp_df[ ,1:3], lcases = log(tmp_df$cellcount_total) )
    tmp_df
    log.fit = glm(cellcount_type ~ samplecondition + offset(lcases), family=poisson, data=tmp_df)
    summary(log.fit)
    fitted(log.fit)
    tmp_anova <- anova(log.fit,dispersion = NULL, test = "LRT") 
    Pr_df[celltype, paste("CONF_vs_",of_interest, sep="")] <- tmp_anova$Pr[2]
  }
}

pv_YFP_df <- Pr_df 
rm(Pr_df)
pv_YFP_df # p-value for YFP+ epithelial cells from likelihood ratio test.

cellcounttable <- red_df
cellcounttable$samplecondition <- factor(cellcounttable$samplecondition, levels = c("CONF", "R2KR", "R2P3"))
cellcounttable$celltype <- factor(cellcounttable$celltype, levels = c("stem_cell", "TA_cell", "Paneth_cell", "Goblet_cell", "Tuft_cell", "Enteroendocrine_cell", "EP_Ent"))
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  datac <- rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}
celltyperatio = cellcounttable[,1] / cellcounttable[,2]
cellcounttable = cbind(cellcounttable,celltyperatio)
cellcount_summary <- summarySE(cellcounttable, measurevar = "celltyperatio", groupvars = c("celltype","samplecondition"))
Csummary <- cellcount_summary

ggplot(Csummary, aes(x=celltype, y=celltyperatio, fill=samplecondition)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=celltyperatio-se, ymax=celltyperatio+se),
                width=.2,                    
                position=position_dodge(.9)) +
  xlab("Cell type") +
  ylab("Cell type fraction") +
  guides(fill=guide_legend(title="Condition")) +
  coord_fixed(ratio = 10) +
  theme(
    axis.text.x = element_text(angle = 90),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("Cell type fractions in RFP+ epithelial cells")

df <- red_df
Pr_df <- matrix(0, nrow = length(celltypes), ncol = 2)
colnames(Pr_df) <- c("CONF_vs_R2KR", "CONF_vs_R2P3")
celltypes <- c("stem_cell", "TA_cell", "Paneth_cell", "Goblet_cell", "Tuft_cell", "Enteroendocrine_cell", "EP_Ent")
rownames(Pr_df) <- celltypes
for (celltype in celltypes) {
  for (of_interest in c("R2KR", "R2P3")) {
    comparison = list("R2KR" = c("CONF", "R2KR"), "R2P3" = c("CONF", "R2P3"))
    tmp_df <- df[( df$celltype == celltype & df$samplecondition %in% comparison[[of_interest]] ), ]
    tmp_df <- cbind(tmp_df[ ,1:3], lcases = log(tmp_df$cellcount_total) )
    tmp_df
    log.fit = glm(cellcount_type ~ samplecondition + offset(lcases), family=poisson, data=tmp_df)
    summary(log.fit)
    fitted(log.fit)
    tmp_anova <- anova(log.fit,dispersion = NULL, test = "LRT") 
    Pr_df[celltype, paste("CONF_vs_",of_interest, sep="")] <- tmp_anova$Pr[2]
  }
}

pv_RFP_df <- Pr_df 
rm(Pr_df)
pv_RFP_df # p-value for RFP+ epithelial cells from likelihood ratio test.
