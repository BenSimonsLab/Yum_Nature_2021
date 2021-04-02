##### Calculating the fractions of immune cell types in Red2Onco models and Confetti control.
##### Author: Seungmin Han (sh906@cam.ac.uk)
##### Last Update: 31/03/2021

library(scater)
library(scran)
library(Seurat)
library(plyr)

result_path <-"../data"
setwd(result_file)

load("R2_sce_Immune.RData")

immune_df <- data.frame(cellcount_type=integer(),
                        cellcount_total=integer(),
                        samplename=character(),
                        celltype=character(),
                        stringsAsFactors=FALSE)

samples <- c("CONF_2", "CONF_3", "R2KR_1", "R2KR_2", "R2KR_3", "R2P3_1", "R2P3_2")
celltypes <- c("B_cell","DC","Monocyte","MP_1","MP_2","Plasma_cell","T_cell")

for (i in celltypes) {
  for (j in samples) {
    tmp_idx <- (colData(R2_sce_Immune)$sample_name == j) 
    tmp_idx2 <- tmp_idx & (colData(R2_sce_Immune)$celltype == i)
    tmp_df <- data.frame(cellcount_type = sum(tmp_idx2), cellcount_total = sum(tmp_idx), samplecondition = strsplit(j,"_")[[1]][1], celltype = i, stringsAsFactors=FALSE)
    immune_df <- rbind(immune_df, tmp_df)
  }
  
}

cellcounttable <- immune_df
cellcounttable$samplecondition <- factor(cellcounttable$samplecondition, levels = c("CONF", "R2KR", "R2P3"))
cellcounttable$celltype <- factor(cellcounttable$celltype, levels = c("B_cell","DC","Monocyte","MP_1","MP_2","Plasma_cell","T_cell"))
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
  theme(
    axis.text.x = element_text(angle = 90),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("Cell type fraction in immune cells")

df <- immune_df

Pr_df <- matrix(0, nrow = length(celltypes), ncol = 2)
colnames(Pr_df) <- c("CONF_vs_R2KR", "CONF_vs_R2P3")
celltypes <- c("B_cell","DC","Monocyte","MP_1","MP_2","Plasma_cell","T_cell")
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

Pr_df # p-value from likelihood ratio test
