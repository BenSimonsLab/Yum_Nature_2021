##### Calculating AUC scores of BMP signaling for epithelial cells
##### Author: Seungmin Han (sh906@cam.ac.uk)
##### Last Update: 31/03/2021

library(scater)
library(scran)
library(ggplot2)
library(Seurat)
library(stats)
library(plyr)
library(RColorBrewer)
library(AUCell)
library(GSEABase)

result_path <-"../data"
setwd(result_file)
load("R2_sce_epi.RData")

geneset_file1 = '../data/Qi_BMP_targetgenes.gmt'
geneSets1 <- getGmt(geneset_file1)

tmp_sce <- R2_sce_epi
geneSets_BMP <- GeneSetCollection(c(geneSets1))
geneSets_BMP <- subsetGeneSets(geneSets_BMP, toupper(rownames(tmp_sce))) 
rownames(tmp_sce) <- toupper(rownames(tmp_sce))

cells_rankings <- AUCell_buildRankings(assays(tmp_sce)$counts, nCores=10, plotStats=TRUE)
cells_AUC_BMP <- AUCell_calcAUC(geneSets_BMP, cells_rankings, aucMaxRank = ceiling(0.3 * nrow(cells_rankings)))    

tmp_sce <- R2_sce_epi
set.seed(123)
par(mfrow=c(1,2)) 
cells_assignment_BMP <- AUCell_exploreThresholds(cells_AUC_BMP, plotHist=TRUE, assign=TRUE)

AUCthr <- c(0.1653310)
samples <- c("bothcolor-CONF_2", "bothcolor-CONF_3", "yellow-R2KR_1", "yellow-R2KR_2", "yellow-R2KR_3", "red-R2KR_1", "red-R2KR_2", "red-R2KR_3", "yellow-R2P3_1", "yellow-R2P3_2", "red-R2P3_1", "red-R2P3_2")
geneSets <- geneSets_BMP

colData(tmp_sce)["tmp_sample_name"] <- vector( mode = "character", length = dim(tmp_sce)[2] )
tmp_idx <- colData(tmp_sce)$RFP_idx == 1
tmp_idx_con <- colData(tmp_sce)$sample_name %in% c("R2KR_1","R2KR_2","R2KR_3","R2P3_1","R2P3_2")
colData(tmp_sce)$tmp_sample_name[tmp_idx & tmp_idx_con] <- "red"
colData(tmp_sce)$tmp_sample_name[(!tmp_idx) & tmp_idx_con] <- "yellow"
colData(tmp_sce)$tmp_sample_name[!tmp_idx_con] <- "bothcolor"

colData(tmp_sce)$tmp_sample_name <- paste(colData(tmp_sce)$tmp_sample_name, colData(tmp_sce)$sample_name, sep="-")

AUC_df <- data.frame(cellcount_AUC=integer(),
                        cellcount_total=integer(),
                        samplename=character(),
                        geneset=character(),
                        stringsAsFactors=FALSE)

names(AUCthr) <- rownames(cbind(nGenes(geneSets)))
Genesets_name <- rownames(cbind(nGenes(geneSets)))

for (i in Genesets_name) {
  for (j in samples) {
    tmp_idx <- (colData(tmp_sce)$tmp_sample_name == j)
    tmp_idx2 <- tmp_idx & unname( assays(cells_AUC_BMP)[[1]][i,] > unname(AUCthr[i]) )
    tmp_df <- data.frame(cellcount_AUC = sum(tmp_idx2), cellcount_total = sum(tmp_idx), samplecondition = strsplit(j,"_")[[1]][1], geneset = i, stringsAsFactors=FALSE)
    AUC_df <- rbind(AUC_df, tmp_df)
  }
}

cellcounttable <- AUC_df
cellcounttable$samplecondition <- factor(cellcounttable$samplecondition, levels = c("bothcolor-CONF", "yellow-R2KR","red-R2KR", "yellow-R2P3", "red-R2P3"))
cellcounttable$geneset <- factor( cellcounttable$geneset, levels = rownames(cbind(nGenes(geneSets))) )
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
cellcount_summary <- summarySE(cellcounttable, measurevar = "celltyperatio", groupvars = c("geneset","samplecondition"))
Csummary <- cellcount_summary

ggplot(Csummary, aes(x=geneset, y=celltyperatio, fill=samplecondition)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=celltyperatio-se, ymax=celltyperatio+se),
                width=.2,                    
                position=position_dodge(.9)) +
  ggtitle("GSEA for BMP target genes in epithelial cells") + 
  xlab("Gene set") +
  ylab("Fraction of cells with active gene set") +
  guides(fill=guide_legend(title="Condition")) +
  theme(
    axis.text.x = element_blank()
    
  )

df <- AUC_df

Pr_df <- matrix(0, nrow = length(table(df$geneset)), ncol = 4)
colnames(Pr_df) <- c("CONF vs Y_R2KR", "CONF vs R_R2KR", "CONF vs Y_R2P3", "CONF vs R_R2P3")
genesets <- names(table(df$geneset))
rownames(Pr_df) <- genesets

for (geneset in genesets) {
  for (of_interest in c("Y_R2KR", "R_R2KR", "Y_R2P3","R_R2P3")) {
    comparison = list("Y_R2KR" = c("bothcolor-CONF", "yellow-R2KR"), "R_R2KR" = c("bothcolor-CONF", "red-R2KR"),"Y_R2P3" = c("bothcolor-CONF", "yellow-R2P3"), "R_R2P3" = c("bothcolor-CONF", "red-R2P3"))
    tmp_df <- df[( df$samplecondition %in% comparison[[of_interest]] ), ]
    tmp_df <- cbind(tmp_df[ ,1:3], lcases = log(tmp_df$cellcount_total) )
    tmp_df

    log.fit = glm(cellcount_AUC ~ samplecondition + offset(lcases), family=poisson, data=tmp_df)
    summary(log.fit)
    fitted(log.fit)
    
    tmp_anova <- anova(log.fit,dispersion = NULL, test = "LRT") 
    Pr_df[geneset, paste("CONF vs ",of_interest, sep="")] <- tmp_anova$Pr[2]
  }
}

Pr_df # p-value from likelihood ratio test
