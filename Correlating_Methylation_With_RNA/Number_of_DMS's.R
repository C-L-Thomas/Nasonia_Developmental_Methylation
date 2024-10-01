library(dplyr)
library(tidyverse)

#Read in Differentially Methylated Sites

EvLDMSs <- read_delim("~/Embryo_vs_Larva_DMS_w10%diff_with_geneID.txt",
                      "\t", escape_double = FALSE, trim_ws = TRUE)

colnames(EvLDMSs) <- c("chr", "start","end","p-value","qvalue","meth_diff","chr2","start2","end2","gene")

EvLDEGs <- read_delim("~/EmbryovsLarva_Gene_results.txt",
                      "\t", escape_double = FALSE, trim_ws = TRUE)

colnames(EvLDEGs) <- c("baseMean", "log2lfc","lfcSE","stat","p-value","FDR","gene")


EvL_gene_quantity <- as.data.frame(table(EvLDMSs$gene))
colnames(EvL_gene_quantity) <- c("gene","DMS")
nrow(EvL_gene_quantity)
EvLGeneList <- EvLDEGs[EvLDEGs$gene %in% EvL_gene_quantity$gene,]
EvLGeneList
EvL <- merge(EvL_gene_quantity,EvLGeneList,by="gene")
EvL$stage <- "embryo_vs_larva"

LvPDMSs <- read_delim("~/Larva_vs_Prepupae_DMS_w10%diff_with_geneID.txt",
                      "\t", escape_double = FALSE, trim_ws = TRUE)

colnames(LvPDMSs) <- c("chr", "start","end","p-value","qvalue","meth_diff","chr2","start2","end2","gene")

LvPDEGs <- read_delim("~/LarvavsPrepupae_Gene_results.txt",
                      "\t", escape_double = FALSE, trim_ws = TRUE)

colnames(LvPDEGs) <- c("baseMean", "log2lfc","lfcSE","stat","p-value","FDR","gene")


LvP_gene_quantity <- as.data.frame(table(LvPDMSs$gene))
colnames(LvP_gene_quantity) <- c("gene","DMS")
nrow(LvP_gene_quantity)
LvPGeneList <- LvPDEGs[LvPDEGs$gene %in% LvP_gene_quantity$gene,]
LvPGeneList
LvP <- merge(LvP_gene_quantity,LvPGeneList,by="gene")
LvP$stage <- "larva_vs_prepupa"

PvPDMSs <- read_delim("~/Prepupae_vs_Pupae_DMS_w10%diff_with_geneID.txt",
                      "\t", escape_double = FALSE, trim_ws = TRUE)

colnames(PvPDMSs) <- c("chr", "start","end","p-value","qvalue","meth_diff","chr2","start2","end2","gene")

PvPDEGs <- read_delim("~/PrepupaevsPupae_Gene_results.txt",
                      "\t", escape_double = FALSE, trim_ws = TRUE)

colnames(PvPDEGs) <- c("baseMean", "log2lfc","lfcSE","stat","p-value","FDR","gene")


PvP_gene_quantity <- as.data.frame(table(PvPDMSs$gene))
colnames(PvP_gene_quantity) <- c("gene","DMS")
nrow(PvP_gene_quantity)
PvPGeneList <- PvPDEGs[PvPDEGs$gene %in% PvP_gene_quantity$gene,]
PvPGeneList
PvP <- merge(PvP_gene_quantity,PvPGeneList,by="gene")
PvP$stage <- "prepupa_vs_pupa"

PvADMSs <- read_delim("~/Pupae_vs_Adult_DMS_w10%diff_with_geneID.txt",
                      "\t", escape_double = FALSE, trim_ws = TRUE)

colnames(PvADMSs) <- c("chr", "start","end","p-value","qvalue","meth_diff","chr2","start2","end2","gene")

PvADEGs <- read_delim("~/PupaevsAdult_Gene_results.txt",
                      "\t", escape_double = FALSE, trim_ws = TRUE)

colnames(PvADEGs) <- c("baseMean", "log2lfc","lfcSE","stat","p-value","FDR","gene")


PvA_gene_quantity <- as.data.frame(table(PvADMSs$gene))
colnames(PvA_gene_quantity) <- c("gene","DMS")
nrow(PvA_gene_quantity)
PvAGeneList <- PvADEGs[PvADEGs$gene %in% PvA_gene_quantity$gene,]
PvAGeneList
PvA <- merge(PvA_gene_quantity,PvAGeneList,by="gene")
PvA$stage <- "pupa_vs_adult"

results <- rbind(EvL,LvP,PvP,PvA)

################################################################################
library(glm2)

hist(results$log2lfc) #is normal
model <- glm2(log2lfc~DMS*stage, family = gaussian, data=results,
              model = TRUE, method = "glm.fit2")

library(car)
#Report the interaction stats from Anova. You are doing a GLM with a quasipoisson distribution
Anova(model)
#summary (model2)

model$deviance/model$df.residual #Overdispersion > 1.2 bad

#pseudo R2
100*((model$null.deviance-model$deviance)/model$null.deviance) #pseudo R2
