library(dplyr)
library(tidyverse)

#Read in Differentially Methylated Sites

EvLDMSs <- read_delim("~/Embryo_vs_Larva_10%_sig_TFBS_distance.txt",
                      "\t", escape_double = FALSE, trim_ws = TRUE)

EvLDEGs <- read_delim("~/EmbryovsLarva_Gene_results.txt",
                      "\t", escape_double = FALSE, trim_ws = TRUE)

colnames(EvLDEGs) <- c("baseMean", "log2lfc","lfcSE","stat","p-value","FDR","gene")

EvLdf<- data.frame(EvLDMSs)

EvLclosest <- EvLdf  %>%  group_by(symbol) %>%
  filter(disttss == min(disttss))

#EvLclosest <- EvLdf #Quick way to boycot above

nrow(EvLclosest)
colnames(EvLclosest) <- c("disttss", "strand", "location", "genomic_accession", "gene","TF" )

EvL_DEGs_in_Closest_TFBS <- EvLDEGs[EvLDEGs$gene %in% EvLclosest$gene,]
EvL_DEGs_in_Closest_TFBS
EvL_TFBS_data <- merge(EvL_DEGs_in_Closest_TFBS,EvLclosest,by="gene")
EvL_TFBS_data$stage <- "embryo_vs_larva"

LvPDMSs <- read_delim("~/Larva_vs_Prepupae_10%_sig_TFBS_distance.txt",
                      "\t", escape_double = FALSE, trim_ws = TRUE)

LvPDEGs <- read_delim("~/LarvavsPrepupae_Gene_results.txt",
                      "\t", escape_double = FALSE, trim_ws = TRUE)

colnames(LvPDEGs) <- c("baseMean", "log2lfc","lfcSE","stat","p-value","FDR","gene")

LvPdf<- data.frame(LvPDMSs)

LvPclosest <- LvPdf  %>%  group_by(symbol) %>%
  filter(disttss == min(disttss))

# LvPclosest <- LvPdf #Quick way to boycot above

nrow(LvPclosest)
colnames(LvPclosest) <- c("disttss", "strand", "location", "genomic_accession", "gene","TF" )

LvP_DEGs_in_Closest_TFBS <- LvPDEGs[LvPDEGs$gene %in% LvPclosest$gene,]
LvP_DEGs_in_Closest_TFBS
LvP_TFBS_data <- merge(LvP_DEGs_in_Closest_TFBS,LvPclosest,by="gene")
LvP_TFBS_data$stage <- "larva_vs_prepupa"

PvPDMSs <- read_delim("~/Prepupae_vs_Pupae_10%_sig_TFBS_distance.txt",
                      "\t", escape_double = FALSE, trim_ws = TRUE)

PvPDEGs <- read_delim("~/PrepupaevsPupae_Gene_results.txt",
                      "\t", escape_double = FALSE, trim_ws = TRUE)

colnames(PvPDEGs) <- c("baseMean", "log2lfc","lfcSE","stat","p-value","FDR","gene")

PvPdf<- data.frame(PvPDMSs)

PvPclosest <- PvPdf  %>%  group_by(symbol) %>%
  filter(disttss == min(disttss))

#PvPclosest <- PvPdf #Quick way to boycot above

nrow(PvPclosest)
colnames(PvPclosest) <- c("disttss", "strand", "location", "genomic_accession", "gene","TF" )

PvP_DEGs_in_Closest_TFBS <- PvPDEGs[PvPDEGs$gene %in% PvPclosest$gene,]
PvP_DEGs_in_Closest_TFBS
PvP_TFBS_data <- merge(PvP_DEGs_in_Closest_TFBS,PvPclosest,by="gene")
PvP_TFBS_data$stage <- "prepupa_vs_pupa"

PvADMSs <- read_delim("~/Pupae_vs_Adult_10%_sig_TFBS_distance.txt",
                      "\t", escape_double = FALSE, trim_ws = TRUE)

PvADEGs <- read_delim("~/PupaevsAdult_Gene_results.txt",
                      "\t", escape_double = FALSE, trim_ws = TRUE)

colnames(PvADEGs) <- c("baseMean", "log2lfc","lfcSE","stat","p-value","FDR","gene")

PvAdf<- data.frame(PvADMSs)

PvAclosest <- PvAdf  %>%  group_by(symbol) %>%
  filter(disttss == min(disttss))

# PvAclosest <- PvAdf #Quick way to boycot above

nrow(PvAclosest)
colnames(PvAclosest) <- c("disttss", "strand", "location", "genomic_accession", "gene","TF" )

PvA_DEGs_in_Closest_TFBS <- PvADEGs[PvADEGs$gene %in% PvAclosest$gene,]
PvA_DEGs_in_Closest_TFBS
PvA_TFBS_data <- merge(PvA_DEGs_in_Closest_TFBS,PvAclosest,by="gene")
PvA_TFBS_data$stage <- "pupa_vs_adult"

results <- rbind(EvL_TFBS_data,LvP_TFBS_data,PvP_TFBS_data,PvA_TFBS_data)

table(results$TF)

################################################################################
hist(results$log2lfc) #is normal
model <- glm(log2lfc~disttss*stage*TF, family = gaussian, data=results,
               model = TRUE, method = "glm.fit2")

library(car)
#Report the interaction stats from Anova. You are doing a GLM with a quasipoisson distribution
Anova(model)
#summary (model2)

model$deviance/model$df.residual #Overdispersion > 1.2 bad

#pseudo R2
100*((model$null.deviance-model$deviance)/model$null.deviance) #pseudo R2

emtrends(model, pairwise~DStage, var = "Percentage")
