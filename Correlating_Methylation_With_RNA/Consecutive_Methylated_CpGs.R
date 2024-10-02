library(dplyr)

source("~/Consecutive_CpG_Finder.R")

tab <- read.table("~/Embryo_vs_Larva_10%_sig.txt")

streaks <- consecutive_CpG_finder(tab, include_non_consecutive = FALSE)
View(streaks)

setwd("~/Closest_TSS")

fet <- read.table("GCF_009193385.2_Nvit_psr_1.1_feature_table.txt", sep = "\t",header = F)
colnames(fet) <- c("feature", "class", "assembly","assembly_unit","seq_type","chromosome",
                   "genomic_accession","start","end","strand","product_accession","non-redundant_refseq",
                   "related_accession","name","symbol","GeneID","locus_tag","feature_interval_length",
                   "product_length","attributes")
bed <- read.table("~/GCF_009193385.2_Nvit_psr_1.1_genomic.bed.txt")
colnames(bed) <- c("chr","start","end","name","score","strand","thickStart","thickEnd","itemRgb","blockCount",
                   "blockSizes","blockStarts")

bed$transcript<-gsub("rna-","",as.character(bed$name))

combined <- sqldf("SELECT trans.transcript,
                trans.thickStart,
                gene.genomic_accession,
                gene.product_accession,
                gene.symbol
                FROM bed AS trans
                LEFT JOIN fet AS gene
                ON trans.transcript = gene.product_accession")

head(combined)
cols <- c(1,2,3,5)
combined <- combined[cols]
colnames(combined) <- c("transcript","TSS","chr","gene")

outputEvL <- Assigning_Closest_TSS(combined,streaks)

tab <- read.table("~/Larva_vs_Prepupae_10%_sig.txt")
streaks <- consecutive_CpG_finder(tab, include_non_consecutive = FALSE)
outputLvP <- Assigning_Closest_TSS(combined,streaks)

tab <- read.table("~/Prepupae_vs_Pupae_10%_sig.txt")
streaks <- consecutive_CpG_finder(tab, include_non_consecutive = FALSE)
outputPvP <- Assigning_Closest_TSS(combined,streaks)

tab <- read.table("~/Pupae_vs_Adult_10%_sig.txt")
streaks <- consecutive_CpG_finder(tab, include_non_consecutive = FALSE)
outputPvA <- Assigning_Closest_TSS(combined,streaks)

EvLDEGs <- read_delim("~/EmbryovsLarva_Gene_results.txt",
                      "\t", escape_double = FALSE, trim_ws = TRUE)

colnames(EvLDEGs) <- c("baseMean", "log2lfc","lfcSE","stat","p-value","FDR","gene")

#output <- output  %>%  group_by(gene) %>% #This makes fit worse
#  filter(distance == min(distance))


EvL_Streaks <- EvLDEGs[EvLDEGs$gene %in% outputEvL$gene,]

EvL_data <- merge(EvL_Streaks,outputEvL,by="gene")
EvL_data$stage <- "embryo_vs_larva"
EvL_data

LvPDEGs <- read_delim("~/LarvavsPrepupae_Gene_results.txt",
                      "\t", escape_double = FALSE, trim_ws = TRUE)

colnames(LvPDEGs) <- c("baseMean", "log2lfc","lfcSE","stat","p-value","FDR","gene")

#output <- output  %>%  group_by(gene) %>% #This makes fit worse
#  filter(distance == min(distance))


LvP_Streaks <- LvPDEGs[LvPDEGs$gene %in% outputLvP$gene,]

LvP_data <- merge(LvP_Streaks,outputLvP,by="gene")
LvP_data$stage <- "larva_vs_prepupa"
LvP_data

PvPDEGs <- read_delim("~/PrepupaevsPupae_Gene_results.txt",
                      "\t", escape_double = FALSE, trim_ws = TRUE)

colnames(PvPDEGs) <- c("baseMean", "log2lfc","lfcSE","stat","p-value","FDR","gene")

#output <- output  %>%  group_by(gene) %>% #This makes fit worse
#  filter(distance == min(distance))


PvP_Streaks <- PvPDEGs[PvPDEGs$gene %in% outputPvP$gene,]

PvP_data <- merge(PvP_Streaks,outputPvP,by="gene")
PvP_data$stage <- "prepupa_vs_pupa"
PvP_data

PvADEGs <- read_delim("~/PupaevsAdult_Gene_results.txt",
                      "\t", escape_double = FALSE, trim_ws = TRUE)

colnames(PvADEGs) <- c("baseMean", "log2lfc","lfcSE","stat","p-value","FDR","gene")

#output <- output  %>%  group_by(gene) %>% #This makes fit worse
#  filter(distance == min(distance))


PvA_Streaks <- PvADEGs[PvADEGs$gene %in% outputPvA$gene,]

PvA_data <- merge(PvA_Streaks,outputPvA,by="gene")
PvA_data$stage <- "pupa_vs_adult"
PvA_data

results <- rbind(EvL_data,LvP_data,PvP_data,PvA_data)

################################################################################
mean(results$log2lfc)
var(results$log2lfc)

hist(results$log2lfc) #is normal
model <- glm(log2lfc~distance*stage*consecutive, family = gaussian, data=results, 
              model = TRUE, method = "glm.fit2")

library(car)
#Report the interaction stats from Anova. You are doing a GLM with a quasipoisson distribution
Anova(model)
#summary (model2)

model$deviance/model$df.residual #Overdispersion > 1.2 bad

#pseudo R2
100*((model$null.deviance-model$deviance)/model$null.deviance) #pseudo R2

emmeans(model, pairwise~stage*consecutive )
