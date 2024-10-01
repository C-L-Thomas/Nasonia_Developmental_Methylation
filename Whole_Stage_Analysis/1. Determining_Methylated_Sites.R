library(methylKit)
library(grid)
library(readr)
library(ggplot2)
library(sqldf)

file.list <- list("DADU.bismark.cov.gz","DEMB.bismark.cov.gz","DLAR.bismark.cov.gz","DPRE.bismark.cov.gz","DPUP.bismark.cov.gz")

treat_conditions <- c(0,1,2,3,4)

raw_data <- methRead(file.list,
                     sample.id = list("ADU","EMB","LAR","PRE","PUP"),
                     treatment = treat_conditions,
                     assembly="Nvit_PSR_1.1",
                     dbtype = NA,
                     pipeline = "bismarkCoverage",
                     header = T, sep = "\t", mincov=1,
                     context="CpG")

#Do I want to set min and max thresholds, eg 10 reads, remove top 0.1%?
meth_all_data <- unite(raw_data, destrand=FALSE, min.per.group = 1L)

nrow(meth_all_data)
tail(meth_all_data, 100)

uniteDB.all <- getData(meth_all_data)

uniteDB.all$rownums <- row.names(uniteDB.all)

head(uniteDB.all)

adu <- uniteDB.all[,c(1,2,5,6)]
emb <- uniteDB.all[,c(1,2,8,9)]
lar <- uniteDB.all[,c(1,2,11,12)]
pre <- uniteDB.all[,c(1,2,14,15)]
pup <- uniteDB.all[,c(1,2,17,18)]

head(adu)

#These are averages from the individual lambda sites
adubt <- function(a, b, p = 0.0053) {binom.test(a, b, 0.0053, alternative="greater")$p.value}
embt <- function(a, b, p = 0.009) {binom.test(a, b, 0.009, alternative="greater")$p.value}
larbt <- function(a, b, p = 0.0087) {binom.test(a, b, 0.0087, alternative="greater")$p.value}
prebt <- function(a, b, p = 0.003) {binom.test(a, b, 0.003, alternative="greater")$p.value}
pupbt <- function(a, b, p = 0.005) {binom.test(a, b, 0.005, alternative="greater")$p.value}

colnames(adu) <- c("chr", "start","coverage","c_count")
adu <- adu[complete.cases(adu),]
adu$pVal <- mapply(adubt, adu$c_count, adu$coverage)
adu$FDR <- p.adjust(adu$pVal, method = "BH", n = length(adu$pVal))
adu$FDR <- as.numeric(adu$FDR)
write.table(adu, file="Adult_Full_output.txt", sep="\t", row.names=F, quote=F)
adu_meth_sites <- subset(adu, adu$FDR < 0.05)
nrow(adu_meth_sites)
write.table(adu_meth_sites, file="Adult_Significant_output.txt", sep="\t", row.names=F, quote=F)

colnames(emb) <- c("chr", "start","coverage","c_count")
emb <- emb[complete.cases(emb),]
emb$pVal <- mapply(embt, emb$c_count, emb$coverage)
emb$FDR <- p.adjust(emb$pVal, method = "BH", n = length(emb$pVal))
emb$FDR <- as.numeric(emb$FDR)
write.table(emb, file="Embryo_Full_output.txt", sep="\t", row.names=F, quote=F)
emb_meth_sites <- subset(emb, emb$FDR < 0.05)
nrow(emb_meth_sites)

write.table(emb_meth_sites, file="Embryo_Significant_output.txt",
            sep="\t", row.names=F, quote=F)

colnames(lar) <- c("chr", "start","coverage","c_count")
lar <- lar[complete.cases(lar),]
lar$pVal <- mapply(larbt, lar$c_count, lar$coverage)
lar$FDR <- p.adjust(lar$pVal, method = "BH", n = length(lar$pVal))
lar$FDR <- as.numeric(lar$FDR)
write.table(lar, file="Larva_Full_output.txt", sep="\t", row.names=F, quote=F)
lar_meth_sites <- subset(lar, lar$FDR < 0.05)
nrow(lar_meth_sites)

write.table(lar_meth_sites, file="Larva_Significant_output.txt",
            sep="\t", row.names=F, quote=F)

colnames(pre) <- c("chr", "start","coverage","c_count")
pre <- pre[complete.cases(pre),]
pre$pVal <- mapply(prebt, pre$c_count, pre$coverage)
pre$FDR <- p.adjust(pre$pVal, method = "BH", n = length(pre$pVal))
pre$FDR <- as.numeric(pre$FDR)
write.table(pre, file="Prepupae_Full_output.txt", sep="\t", row.names=F, quote=F)
pre_meth_sites <- subset(pre, pre$FDR < 0.05)
nrow(pre_meth_sites)

write.table(pre_meth_sites, file="Prepupae_Significant_output.txt",
            sep="\t", row.names=F, quote=F)

colnames(pup) <- c("chr", "start","coverage","c_count")
pup <- pup[complete.cases(pup),]
pup$pVal <- mapply(pupbt, pup$c_count, pup$coverage)
pup$FDR <- p.adjust(pup$pVal, method = "BH", n = length(pup$pVal))
pup$FDR <- as.numeric(pup$FDR)
write.table(pup, file="Pupae_Full_output.txt", sep="\t", row.names=F, quote=F)
pup_meth_sites <- subset(pup, pup$FDR < 0.05)
nrow(pup_meth_sites)

write.table(pup_meth_sites, file="Pupae_Significant_output.txt",
            sep="\t", row.names=F, quote=F)

save.image(file = "Stage_Level_Methylation.RData")
