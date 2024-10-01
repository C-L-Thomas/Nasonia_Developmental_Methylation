library(methylKit)
library(grid)
library(readr)
library(ggplot2)
library(sqldf)

genome_annotation<-read.delim("~/Nvit_1.1_PSR_ALL_ANNOTATIONS_numbered_exons.txt", sep="\t",header = T)
colnames(genome_annotation) <- c("chr","feature","start","end","strand","geneID","feature_numb")

EvL<-read.delim("~/Embryo_vs_Larva_10%_sig.txt", sep="\t",header = T)

EvLinput <- EvL[,1:2]
colnames(EvLinput) <- c("chr2","start")

nrow(EvLinput)

head(EvLinput)
head(genome_annotation)

EvL_features <- sqldf("SELECT stage.chr2,
                stage.start,
                ga.chr,
                ga.feature,
                ga.start,
                ga.end,
                ga.geneID,
                ga.feature_numb
                FROM EvLinput AS stage
                LEFT JOIN genome_annotation AS ga
                ON stage.chr2 = ga.chr
                AND (stage.start >= ga.start AND stage.start <= ga.end)")

save.image(file = "testsave.RData")

write.table(EvL_features, file="Embryo_vs_Larva_Significant_Features.txt",
            sep="\t", row.names=F, quote=F)


head(EvL_features)

EvL_intron <- subset(EvL_features, EvL_features[,4] =="intron")
nrow(EvL_intron)
write.table(EvL_intron, file="Embryo_vs_Larva_Significant_Introns.txt",
            sep="\t", row.names=F, quote=F)

EvL_gene <- subset(EvL_features, EvL_features[,4] =="gene")
nrow(EvL_gene)
write.table(EvL_gene, file="Embryo_vs_Larva_Significant_Genes.txt",
            sep="\t", row.names=F, quote=F)

EvL_exon <- subset(EvL_features, EvL_features[,4] =="exon")
nrow(EvL_exon)
write.table(EvL_exon, file="Embryo_vs_Larva_Significant_Exons.txt",
            sep="\t", row.names=F, quote=F)

EvL_pseudogene <- subset(EvL_features, EvL_features[,4] =="pseudogene")
nrow(EvL_pseudogene)
write.table(EvL_pseudogene, file="Embryo_vs_Larva_Significant_pseudogenes.txt",
            sep="\t", row.names=F, quote=F)

EvL_intergenic <- subset(EvL_features, EvL_features[,4] =="intergenic")
nrow(EvL_intergenic)
write.table(EvL_intergenic, file="Embryo_vs_Larva_Significant_intergenics.txt",
            sep="\t", row.names=F, quote=F)

EvL_lnc_RNA <- subset(EvL_features, EvL_features[,4] =="lnc_RNA")
nrow(EvL_lnc_RNA)
write.table(EvL_lnc_RNA, file="Embryo_vs_Larva_Significant_lnc_RNAs.txt",
            sep="\t", row.names=F, quote=F)

EvL_five_prime_UTR <- subset(EvL_features, EvL_features[,4] =="five_prime_UTR")
nrow(EvL_five_prime_UTR)
write.table(EvL_five_prime_UTR, file="Embryo_vs_Larva_Significant_five_prime_UTRs.txt",
            sep="\t", row.names=F, quote=F)

EvL_three_prime_UTR <- subset(EvL_features, EvL_features[,4] =="three_prime_UTR")
nrow(EvL_three_prime_UTR)
write.table(EvL_three_prime_UTR, file="Embryo_vs_Larva_Significant_three_prime_UTRs.txt",
            sep="\t", row.names=F, quote=F)

EvL_promoter <- subset(EvL_features, EvL_features[,4] =="promoter")
nrow(EvL_promoter)
write.table(EvL_promoter, file="Embryo_vs_Larva_Significant_promoters.txt",
            sep="\t", row.names=F, quote=F)


LvP<-read.delim("~/Larva_vs_Prepupae_10%_sig.txt", sep="\t",header = T)
head(LvP)
LvP <- LvP[,1:2]
colnames(LvP) <- c("chr2","start")

head(LvP)
nrow(LvP)

larva_features <- sqldf("SELECT stage.chr2,
                stage.start,
                ga.chr,
                ga.feature,
                ga.start,
                ga.end,
                ga.geneID,
                ga.feature_numb
                FROM LvP AS stage
                LEFT JOIN genome_annotation AS ga
                ON stage.chr2 = ga.chr
                AND (stage.start >= ga.start AND stage.start <= ga.end)")

nrow(larva_features)

write.table(larva_features, file="Larva_vs_Prepupa_Features.txt",
            sep="\t", row.names=F, quote=F)

l_intron <- subset(larva_features, larva_features[,4] =="intron")
nrow(l_intron)
write.table(l_intron, file="Larva_vs_Prepupa_Introns.txt",
            sep="\t", row.names=F, quote=F)

l_gene <- subset(larva_features, larva_features[,4] =="gene")
nrow(l_gene)
write.table(l_gene, file="Larva_vs_Prepupa_Genes.txt",
            sep="\t", row.names=F, quote=F)

exon <- subset(larva_features, larva_features[,4] =="exon")
nrow(exon)
write.table(exon, file="Larva_vs_Prepupa_Exons.txt",
            sep="\t", row.names=F, quote=F)

pseudogene <- subset(larva_features, larva_features[,4] =="pseudogene")
nrow(pseudogene)
write.table(pseudogene, file="Larva_vs_Prepupa_pseudogenes.txt",
            sep="\t", row.names=F, quote=F)

intergenic <- subset(larva_features, larva_features[,4] =="intergenic")
nrow(intergenic)
write.table(intergenic, file="Larva_vs_Prepupa_intergenics.txt",
            sep="\t", row.names=F, quote=F)

lnc_RNA <- subset(larva_features, larva_features[,4] =="lnc_RNA")
nrow(lnc_RNA)
write.table(lnc_RNA, file="Larva_vs_Prepupa_lnc_RNAs.txt",
            sep="\t", row.names=F, quote=F)

five_prime_UTR <- subset(larva_features, larva_features[,4] =="five_prime_UTR")
nrow(five_prime_UTR)
write.table(five_prime_UTR, file="Larva_vs_Prepupa_five_prime_UTRs.txt",
            sep="\t", row.names=F, quote=F)

three_prime_UTR <- subset(larva_features, larva_features[,4] =="three_prime_UTR")
nrow(three_prime_UTR)
write.table(three_prime_UTR, file="Larva_vs_Prepupa_three_prime_UTRs.txt",
            sep="\t", row.names=F, quote=F)

promoter <- subset(larva_features, larva_features[,4] =="promoter")
nrow(promoter)
write.table(promoter, file="Larva_vs_Prepupa_promoters.txt",
            sep="\t", row.names=F, quote=F)

PvP<-read.delim("~/Prepupae_vs_Pupae_10%_sig.txt", sep="\t",header = T)

PvP <- PvP[,1:2]
colnames(PvP) <- c("chr2","start")

prepupae_features <- sqldf("SELECT stage.chr2,
                stage.start,
                ga.chr,
                ga.feature,
                ga.start,
                ga.end,
                ga.geneID,
                ga.feature_numb
                FROM PvP AS stage
                LEFT JOIN genome_annotation AS ga
                ON stage.chr2 = ga.chr
                AND (stage.start >= ga.start AND stage.start <= ga.end)")

nrow(prepupae_features)

write.table(prepupae_features, file="Prepupae_vs_Pupae_Features.txt",
            sep="\t", row.names=F, quote=F)

pr_intron <- subset(prepupae_features, prepupae_features[,4] =="intron")
nrow(pr_intron)
write.table(pr_intron, file="Prepupae_vs_Pupae_Introns.txt",
            sep="\t", row.names=F, quote=F)

pr_exon <- subset(prepupae_features, prepupae_features[,4] =="exon")
nrow(pr_exon)
write.table(pr_exon, file="Prepupae_vs_Pupae_Exons.txt",
            sep="\t", row.names=F, quote=F)

pr_gene <- subset(prepupae_features, prepupae_features[,4] =="gene")
nrow(pr_gene)
write.table(pr_gene, file="Prepupae_vs_Pupae_Genes.txt",
            sep="\t", row.names=F, quote=F)

pseudogene <- subset(prepupae_features, prepupae_features[,4] =="pseudogene")
nrow(pseudogene)
write.table(pseudogene, file="Prepupae_vs_Pupae_pseudogenes.txt",
            sep="\t", row.names=F, quote=F)

intergenic <- subset(prepupae_features, prepupae_features[,4] =="intergenic")
nrow(intergenic)
write.table(intergenic, file="Prepupae_vs_Pupae_intergenics.txt",
            sep="\t", row.names=F, quote=F)

lnc_RNA <- subset(prepupae_features, prepupae_features[,4] =="lnc_RNA")
nrow(lnc_RNA)
write.table(lnc_RNA, file="Prepupae_vs_Pupae_lnc_RNAs.txt",
            sep="\t", row.names=F, quote=F)

five_prime_UTR <- subset(prepupae_features, prepupae_features[,4] =="five_prime_UTR")
nrow(five_prime_UTR)
write.table(five_prime_UTR, file="Prepupae_vs_Pupae_five_prime_UTRs.txt",
            sep="\t", row.names=F, quote=F)

three_prime_UTR <- subset(prepupae_features, prepupae_features[,4] =="three_prime_UTR")
nrow(three_prime_UTR)
write.table(three_prime_UTR, file="Prepupae_vs_Pupae_three_prime_UTRs.txt",
            sep="\t", row.names=F, quote=F)

promoter <- subset(prepupae_features, prepupae_features[,4] =="promoter")
nrow(promoter)
write.table(promoter, file="Prepupae_vs_Pupae_promoters.txt",
            sep="\t", row.names=F, quote=F)

PvA<-read.delim("~/Pupae_vs_Adult_10%_sig.txt", sep="\t",header = T)

PvA <- PvA[,1:2]
colnames(PvA) <- c("chr2","start")

pupae_features <- sqldf("SELECT stage.chr2,
                stage.start,
                ga.chr,
                ga.feature,
                ga.start,
                ga.end,
                ga.geneID,
                ga.feature_numb
                FROM PvA AS stage
                LEFT JOIN genome_annotation AS ga
                ON stage.chr2 = ga.chr
                AND (stage.start >= ga.start AND stage.start <= ga.end)")

nrow(pupae_features)

write.table(pupae_features, file="Pupae_vs_Adult_Features.txt",
            sep="\t", row.names=F, quote=F)

pu_intron <- subset(pupae_features, pupae_features[,4] =="intron")
nrow(pu_intron)
write.table(pu_intron, file="Pupae_vs_Adult_Introns.txt",
            sep="\t", row.names=F, quote=F)

pu_gene <- subset(pupae_features, pupae_features[,4] =="gene")
nrow(pu_gene)
write.table(pu_gene, file="Pupae_vs_Adult_Genes.txt",
            sep="\t", row.names=F, quote=F)

pu_exon <- subset(pupae_features, pupae_features[,4] =="exon")
nrow(pu_exon)
write.table(pu_exon, file="Pupae_vs_Adult_Exons.txt",
            sep="\t", row.names=F, quote=F)

pseudogene <- subset(pupae_features, pupae_features[,4] =="pseudogene")
nrow(pseudogene)
write.table(pseudogene, file="Pupae_vs_Adult_pseudogenes.txt",
            sep="\t", row.names=F, quote=F)

intergenic <- subset(pupae_features, pupae_features[,4] =="intergenic")
nrow(intergenic)
write.table(intergenic, file="Pupae_vs_Adult_intergenics.txt",
            sep="\t", row.names=F, quote=F)

lnc_RNA <- subset(pupae_features, pupae_features[,4] =="lnc_RNA")
nrow(lnc_RNA)
write.table(lnc_RNA, file="Pupae_vs_Adult_lnc_RNAs.txt",
            sep="\t", row.names=F, quote=F)

five_prime_UTR <- subset(pupae_features, pupae_features[,4] =="five_prime_UTR")
nrow(five_prime_UTR)
write.table(five_prime_UTR, file="Pupae_vs_Adult_five_prime_UTRs.txt",
            sep="\t", row.names=F, quote=F)

three_prime_UTR <- subset(pupae_features, pupae_features[,4] =="three_prime_UTR")
nrow(three_prime_UTR)
write.table(three_prime_UTR, file="Pupae_vs_Adult_three_prime_UTRs.txt",
            sep="\t", row.names=F, quote=F)

promoter <- subset(pupae_features, pupae_features[,4] =="promoter")
nrow(promoter)
write.table(promoter, file="Pupae_vs_Adult_promoters.txt",
            sep="\t", row.names=F, quote=F)
