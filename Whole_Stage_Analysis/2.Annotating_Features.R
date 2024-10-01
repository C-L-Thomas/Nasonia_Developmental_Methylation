library(methylKit)
library(grid)
library(readr)
library(ggplot2)
library(sqldf)

load(file = "Stage_Level_Methylation.RData")

genome_annotation<-read.delim("~/Nvit_1.1_PSR_ALL_ANNOTATIONS_numbered_exons.txt", sep="\t",header = T) #A file located in Nasonia Genome folder
colnames(genome_annotation) <- c("chr","feature","start","end","strand","geneID","feature_numb")

head(emb_meth_sites)

embryoinput <- emb_meth_sites[,1:2]

nrow(embryoinput)

embryo_features <- sqldf("SELECT stage.chr,
                stage.start,
                ga.chr,
                ga.feature,
                ga.start,
                ga.end,
                ga.geneID,
                ga.feature_numb
                FROM embryoinput AS stage
                LEFT JOIN genome_annotation AS ga
                ON stage.chr = ga.chr
                AND (stage.start >= ga.start AND stage.start <= ga.end)")

write.table(embryo_features, file="embryo_Features.txt",
            sep="\t", row.names=F, quote=F)


e_intron <- subset(embryo_features, embryo_features[,4] =="intron")
nrow(e_intron)
write.table(e_intron, file="embryo_Introns.txt",
            sep="\t", row.names=F, quote=F)

e_gene <- subset(embryo_features, embryo_features[,4] =="gene")
nrow(e_gene)
write.table(e_gene, file="embryo_Genes.txt",
            sep="\t", row.names=F, quote=F)

e_exon <- subset(embryo_features, embryo_features[,4] =="exon")
nrow(e_exon)
write.table(e_exon, file="embryo_Exons.txt",
            sep="\t", row.names=F, quote=F)

pseudogene <- subset(embryo_features, embryo_features[,4] =="pseudogene")
nrow(pseudogene)
write.table(pseudogene, file="embryo_pseudogenes.txt",
            sep="\t", row.names=F, quote=F)

intergenic <- subset(embryo_features, embryo_features[,4] =="intergenic")
nrow(intergenic)
write.table(intergenic, file="embryo_intergenics.txt",
            sep="\t", row.names=F, quote=F)

lnc_RNA <- subset(embryo_features, embryo_features[,4] =="lnc_RNA")
nrow(lnc_RNA)
write.table(lnc_RNA, file="embryo_lnc_RNAs.txt",
            sep="\t", row.names=F, quote=F)

five_prime_UTR <- subset(embryo_features, embryo_features[,4] =="five_prime_UTR")
nrow(five_prime_UTR)
write.table(five_prime_UTR, file="embryo_five_prime_UTRs.txt",
            sep="\t", row.names=F, quote=F)

three_prime_UTR <- subset(embryo_features, embryo_features[,4] =="three_prime_UTR")
nrow(three_prime_UTR)
write.table(three_prime_UTR, file="embryo_three_prime_UTRs.txt",
            sep="\t", row.names=F, quote=F)

promoter <- subset(embryo_features, embryo_features[,4] =="promoter")
nrow(promoter)
write.table(promoter, file="embryo_promoters.txt",
            sep="\t", row.names=F, quote=F)

larvainput <- lar_meth_sites[,1:2]

nrow(larvainput)

larva_features <- sqldf("SELECT stage.chr,
                stage.start,
                ga.chr,
                ga.feature,
                ga.start,
                ga.end,
                ga.geneID,
                ga.feature_numb
                FROM larvainput AS stage
                LEFT JOIN genome_annotation AS ga
                ON stage.chr = ga.chr
                AND (stage.start >= ga.start AND stage.start <= ga.end)")

nrow(larva_features)

write.table(larva_features, file="larva_Features.txt",
            sep="\t", row.names=F, quote=F)

l_intron <- subset(larva_features, larva_features[,4] =="intron")
nrow(l_intron)
write.table(l_intron, file="larva_Introns.txt",
            sep="\t", row.names=F, quote=F)

l_gene <- subset(larva_features, larva_features[,4] =="gene")
nrow(l_gene)
write.table(l_gene, file="larva_Genes.txt",
            sep="\t", row.names=F, quote=F)

exon <- subset(larva_features, larva_features[,4] =="exon")
nrow(exon)
write.table(exon, file="larva_Exons.txt",
            sep="\t", row.names=F, quote=F)

pseudogene <- subset(larva_features, larva_features[,4] =="pseudogene")
nrow(pseudogene)
write.table(pseudogene, file="larva_pseudogenes.txt",
            sep="\t", row.names=F, quote=F)

intergenic <- subset(larva_features, larva_features[,4] =="intergenic")
nrow(intergenic)
write.table(intergenic, file="larva_intergenics.txt",
            sep="\t", row.names=F, quote=F)

lnc_RNA <- subset(larva_features, larva_features[,4] =="lnc_RNA")
nrow(lnc_RNA)
write.table(lnc_RNA, file="larva_lnc_RNAs.txt",
            sep="\t", row.names=F, quote=F)

five_prime_UTR <- subset(larva_features, larva_features[,4] =="five_prime_UTR")
nrow(five_prime_UTR)
write.table(five_prime_UTR, file="larva_five_prime_UTRs.txt",
            sep="\t", row.names=F, quote=F)

three_prime_UTR <- subset(larva_features, larva_features[,4] =="three_prime_UTR")
nrow(three_prime_UTR)
write.table(three_prime_UTR, file="larva_three_prime_UTRs.txt",
            sep="\t", row.names=F, quote=F)

promoter <- subset(larva_features, larva_features[,4] =="promoter")
nrow(promoter)
write.table(promoter, file="larva_promoters.txt",
            sep="\t", row.names=F, quote=F)

prepupaeinput <- pre_meth_sites[,1:2]

nrow(prepupaeinput)

prepupae_features <- sqldf("SELECT stage.chr,
                stage.start,
                ga.chr,
                ga.feature,
                ga.start,
                ga.end,
                ga.geneID,
                ga.feature_numb
                FROM prepupaeinput AS stage
                LEFT JOIN genome_annotation AS ga
                ON stage.chr = ga.chr
                AND (stage.start >= ga.start AND stage.start <= ga.end)")

nrow(prepupae_features)

write.table(prepupae_features, file="prepupae_Features.txt",
            sep="\t", row.names=F, quote=F)

pr_intron <- subset(prepupae_features, prepupae_features[,4] =="intron")
nrow(pr_intron)
write.table(pr_intron, file="prepupae_Introns.txt",
            sep="\t", row.names=F, quote=F)

pr_exon <- subset(prepupae_features, prepupae_features[,4] =="exon")
nrow(pr_exon)
write.table(pr_exon, file="prepupae_Exons.txt",
            sep="\t", row.names=F, quote=F)

pr_gene <- subset(prepupae_features, prepupae_features[,4] =="gene")
nrow(pr_gene)
write.table(pr_gene, file="prepupae_Genes.txt",
            sep="\t", row.names=F, quote=F)

pseudogene <- subset(prepupae_features, prepupae_features[,4] =="pseudogene")
nrow(pseudogene)
write.table(pseudogene, file="prepupae_pseudogenes.txt",
            sep="\t", row.names=F, quote=F)

intergenic <- subset(prepupae_features, prepupae_features[,4] =="intergenic")
nrow(intergenic)
write.table(intergenic, file="prepupae_intergenics.txt",
            sep="\t", row.names=F, quote=F)

lnc_RNA <- subset(prepupae_features, prepupae_features[,4] =="lnc_RNA")
nrow(lnc_RNA)
write.table(lnc_RNA, file="prepupae_lnc_RNAs.txt",
            sep="\t", row.names=F, quote=F)

five_prime_UTR <- subset(prepupae_features, prepupae_features[,4] =="five_prime_UTR")
nrow(five_prime_UTR)
write.table(five_prime_UTR, file="prepupae_five_prime_UTRs.txt",
            sep="\t", row.names=F, quote=F)

three_prime_UTR <- subset(prepupae_features, prepupae_features[,4] =="three_prime_UTR")
nrow(three_prime_UTR)
write.table(three_prime_UTR, file="prepupae_three_prime_UTRs.txt",
            sep="\t", row.names=F, quote=F)

promoter <- subset(prepupae_features, prepupae_features[,4] =="promoter")
nrow(promoter)
write.table(promoter, file="prepupae_promoters.txt",
            sep="\t", row.names=F, quote=F)

pupaeinput <- pup_meth_sites[,1:2]

nrow(pupaeinput)

pupae_features <- sqldf("SELECT stage.chr,
                stage.start,
                ga.chr,
                ga.feature,
                ga.start,
                ga.end,
                ga.geneID,
                ga.feature_numb
                FROM pupaeinput AS stage
                LEFT JOIN genome_annotation AS ga
                ON stage.chr = ga.chr
                AND (stage.start >= ga.start AND stage.start <= ga.end)")

nrow(pupae_features)

write.table(pupae_features, file="pupae_Features.txt",
            sep="\t", row.names=F, quote=F)

pu_intron <- subset(pupae_features, pupae_features[,4] =="intron")
nrow(pu_intron)
write.table(pu_intron, file="pupae_Introns.txt",
            sep="\t", row.names=F, quote=F)

pu_gene <- subset(pupae_features, pupae_features[,4] =="gene")
nrow(pu_gene)
write.table(pu_gene, file="pupae_Genes.txt",
            sep="\t", row.names=F, quote=F)

pu_exon <- subset(pupae_features, pupae_features[,4] =="exon")
nrow(pu_exon)
write.table(pu_exon, file="pupae_Exons.txt",
            sep="\t", row.names=F, quote=F)

pseudogene <- subset(pupae_features, pupae_features[,4] =="pseudogene")
nrow(pseudogene)
write.table(pseudogene, file="pupae_pseudogenes.txt",
            sep="\t", row.names=F, quote=F)

intergenic <- subset(pupae_features, pupae_features[,4] =="intergenic")
nrow(intergenic)
write.table(intergenic, file="pupae_intergenics.txt",
            sep="\t", row.names=F, quote=F)

lnc_RNA <- subset(pupae_features, pupae_features[,4] =="lnc_RNA")
nrow(lnc_RNA)
write.table(lnc_RNA, file="pupae_lnc_RNAs.txt",
            sep="\t", row.names=F, quote=F)

five_prime_UTR <- subset(pupae_features, pupae_features[,4] =="five_prime_UTR")
nrow(five_prime_UTR)
write.table(five_prime_UTR, file="pupae_five_prime_UTRs.txt",
            sep="\t", row.names=F, quote=F)

three_prime_UTR <- subset(pupae_features, pupae_features[,4] =="three_prime_UTR")
nrow(three_prime_UTR)
write.table(three_prime_UTR, file="pupae_three_prime_UTRs.txt",
            sep="\t", row.names=F, quote=F)

promoter <- subset(pupae_features, pupae_features[,4] =="promoter")
nrow(promoter)
write.table(promoter, file="pupae_promoters.txt",
            sep="\t", row.names=F, quote=F)

adultinput <- adu_meth_sites[,1:2]

nrow(adultinput)

adult_features <- sqldf("SELECT stage.chr,
                stage.start,
                ga.chr,
                ga.feature,
                ga.start,
                ga.end,
                ga.geneID,
                ga.feature_numb
                FROM adultinput AS stage
                LEFT JOIN genome_annotation AS ga
                ON stage.chr = ga.chr
                AND (stage.start >= ga.start AND stage.start <= ga.end)")

nrow(adult_features)

write.table(adult_features, file="adult_Features.txt",
            sep="\t", row.names=F, quote=F)

a_intron <- subset(adult_features, adult_features[,4] =="intron")
nrow(a_intron)
write.table(a_intron, file="adult_Introns.txt",
            sep="\t", row.names=F, quote=F)

a_exon <- subset(adult_features, adult_features[,4] =="exon")
nrow(a_exon)
write.table(a_exon, file="adult_Exons.txt",
            sep="\t", row.names=F, quote=F)

a_gene <- subset(adult_features, adult_features[,4] =="gene")
nrow(a_gene)
write.table(a_gene, file="adult_Genes.txt",
            sep="\t", row.names=F, quote=F)

pseudogene <- subset(adult_features, adult_features[,4] =="pseudogene")
nrow(pseudogene)
write.table(pseudogene, file="adult_pseudogenes.txt",
            sep="\t", row.names=F, quote=F)

intergenic <- subset(adult_features, adult_features[,4] =="intergenic")
nrow(intergenic)
write.table(intergenic, file="adult_intergenics.txt",
            sep="\t", row.names=F, quote=F)

lnc_RNA <- subset(adult_features, adult_features[,4] =="lnc_RNA")
nrow(lnc_RNA)
write.table(lnc_RNA, file="adult_lnc_RNAs.txt",
            sep="\t", row.names=F, quote=F)

five_prime_UTR <- subset(adult_features, adult_features[,4] =="five_prime_UTR")
nrow(five_prime_UTR)
write.table(five_prime_UTR, file="adult_five_prime_UTRs.txt",
            sep="\t", row.names=F, quote=F)

three_prime_UTR <- subset(adult_features, adult_features[,4] =="three_prime_UTR")
nrow(three_prime_UTR)
write.table(three_prime_UTR, file="adult_three_prime_UTRs.txt",
            sep="\t", row.names=F, quote=F)

promoter <- subset(adult_features, adult_features[,4] =="promoter")
nrow(promoter)
write.table(promoter, file="adult_promoters.txt",
            sep="\t", row.names=F, quote=F)
