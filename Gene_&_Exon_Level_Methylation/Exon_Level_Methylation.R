library(methylKit)
library(grid)
library(readr)
library(ggplot2)
library(sqldf)
library(dplyr)
library(doBy)
library(genomation)

file.list <- list("Embryo1.cov","Embryo2.cov","Larva1.cov", "Larva2.cov","Larva3.cov","Prepupa1.cov","Prepupa2.cov","Prepupa3.cov",
                  "Pupa1.cov","Pupa2.cov","Pupa3.cov","Adult1.cov","Adult2.cov","Adult3.cov")

treat_conditions <- c(0,0,1,1,1,2,2,2,3,3,3,4,4,4)

raw_data <- methRead(file.list,
                     sample.id = list("EMB1","EMB2","LAR1","LAR2","LAR3","PRE1","PRE2","PRE3","PUP1","PUP2","PUP3","ADU1","ADU2","ADU3"),
                     treatment = treat_conditions,
                     assembly="Nvit_PSR_1.1",
                     dbtype = NA,
                     pipeline = "bismarkCoverage",
                     header = T, sep = "\t", mincov=1,
                     context = 'CpG')

filtered_data <- filterByCoverage(raw_data,lo.count=10,lo.perc=NULL,
                                  hi.count=NULL,hi.perc=99.9)

meth_all_data <- unite(filtered_data, destrand=FALSE)

#Input Gene Files
Exon_Regions <- readGeneric("~/Nvit_psr_1.1/Beds/exons.bed", header = FALSE)
Exon_coverage <- regionCounts(meth_all_data,Gene_Regions)
Exons <- getData(Exon_coverage)

head(gene_coverage)
library(tidyverse)
Exons$chr <- str_replace(Exons$chr, "_",".")

Exon_annotation <- read.delim("~Nvit_psr_1.1/Beds/exon.bed", header = FALSE)
colnames(Exon_annotation) <- c("chrom","type","Start","End","strand","geneID","number")
head(Exon_annotation)


output <- sqldf("SELECT out.chr,
                out.start,
                out.end,
                out.strand,
                out.coverage1,
                out.numCs1,
                out.numTs1,
                out.coverage2,
                out.numCs2,
                out.numTs2,
                out.coverage3,
                out.numCs3,
                out.numTs3,
                out.coverage4,
                out.numCs4,
                out.numTs4,
                out.coverage5,
                out.numCs5,
                out.numTs5,
                out.coverage6,
                out.numCs6,
                out.numTs6,
                out.coverage7,
                out.numCs7,
                out.numTs7,
		out.coverage8,
                out.numCs8,
                out.numTs8,
                out.coverage9,
                out.numCs9,
                out.numTs9,
                out.coverage10,
                out.numCs10,
                out.numTs10,
                out.coverage11,
                out.numCs11,
                out.numTs11,
                out.coverage12,
                out.numCs12,
                out.numTs12,
                out.coverage13,
                out.numCs13,
                out.numTs13,
                out.coverage14,
                out.numCs14,
                out.numTs14,
                ga.chrom,
                ga.Start,
                ga.End,
                ga.geneID,
		ga.number
                FROM Exons AS out
                LEFT JOIN Exon_annotation AS ga
                ON out.chr = ga.chrom
                AND out.start = ga.Start
		AND out.end = ga.End")

write.table(output, "Total_Exons.txt", quote=F, row.names = F, sep = '\t' )

#______________________________________________________________________________
#Reorganise All Analyses

Embryo_vs_Larva=reorganize(Exon_coverage,sample.ids=c("EMB1","EMB2","LAR1","LAR2","LAR3"),
                                treatment=c(0,0,1,1,1))

Larva_vs_Prepupae=reorganize(Exon_coverage,sample.ids=c("LAR1","LAR2","LAR3","PRE1","PRE2","PRE3"),
                                treatment=c(1,1,1,2,2,2))

Prepupae_vs_Pupae=reorganize(Exon_coverage,sample.ids=c("PRE1","PRE2","PRE3","PUP1","PUP2","PUP3"),
                                treatment=c(2,2,2,3,3,3))

Pupae_vs_Adult=reorganize(Exon_coverage,sample.ids=c("PUP1","PUP2","PUP3","ADU1","ADU2","ADU3"),
                                treatment=c(3,3,3,4,4,4))


#______________________________________________________________________________
#Differential

library(tidyverse)

diff_Embryo_vs_Larva <- calculateDiffMeth(Embryo_vs_Larva, adjust = c("BH"))
nrow(diff_Embryo_vs_Larva)
head(diff_Embryo_vs_Larva)
write.table(diff_Embryo_vs_Larva, file="Embryo_vs_Larva_Complete_Exon_Result.txt", quote=F, row.names = F, sep = '\t')
diff_Embryo_vs_Larva$chr <- str_replace(diff_Embryo_vs_Larva$chr, "_",".")

EvL <- getData(diff_Embryo_vs_Larva)
colnames(EvL) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")
head(EvL)

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chrom,
                ga.Start,
                ga.End,
                ga.geneID,
		ga.number
                FROM EvL AS stage
                LEFT JOIN Exon_annotation AS ga
                ON stage.chr = ga.chrom
                AND (stage.start = ga.start AND stage.end = ga.end)")

head(output)

cols <- c(1,2,3,4,5,6,10,11)
total_diff_Embryo_vs_Larva <- output[,c(cols)]

write.table(total_diff_Embryo_vs_Larva, file="Embryo_vs_Larva_Complete_Exon_Result_w_gene.txt", quote=F, row.names = F, sep = '\t')

diff_Embryo_vs_Larva_10 <- getMethylDiff(diff_Embryo_vs_Larva, difference=10, qvalue=0.0125)
write.table(diff_Embryo_vs_Larva_10, file="Embryo_vs_Larva_Exon_10%_sig.txt", quote=F, row.names = F, sep = '\t')
EvL <- getData(diff_Embryo_vs_Larva_10)
nrow(EvL)
colnames(EvL) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chrom,
                ga.Start,
                ga.End,
                ga.geneID,
                ga.number
                FROM EvL AS stage
                LEFT JOIN Exon_annotation AS ga
                ON stage.chr = ga.chrom
                AND (stage.start = ga.start AND stage.end = ga.end)")

head(output)

cols <- c(1,2,3,4,5,6,10,11)
total_diff_Embryo_vs_Larva_10 <- output[,c(cols)]
write.table(diff_Embryo_vs_Larva_10, file="Embryo_vs_Larva_Exon_10%_sig.txt", quote=F, row.names = F, sep = '\t')
write.table(total_diff_Embryo_vs_Larva_10, file="Embryo_vs_Larva_Exon_10%_sig_w_gene.txt", quote=F, row.names = F, sep = '\t')

diff_Embryo_vs_Larva_hyper=getMethylDiff(diff_Embryo_vs_Larva,difference=10,qvalue=0.0125,type="hyper")
EvL <- getData(diff_Embryo_vs_Larva_hyper)
head(EvL)

colnames(EvL) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chrom,
                ga.Start,
                ga.End,
                ga.geneID,
                ga.number
                FROM EvL AS stage
                LEFT JOIN Exon_annotation AS ga
                ON stage.chr = ga.chrom
                AND (stage.start = ga.start AND stage.end = ga.end)")

cols <- c(1,2,3,4,5,6,10,11)
hyper_diff_Embryo_vs_Larva <- output[,c(cols)]
write.table(diff_Embryo_vs_Larva_hyper, file="Embryo_vs_Larva_Larva_Hypermethylated_Exons_10%_sig.txt", quote=F, row.names = T, sep = '\t')
write.table(hyper_diff_Embryo_vs_Larva, file="Embryo_vs_Larva_Larva_Hypermethylated_Exons_10%_sig_w_gene.txt", quote=F, row.names = T, sep = '\t')

diff_Embryo_vs_Larva_hypo=getMethylDiff(diff_Embryo_vs_Larva,difference=10,qvalue=0.0125,type="hypo")
EvL <- getData(diff_Embryo_vs_Larva_hypo)
head(EvL)

colnames(EvL) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chrom,
                ga.Start,
                ga.End,
                ga.geneID,
                ga.number
                FROM EvL AS stage
                LEFT JOIN Exon_annotation AS ga
                ON stage.chr = ga.chrom
                AND (stage.start = ga.start AND stage.end = ga.end)")

cols <- c(1,2,3,4,5,6,10,11)
hypo_diff_Embryo_vs_Larva <- output[,c(cols)]
write.table(diff_Embryo_vs_Larva_hypo, file="Embryo_vs_Larva_Larva_Hypomethylated_Exons_10%_sig_w_gene.txt", quote=F, row.names = T, sep = '\t')
write.table(hypo_diff_Embryo_vs_Larva, file="Embryo_vs_Larva_Larva_Hypomethylated_Exons_10%_sig.txt", quote=F, row.names = T, sep = '\t')

################################################################################################################

diff_Larva_vs_Prepupae <- calculateDiffMeth(Larva_vs_Prepupae, adjust = c("BH"))
nrow(diff_Larva_vs_Prepupae)
head(diff_Larva_vs_Prepupae)
diff_Larva_vs_Prepupae$chr <- str_replace(diff_Larva_vs_Prepupae$chr, "_",".")

LvP <- getData(diff_Larva_vs_Prepupae)
colnames(LvP) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")
head(LvP)

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chrom,
                ga.Start,
                ga.End,
                ga.geneID,
                ga.number
                FROM LvP AS stage
                LEFT JOIN Exon_annotation AS ga
                ON stage.chr = ga.chrom
                AND (stage.start = ga.start AND stage.end = ga.end)")

head(output)

cols <- c(1,2,3,4,5,6,10,11)
total_diff_Larva_vs_Prepupae <- output[,c(cols)]
write.table(diff_Larva_vs_Prepupae, file="Larva_vs_Prepupae_Complete_Exon_Result.txt", quote=F, row.names = F, sep = '\t')
write.table(total_diff_Larva_vs_Prepupae, file="Larva_vs_Prepupae_Complete_Exon_Result_w_gene.txt", quote=F, row.names = F, sep = '\t')

diff_Larva_vs_Prepupae_10 <- getMethylDiff(diff_Larva_vs_Prepupae, difference=10, qvalue=0.0125)

LvP <- getData(diff_Larva_vs_Prepupae_10)
nrow(LvP)
colnames(LvP) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chrom,
                ga.Start,
                ga.End,
                ga.geneID,
                ga.number
                FROM LvP AS stage
                LEFT JOIN Exon_annotation AS ga
                ON stage.chr = ga.chrom
                AND (stage.start = ga.start AND stage.end = ga.end)")

head(output)

cols <- c(1,2,3,4,5,6,10,11)
total_diff_Larva_vs_Prepupae_10 <- output[,c(cols)]
write.table(diff_Larva_vs_Prepupae_10, file="Larva_vs_Prepupae_Exon_10%_sig.txt", quote=F, row.names = F, sep = '\t')
write.table(total_diff_Larva_vs_Prepupae_10, file="Larva_vs_Prepupae_Exon_10%_sig_w_gene.txt", quote=F, row.names = F, sep = '\t')

diff_Larva_vs_Prepupae_hyper=getMethylDiff(diff_Larva_vs_Prepupae,difference=10,qvalue=0.0125,type="hyper")
LvP <- getData(diff_Larva_vs_Prepupae_hyper)
head(LvP)

colnames(LvP) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chrom,
                ga.Start,
                ga.End,
                ga.geneID,
                ga.number
                FROM LvP AS stage
                LEFT JOIN Exon_annotation AS ga
                ON stage.chr = ga.chrom
                AND (stage.start = ga.start AND stage.end = ga.end)")

cols <- c(1,2,3,4,5,6,10,11)
hyper_diff_Larva_vs_Prepupae <- output[,c(cols)]
write.table(diff_Larva_vs_Prepupae_hyper, file="Larva_vs_Prepupae_Prepupae_Hypermethylated_Exons_10%_sig.txt", quote=F, row.names = T, sep = '\t')
write.table(hyper_diff_Larva_vs_Prepupae, file="Larva_vs_Prepupae_Prepupae_Hypermethylated_Exons_10%_sig_w_gene.txt", quote=F, row.names = T, sep = '\t')

diff_Larva_vs_Prepupae_hypo=getMethylDiff(diff_Larva_vs_Prepupae,difference=10,qvalue=0.0125,type="hypo")
LvP <- getData(diff_Larva_vs_Prepupae_hypo)
head(LvP)

colnames(LvP) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chrom,
                ga.Start,
                ga.End,
                ga.geneID,
                ga.number
                FROM LvP AS stage
                LEFT JOIN Exon_annotation AS ga
                ON stage.chr = ga.chrom
                AND (stage.start = ga.start AND stage.end = ga.end)")

cols <- c(1,2,3,4,5,6,10,11)
hypo_diff_Larva_vs_Prepupae <- output[,c(cols)]
write.table(diff_Larva_vs_Prepupae_hypo, file="Larva_vs_Prepupae_Prepupae_Hypomethylated_Exons_10%_sig.txt", quote=F, row.names = T, sep = '\t')
write.table(hypo_diff_Larva_vs_Prepupae, file="Larva_vs_Prepupae_Prepupae_Hypomethylated_Exons_10%_sig_w_gene.txt", quote=F, row.names = T, sep = '\t')

###############################################################################################################################

diff_Prepupae_vs_Pupae <- calculateDiffMeth(Prepupae_vs_Pupae, adjust = c("BH"))
nrow(diff_Prepupae_vs_Pupae)
head(diff_Prepupae_vs_Pupae)
diff_Prepupae_vs_Pupae$chr <- str_replace(diff_Prepupae_vs_Pupae$chr, "_",".")

PvP <- getData(diff_Prepupae_vs_Pupae)
colnames(PvP) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")
head(PvP)

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chrom,
                ga.Start,
                ga.End,
                ga.geneID,
                ga.number
                FROM PvP AS stage
                LEFT JOIN Exon_annotation AS ga
                ON stage.chr = ga.chrom
                AND (stage.start = ga.start AND stage.end = ga.end)")

head(output)

cols <- c(1,2,3,4,5,6,10,11)
total_diff_Prepupae_vs_Pupae <- output[,c(cols)]
write.table(diff_Prepupae_vs_Pupae, file="Prepupae_vs_Pupae_Complete_Exon_Result.txt", quote=F, row.names = F, sep = '\t')
write.table(total_diff_Prepupae_vs_Pupae, file="Prepupae_vs_Pupae_Complete_Exon_Result_w_gene.txt", quote=F, row.names = F, sep = '\t')

diff_Prepupae_vs_Pupae_10 <- getMethylDiff(diff_Prepupae_vs_Pupae, difference=10, qvalue=0.0125)

PvP <- getData(diff_Prepupae_vs_Pupae_10)
nrow(PvP)
colnames(PvP) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chrom,
                ga.Start,
                ga.End,
                ga.geneID,
                ga.number
                FROM PvP AS stage
                LEFT JOIN Exon_annotation AS ga
                ON stage.chr = ga.chrom
                AND (stage.start = ga.start AND stage.end = ga.end)")

head(output)

cols <- c(1,2,3,4,5,6,10,11)
total_diff_Prepupae_vs_Pupae_10 <- output[,c(cols)]
write.table(diff_Prepupae_vs_Pupae_10, file="Prepupae_vs_Pupae_Exon_10%_sig.txt", quote=F, row.names = F, sep = '\t')
write.table(total_diff_Prepupae_vs_Pupae_10, file="Prepupae_vs_Pupae_Exon_10%_sig_w_gene.txt", quote=F, row.names = F, sep = '\t')

diff_Prepupae_vs_Pupae_hyper=getMethylDiff(diff_Prepupae_vs_Pupae,difference=10,qvalue=0.0125,type="hyper")
PvP <- getData(diff_Prepupae_vs_Pupae_hyper)
head(PvP)

colnames(PvP) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chrom,
                ga.Start,
                ga.End,
                ga.geneID,
                ga.number
                FROM PvP AS stage
                LEFT JOIN Exon_annotation AS ga
                ON stage.chr = ga.chrom
                AND (stage.start = ga.start AND stage.end = ga.end)")

cols <- c(1,2,3,4,5,6,10,11)
hyper_diff_Prepupae_vs_Pupae <- output[,c(cols)]
write.table(diff_Prepupae_vs_Pupae_hyper, file="Prepupae_vs_Pupae_Pupae_Hypermethylated_Exons_10%_sig.txt", quote=F, row.names = T, sep = '\t')
write.table(hyper_diff_Prepupae_vs_Pupae, file="Prepupae_vs_Pupae_Pupae_Hypermethylated_Exons_10%_sig_w_gene.txt", quote=F, row.names = T, sep = '\t')

diff_Prepupae_vs_Pupae_hypo=getMethylDiff(diff_Prepupae_vs_Pupae,difference=10,qvalue=0.0125,type="hypo")
PvP <- getData(diff_Prepupae_vs_Pupae_hypo)
head(PvP)

colnames(PvP) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chrom,
                ga.Start,
                ga.End,
                ga.geneID,
                ga.number
                FROM PvP AS stage
                LEFT JOIN Exon_annotation AS ga
                ON stage.chr = ga.chrom
                AND (stage.start = ga.start AND stage.end = ga.end)")

cols <- c(1,2,3,4,5,6,10,11)
hypo_diff_Prepupae_vs_Pupae <- output[,c(cols)]
write.table(diff_Prepupae_vs_Pupae_hypo, file="Prepupae_vs_Pupae_Pupae_Hypomethylated_Exons_10%_sig.txt", quote=F, row.names = T, sep = '\t')
write.table(hypo_diff_Prepupae_vs_Pupae, file="Prepupae_vs_Pupae_Pupae_Hypomethylated_Exons_10%_sig_w_gene.txt", quote=F, row.names = T, sep = '\t')

#############################################################################################################################################

diff_Pupae_vs_Adult <- calculateDiffMeth(Pupae_vs_Adult, adjust = c("BH"))
nrow(diff_Pupae_vs_Adult)
head(diff_Pupae_vs_Adult)
diff_Pupae_vs_Adult$chr <- str_replace(diff_Pupae_vs_Adult$chr, "_",".")

PvA <- getData(diff_Pupae_vs_Adult)
colnames(PvA) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")
head(PvA)

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chrom,
                ga.Start,
                ga.End,
                ga.geneID,
                ga.number
                FROM PvA AS stage
                LEFT JOIN Exon_annotation AS ga
                ON stage.chr = ga.chrom
                AND (stage.start = ga.start AND stage.end = ga.end)")

head(output)

cols <- c(1,2,3,4,5,6,10,11)
total_diff_Pupae_vs_Adult <- output[,c(cols)]
write.table(diff_Pupae_vs_Adult, file="Pupae_vs_Adult_Complete_Exon_Result.txt", quote=F, row.names = F, sep = '\t')
write.table(total_diff_Pupae_vs_Adult, file="Pupae_vs_Adult_Complete_Exon_Result_w_gene.txt", quote=F, row.names = F, sep = '\t')

diff_Pupae_vs_Adult_10 <- getMethylDiff(diff_Pupae_vs_Adult, difference=10, qvalue=0.0125)

PvA <- getData(diff_Pupae_vs_Adult_10)
nrow(PvA)
colnames(PvA) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chrom,
                ga.Start,
                ga.End,
                ga.geneID,
                ga.number
                FROM PvA AS stage
                LEFT JOIN Exon_annotation AS ga
                ON stage.chr = ga.chrom
                AND (stage.start = ga.start AND stage.end = ga.end)")

head(output)

cols <- c(1,2,3,4,5,6,10,11)
total_diff_Pupae_vs_Adult_10 <- output[,c(cols)]
write.table(diff_Pupae_vs_Adult_10, file="Pupae_vs_Adult_Exon_10%_sig.txt", quote=F, row.names = F, sep = '\t')
write.table(total_diff_Pupae_vs_Adult_10, file="Pupae_vs_Adult_Exon_10%_sig_w_gene.txt", quote=F, row.names = F, sep = '\t')

diff_Pupae_vs_Adult_hyper=getMethylDiff(diff_Pupae_vs_Adult,difference=10,qvalue=0.0125,type="hyper")
PvA <- getData(diff_Pupae_vs_Adult_hyper)
head(PvA)

colnames(PvA) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chrom,
                ga.Start,
                ga.End,
                ga.geneID,
                ga.number
                FROM PvA AS stage
                LEFT JOIN Exon_annotation AS ga
                ON stage.chr = ga.chrom
                AND (stage.start = ga.start AND stage.end = ga.end)")

cols <- c(1,2,3,4,5,6,10,11)
diff_Pupae_vs_Adult_hyper <- output[,c(cols)]
write.table(hyper_diff_Pupae_vs_Adult, file="Pupae_vs_Adult_Adult_Hypermethylated_Exons_10%_sig.txt", quote=F, row.names = T, sep = '\t')
write.table(hyper_diff_Pupae_vs_Adult, file="Pupae_vs_Adult_Adult_Hypermethylated_Exons_10%_sig_w_gene.txt", quote=F, row.names = T, sep = '\t')

diff_Pupae_vs_Adult_hypo=getMethylDiff(diff_Pupae_vs_Adult,difference=10,qvalue=0.0125,type="hypo")
PvA <- getData(diff_Pupae_vs_Adult_hypo)
head(PvA)

colnames(PvA) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chrom,
                ga.Start,
                ga.End,
                ga.geneID,
                ga.number
                FROM PvA AS stage
                LEFT JOIN Exon_annotation AS ga
                ON stage.chr = ga.chrom
                AND (stage.start = ga.start AND stage.end = ga.end)")

cols <- c(1,2,3,4,5,6,10,11)
hypo_diff_Pupae_vs_Adult <- output[,c(cols)]
write.table(diff_Pupae_vs_Adult_hypo, file="Pupae_vs_Adult_Adult_Hypomethylated_Exons_10%_sig.txt", quote=F, row.names = T, sep = '\t')
write.table(hypo_diff_Pupae_vs_Adult, file="Pupae_vs_Adult_Adult_Hypomethylated_Exons_10%_sig_w_gene.txt", quote=F, row.names = T, sep = '\t')
