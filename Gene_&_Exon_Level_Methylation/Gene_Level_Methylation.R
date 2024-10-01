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
head(meth_all_data)

save.image(file = "United.RData")

#load(file = "United.RData")

#Input Gene Files
gene.obj = gffToGRanges("~/GCF_009193385.2_Nvit_psr_1.1_genomic.gff")
gene_bed <- subset(gene.obj, type =="gene" )
gene_coverage=regionCounts(meth_all_data,gene_bed)

genes<-read.csv.sql("~/genes_with_start_and_end.txt",
                                sql ="select * from file", sep="\t",header = T)
colnames(genes) <- c("chrom","Start","End","geneID")

gene_coverage <- getData(gene_coverage)
head(gene_coverage)

write.table(gene_coverage, "Annotation_Table.txt", quote=F, row.names = F, sep = '\t' )

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
                ga.geneID
                FROM gene_coverage AS out
                LEFT JOIN genes AS ga
                ON out.chr = ga.chrom
                AND out.start = ga.Start
                AND out.end = ga.End")

write.table(output, "Annotation_Table.txt", quote=F, row.names = F, sep = '\t' )

output$embryo1 <- (output$numCs1 / output$coverage1)*100
output$embryo2 <- (output$numCs2 / output$coverage2)*100
output$larva1 <- (output$numCs3 / output$coverage3)*100
output$larva2 <- (output$numCs4 / output$coverage4)*100
output$larva3 <- (output$numCs5 / output$coverage5)*100
output$prepupa1 <- (output$numCs6 / output$coverage6)*100
output$prepupa2 <- (output$numCs7 / output$coverage7)*100
output$prepupa3 <- (output$numCs8 / output$coverage8)*100
output$pupa1 <- (output$numCs9 / output$coverage9)*100
output$pupa2 <- (output$numCs10 / output$coverage10)*100
output$pupa3 <- (output$numCs11 / output$coverage11)*100
output$adult1 <- (output$numCs12 / output$coverage12)*100
output$adult2 <- (output$numCs13 / output$coverage13)*100
output$adult3 <- (output$numCs14 / output$coverage14)*100

output <- output[,c(1,2,3,4,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64)]

write.table(output, "Whole_Development_Methylation_Percentages.txt", quote=F, row.names = F, sep = '\t' )

#------------------------------------------------------------------------------

# Generate some quick plots which can be fixed up later
pdf("Correlation.pdf")
getCorrelation(gene_coverage,plot=TRUE)
dev.off()

pdf("Cluster.pdf")
clusterSamples(gene_coverage, dist="correlation", method="ward", plot=TRUE)
dev.off()

pdf("PCA.pdf")
PCASamples(gene_coverage, adj.lim = c(1.5,0.15))
dev.off()

#______________________________________________________________________________
#Reorganise All Analyses

Embryo_vs_Larva=reorganize(gene_coverage,sample.ids=c("EMB1","EMB2","LAR1","LAR2","LAR3"),
                                treatment=c(0,0,1,1,1))

Larva_vs_Prepupae=reorganize(gene_coverage,sample.ids=c("LAR1","LAR2","LAR3","PRE1","PRE2","PRE3"),
                                treatment=c(1,1,1,2,2,2))

Prepupae_vs_Pupae=reorganize(gene_coverage,sample.ids=c("PRE1","PRE2","PRE3","PUP1","PUP2","PUP3"),
                                treatment=c(2,2,2,3,3,3))

Pupae_vs_Adult=reorganize(gene_coverage,sample.ids=c("PUP1","PUP2","PUP3","ADU1","ADU2","ADU3"),
                                treatment=c(3,3,3,4,4,4))

#------------------------------------------------------------------------------

diff_Embryo_vs_Larva <- calculateDiffMeth(Embryo_vs_Larva, adjust = c("BH"))
nrow(diff_Embryo_vs_Larva)
head(diff_Embryo_vs_Larva)

EvL <- getData(diff_Embryo_vs_Larva)
colnames(EvL) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")
head(EvL)

genome_annotation <- read.delim("/data/monoallelic/clt48/Nvit_psr_1.1/genes_with_start_and_end.txt", header = FALSE)
colnames(genome_annotation) <- c("chr","start","end","geneID")
head(genome_annotation)

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM EvL AS stage
                LEFT JOIN genome_annotation AS ga
                ON stage.chr = ga.chr
                AND (stage.start = ga.start AND stage.end = ga.end)")

head(output)

cols <- c(1,2,3,4,5,6,10)
total_diff_Embryo_vs_Larva <- output[,c(cols)]

write.table(total_diff_Embryo_vs_Larva, file="Embryo_vs_Larva_Complete_Gene_Result.txt", quote=F, row.names = F, sep = '\t')

diff_Embryo_vs_Larva_10 <- getMethylDiff(diff_Embryo_vs_Larva, difference=10, qvalue=0.0125)

EvL <- getData(diff_Embryo_vs_Larva_10)
nrow(EvL)
colnames(EvL) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM EvL AS stage
                LEFT JOIN genome_annotation AS ga
                ON stage.chr = ga.chr
                AND (stage.start = ga.start AND stage.end = ga.end)")

cols <- c(1,2,3,4,5,6,10)
total_diff_Embryo_vs_Larva_10 <- output[,c(cols)]
write.table(total_diff_Embryo_vs_Larva_10, file="Embryo_vs_Larva_Gene_10%_sig.txt", quote=F, row.names = F, sep = '\t')
nrow(diff_Embryo_vs_Larva_10)

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
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM EvL AS stage
                LEFT JOIN genome_annotation AS ga
                ON stage.chr = ga.chr
                AND (stage.start = ga.start AND stage.end = ga.end)")

cols <- c(1,2,3,4,5,6,10)
hyper_diff_Embryo_vs_Larva <- output[,c(cols)]

write.table(hyper_diff_Embryo_vs_Larva, file="Embryo_vs_Larva_Larva_Hypermethylated_Genes_10%_sig.txt", quote=F, row.names = T, sep = '\t')

diff_Embryo_vs_Larva_hypo=getMethylDiff(diff_Embryo_vs_Larva,difference=10,qvalue=0.0125,type="hypo")

EvL <- getData(diff_Embryo_vs_Larva_hypo)
nrow(EvL)
colnames(EvL) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM EvL AS stage
                LEFT JOIN genome_annotation AS ga
                ON stage.chr = ga.chr
                AND (stage.start = ga.start AND stage.end = ga.end)")

cols <- c(1,2,3,4,5,6,10)
hypo_diff_Embryo_vs_Larva <- output[,c(cols)]

write.table(hypo_diff_Embryo_vs_Larva, file="Embryo_vs_Larva_Embryo_Hypermethylated_Genes_10%_sig.txt", quote=F, row.names = T, sep = '\t')

#------------------------------------------------------------------------------

diff_Larva_vs_Prepupae <- calculateDiffMeth(Larva_vs_Prepupae, adjust = c("BH"))
nrow(diff_Larva_vs_Prepupae)
head(diff_Larva_vs_Prepupae)

EvL <- getData(diff_Larva_vs_Prepupae)
colnames(EvL) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")
head(EvL)

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM EvL AS stage
                LEFT JOIN genome_annotation AS ga
                ON stage.chr = ga.chr
                AND (stage.start = ga.start AND stage.end = ga.end)")

head(output)

cols <- c(1,2,3,4,5,6,10)
total_diff_Larva_vs_Prepupae <- output[,c(cols)]

write.table(total_diff_Larva_vs_Prepupae, file="Larva_vs_Prepupae_Complete_Gene_Result.txt", quote=F, row.names = F, sep = '\t')

diff_Larva_vs_Prepupae_10 <- getMethylDiff(diff_Larva_vs_Prepupae, difference=10, qvalue=0.0125)

EvL <- getData(diff_Larva_vs_Prepupae_10)
nrow(EvL)
colnames(EvL) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM EvL AS stage
                LEFT JOIN genome_annotation AS ga
                ON stage.chr = ga.chr
                AND (stage.start = ga.start AND stage.end = ga.end)")

cols <- c(1,2,3,4,5,6,10)
total_diff_Larva_vs_Prepupae_10 <- output[,c(cols)]
write.table(total_diff_Larva_vs_Prepupae_10, file="Larva_vs_Prepupae_Gene_10%_sig.txt", quote=F, row.names = F, sep = '\t')
nrow(diff_Larva_vs_Prepupae_10)

diff_Larva_vs_Prepupae_hyper=getMethylDiff(diff_Larva_vs_Prepupae,difference=10,qvalue=0.0125,type="hyper")
EvL <- getData(diff_Larva_vs_Prepupae_hyper)

colnames(EvL) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM EvL AS stage
                LEFT JOIN genome_annotation AS ga
                ON stage.chr = ga.chr
                AND (stage.start = ga.start AND stage.end = ga.end)")

cols <- c(1,2,3,4,5,6,10)
hyper_diff_Larva_vs_Prepupae <- output[,c(cols)]

write.table(hyper_diff_Larva_vs_Prepupae, file="Larva_vs_Prepupae_Prepupae_Hypermethylated_Genes_10%_sig.txt", quote=F, row.names = T, sep = '\t')

diff_Larva_vs_Prepupae_hypo=getMethylDiff(diff_Larva_vs_Prepupae,difference=10,qvalue=0.0125,type="hypo")

EvL <- getData(diff_Larva_vs_Prepupae_hypo)
nrow(EvL)
colnames(EvL) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM EvL AS stage
                LEFT JOIN genome_annotation AS ga
                ON stage.chr = ga.chr
                AND (stage.start = ga.start AND stage.end = ga.end)")

cols <- c(1,2,3,4,5,6,10)
hypo_diff_Larva_vs_Prepupae <- output[,c(cols)]

write.table(hypo_diff_Larva_vs_Prepupae, file="Larva_vs_Prepupae_Larva_Hypermethylated_Genes_10%_sig.txt", quote=F, row.names = T, sep = '\t')

#------------------------------------------------------------------------------

diff_Prepupae_vs_Pupae <- calculateDiffMeth(Prepupae_vs_Pupae, adjust = c("BH"))
nrow(diff_Prepupae_vs_Pupae)
head(diff_Prepupae_vs_Pupae)

EvL <- getData(diff_Prepupae_vs_Pupae)
colnames(EvL) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")
head(EvL)

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM EvL AS stage
                LEFT JOIN genome_annotation AS ga
                ON stage.chr = ga.chr
                AND (stage.start = ga.start AND stage.end = ga.end)")

head(output)

cols <- c(1,2,3,4,5,6,10)
total_diff_Prepupae_vs_Pupae <- output[,c(cols)]

write.table(total_diff_Prepupae_vs_Pupae, file="Prepupae_vs_Pupae_Complete_Gene_Result.txt", quote=F, row.names = F, sep = '\t')

diff_Prepupae_vs_Pupae_10 <- getMethylDiff(diff_Prepupae_vs_Pupae, difference=10, qvalue=0.0125)

EvL <- getData(diff_Prepupae_vs_Pupae_10)
nrow(EvL)
colnames(EvL) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM EvL AS stage
                LEFT JOIN genome_annotation AS ga
                ON stage.chr = ga.chr
                AND (stage.start = ga.start AND stage.end = ga.end)")

cols <- c(1,2,3,4,5,6,10)
total_diff_Prepupae_vs_Pupae_10 <- output[,c(cols)]
write.table(total_diff_Prepupae_vs_Pupae_10, file="Prepupae_vs_Pupae_Gene_10%_sig.txt", quote=F, row.names = F, sep = '\t')
nrow(diff_Prepupae_vs_Pupae_10)

diff_Prepupae_vs_Pupae_hyper=getMethylDiff(diff_Prepupae_vs_Pupae,difference=10,qvalue=0.0125,type="hyper")
EvL <- getData(diff_Prepupae_vs_Pupae_hyper)

colnames(EvL) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM EvL AS stage
                LEFT JOIN genome_annotation AS ga
                ON stage.chr = ga.chr
                AND (stage.start = ga.start AND stage.end = ga.end)")

cols <- c(1,2,3,4,5,6,10)
hyper_diff_Prepupae_vs_Pupae <- output[,c(cols)]

write.table(hyper_diff_Prepupae_vs_Pupae, file="Prepupae_vs_Pupae_Pupae_Hypermethylated_Genes_10%_sig.txt", quote=F, row.names = T, sep = '\t')

diff_Prepupae_vs_Pupae_hypo=getMethylDiff(diff_Prepupae_vs_Pupae,difference=10,qvalue=0.0125,type="hypo")

EvL <- getData(diff_Prepupae_vs_Pupae_hypo)
nrow(EvL)
colnames(EvL) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM EvL AS stage
                LEFT JOIN genome_annotation AS ga
                ON stage.chr = ga.chr
                AND (stage.start = ga.start AND stage.end = ga.end)")

cols <- c(1,2,3,4,5,6,10)
hypo_diff_Prepupae_vs_Pupae <- output[,c(cols)]

write.table(hypo_diff_Prepupae_vs_Pupae, file="Prepupae_vs_Pupae_Prepupae_Hypermethylated_Genes_10%_sig.txt", quote=F, row.names = T, sep = '\t')

#------------------------------------------------------------------------------

diff_Pupae_vs_Adult <- calculateDiffMeth(Pupae_vs_Adult, adjust = c("BH"))
nrow(diff_Pupae_vs_Adult)
head(diff_Pupae_vs_Adult)

EvL <- getData(diff_Pupae_vs_Adult)
colnames(EvL) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")
head(EvL)

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM EvL AS stage
                LEFT JOIN genome_annotation AS ga
                ON stage.chr = ga.chr
                AND (stage.start = ga.start AND stage.end = ga.end)")

head(output)

cols <- c(1,2,3,4,5,6,10)
total_diff_Pupae_vs_Adult <- output[,c(cols)]

write.table(total_diff_Pupae_vs_Adult, file="Pupae_vs_Adult_Complete_Gene_Result.txt", quote=F, row.names = F, sep = '\t')

diff_Pupae_vs_Adult_10 <- getMethylDiff(diff_Pupae_vs_Adult, difference=10, qvalue=0.0125)

EvL <- getData(diff_Pupae_vs_Adult_10)
nrow(EvL)
colnames(EvL) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM EvL AS stage
                LEFT JOIN genome_annotation AS ga
                ON stage.chr = ga.chr
                AND (stage.start = ga.start AND stage.end = ga.end)")

cols <- c(1,2,3,4,5,6,10)
total_diff_Pupae_vs_Adult_10 <- output[,c(cols)]
write.table(total_diff_Pupae_vs_Adult_10, file="Pupae_vs_Adult_Gene_10%_sig.txt", quote=F, row.names = F, sep = '\t')
nrow(diff_Pupae_vs_Adult_10)

diff_Pupae_vs_Adult_hyper=getMethylDiff(diff_Pupae_vs_Adult,difference=10,qvalue=0.0125,type="hyper")
EvL <- getData(diff_Pupae_vs_Adult_hyper)

colnames(EvL) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM EvL AS stage
                LEFT JOIN genome_annotation AS ga
                ON stage.chr = ga.chr
                AND (stage.start = ga.start AND stage.end = ga.end)")

cols <- c(1,2,3,4,5,6,10)
hyper_diff_Pupae_vs_Adult <- output[,c(cols)]

write.table(hyper_diff_Pupae_vs_Adult, file="Pupae_vs_Adult_Adult_Hypermethylated_Genes_10%_sig.txt", quote=F, row.names = T, sep = '\t')

diff_Pupae_vs_Adult_hypo=getMethylDiff(diff_Pupae_vs_Adult,difference=10,qvalue=0.0125,type="hypo")

EvL <- getData(diff_Pupae_vs_Adult_hypo)
nrow(EvL)
colnames(EvL) <- c("chr","start","end","strand","pvalue","qvalue","meth_diff")

output <- sqldf("SELECT stage.chr,
                stage.start,
                stage.end,
                stage.pvalue,
                stage.qvalue,
                stage.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM EvL AS stage
                LEFT JOIN genome_annotation AS ga
                ON stage.chr = ga.chr
                AND (stage.start = ga.start AND stage.end = ga.end)")

cols <- c(1,2,3,4,5,6,10)
hypo_diff_Pupae_vs_Adult <- output[,c(cols)]

write.table(hypo_diff_Pupae_vs_Adult, file="Pupae_vs_Adult_Pupae_Hypermethylated_Genes_10%_sig.txt", quote=F, row.names = T, sep = '\t')
