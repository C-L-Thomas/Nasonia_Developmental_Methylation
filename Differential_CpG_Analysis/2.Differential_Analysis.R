library("methylKit")
library("sqldf")
library("genomation")
library(stringr)

load(file = "Counted_CpGs.RData")

#------------------------------------------------------------------------------
#Reorganise All Analyses

Embryo_vs_Larva=reorganize(subset_methBase,sample.ids=c("EMB5","EMB6","LAR4","LAR5","LAR6"),
                                treatment=c(0,0,1,1,1))

Larva_vs_Prepupae=reorganize(subset_methBase,sample.ids=c("LAR4","LAR5","LAR6","PRE1","PRE2","PRE3"),
                                treatment=c(1,1,1,2,2,2))

Prepupae_vs_Pupae=reorganize(subset_methBase,sample.ids=c("PRE1","PRE2","PRE3","PUP2","PUP3","PUP4"),
                                treatment=c(2,2,2,3,3,3))

Pupae_vs_Adult=reorganize(subset_methBase,sample.ids=c("PUP2","PUP3","PUP4","ADU1","ADU2","ADU4"),
                                treatment=c(3,3,3,4,4,4))

#------------------------------------------------------------------------------
#Embryo vs Larva

diff_Embryo_vs_Larva <- calculateDiffMeth(Embryo_vs_Larva, adjust = c("BH"))

write.table(diff_Embryo_vs_Larva, file="Embryo_vs_Larva_Complete_Result.txt", quote=F, row.names = T, sep = '\t')

diff_Embryo_vs_Larva_10 <- getMethylDiff(diff_Embryo_vs_Larva, difference=10, qvalue=0.0125)
write.table(diff_Embryo_vs_Larva_10, file="Embryo_vs_Larva_10%_sig.txt",quote=F, row.names = T, sep = '\t')
nrow(diff_Embryo_vs_Larva_10)

diff_Embryo_vs_Larva_hyper=getMethylDiff(diff_Embryo_vs_Larva,difference=10,qvalue=0.0125,type="hyper")
write.table(diff_Embryo_vs_Larva_hyper, file="Embryo_vs_Larva_Larva_10%_sig.txt", quote=F, row.names = T, sep = '\t')

diff_Embryo_vs_Larva_hypo=getMethylDiff(diff_Embryo_vs_Larva,difference=10,qvalue=0.0125,type="hypo")
write.table(diff_Embryo_vs_Larva_hypo, file="Embryo_vs_Larva_Embryo_10%_sig.txt", quote=F, row.names = T, sep = '\t')

# Annotating diff meth sites with gene IDs

genome_annotation<-read.csv.sql("/data/monoallelic/clt48/Nvit_psr_1.1/genes_with_start_and_end.txt",
                                sql ="select * from file", sep="\t",header = T)
colnames(genome_annotation) <- c("chr","start","end","geneID")

diff_E_vs_L <- getData(diff_Embryo_vs_Larva_10)

colnames(diff_E_vs_L)[7]<-"meth_diff"

output <- sqldf("SELECT diff.chr,
                diff.start,
                diff.end,
                diff.pvalue,
                diff.qvalue,
                diff.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM diff_E_vs_L AS diff
                LEFT JOIN genome_annotation AS ga
                ON diff.chr = ga.chr
                AND (diff.start >= ga.start AND diff.start <= ga.end)")

write.table(output, file="Embryo_vs_Larva_DMS_w10%diff_with_geneID.txt",
          sep="\t", row.names=F, quote=F)

diff_E_vs_L <- getData(diff_Embryo_vs_Larva_hyper)
colnames(diff_E_vs_L)[7]<-"meth_diff"

output <- sqldf("SELECT diff.chr,
                diff.start,
                diff.end,
                diff.pvalue,
                diff.qvalue,
                diff.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM diff_E_vs_L AS diff
                LEFT JOIN genome_annotation AS ga
                ON diff.chr = ga.chr
                AND (diff.start >= ga.start AND diff.start <= ga.end)")

write.table(output, file="Embryo_vs_Larva_DMS_Hypermethylated_with_geneID.txt",
          sep="\t", row.names=F, quote=F)

diff_E_vs_L <- getData(diff_Embryo_vs_Larva_hypo)
colnames(diff_E_vs_L)[7]<-"meth_diff"

output <- sqldf("SELECT diff.chr,
                diff.start,
                diff.end,
                diff.pvalue,
                diff.qvalue,
                diff.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM diff_E_vs_L AS diff
                LEFT JOIN genome_annotation AS ga
                ON diff.chr = ga.chr
                AND (diff.start >= ga.start AND diff.start <= ga.end)")

write.table(output, file="Embryo_vs_Larva_DMS_Hypomethylated_with_geneID.txt",
          sep="\t", row.names=F, quote=F)


#------------------------------------------------------------------------------
#Larva_vs_Prepupae

diff_Larva_vs_Prepupae <- calculateDiffMeth(Larva_vs_Prepupae, adjust = c("BH"))

write.table(diff_Larva_vs_Prepupae, file="Larva_vs_Prepupae_Complete_Result.txt", quote=F, row.names = T, sep = '\t')

diff_Larva_vs_Prepupae_10 <- getMethylDiff(diff_Larva_vs_Prepupae, difference=10, qvalue=0.0125)
write.table(diff_Larva_vs_Prepupae_10, file="Larva_vs_Prepupae_10%_sig.txt", quote=F, row.names = T, sep = '\t')
nrow(diff_Larva_vs_Prepupae_10)

diff_Larva_vs_Prepupae_hyper=getMethylDiff(diff_Larva_vs_Prepupae,difference=10,qvalue=0.0125,type="hyper")
write.table(diff_Larva_vs_Prepupae_hyper, file="Larva_vs_Prepupae_Prepupae_10%_sig.txt", quote=T, row.names = T, sep = '\t')

diff_Larva_vs_Prepupae_hypo=getMethylDiff(diff_Larva_vs_Prepupae,difference=10,qvalue=0.0125,type="hypo")
write.table(diff_Larva_vs_Prepupae_hypo, file="Larva_vs_Prepupae_Larva_10%_sig.txt", quote=F, row.names = T, sep = '\t')

diff_L_vs_P <- getData(diff_Larva_vs_Prepupae_10)

colnames(diff_L_vs_P)[7]<-"meth_diff"

output <- sqldf("SELECT diff.chr,
                diff.start,
                diff.end,
                diff.pvalue,
                diff.qvalue,
                diff.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM diff_L_vs_P AS diff
                LEFT JOIN genome_annotation AS ga
                ON diff.chr = ga.chr
                AND (diff.start >= ga.start AND diff.start <= ga.end)")

write.table(output, file="Larva_vs_Prepupae_DMS_w10%diff_with_geneID.txt",
          sep="\t", row.names=F, quote=F)

diff_L_vs_P <- getData(diff_Larva_vs_Prepupae_hyper)
colnames(diff_L_vs_P)[7]<-"meth_diff"

output <- sqldf("SELECT diff.chr,
                diff.start,
                diff.end,
                diff.pvalue,
                diff.qvalue,
                diff.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM diff_L_vs_P AS diff
                LEFT JOIN genome_annotation AS ga
                ON diff.chr = ga.chr
                AND (diff.start >= ga.start AND diff.start <= ga.end)")

write.table(output, file="Larva_vs_Prepupae_DMS_Hypermethylated_with_geneID.txt",
          sep="\t", row.names=F, quote=F)

diff_L_vs_P <- getData(diff_Larva_vs_Prepupae_hypo)
colnames(diff_L_vs_P)[7]<-"meth_diff"

output <- sqldf("SELECT diff.chr,
                diff.start,
                diff.end,
                diff.pvalue,
                diff.qvalue,
                diff.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM diff_L_vs_P AS diff
                LEFT JOIN genome_annotation AS ga
                ON diff.chr = ga.chr
                AND (diff.start >= ga.start AND diff.start <= ga.end)")

write.table(output, file="Larva_vs_Prepupae_DMS_Hypomethylated_with_geneID.txt",
          sep="\t", row.names=F, quote=F)

#------------------------------------------------------------------------------
#Prepupae_vs_Pupae                

diff_Prepupae_vs_Pupae <- calculateDiffMeth(Prepupae_vs_Pupae, adjust = c("BH"))

write.table(diff_Prepupae_vs_Pupae, file="Prepupae_vs_Pupae_Complete_Result.txt", quote=F, row.names = T, sep = '\t')

diff_Prepupae_vs_Pupae_10 <- getMethylDiff(diff_Prepupae_vs_Pupae, difference=10, qvalue=0.0125)
write.table(diff_Prepupae_vs_Pupae_10, file="Prepupae_vs_Pupae_10%_sig.txt", quote=F, row.names = T, sep = '\t')

diff_Prepupae_vs_Pupae_hyper=getMethylDiff(diff_Prepupae_vs_Pupae,difference=10,qvalue=0.0125,type="hyper")
write.table(diff_Prepupae_vs_Pupae_hyper, file="Prepupae_vs_Pupae_Pupae_10%_sig.txt", quote=F, row.names = T, sep = '\t')

diff_Prepupae_vs_Pupae_hypo=getMethylDiff(diff_Prepupae_vs_Pupae,difference=10,qvalue=0.0125,type="hypo")
write.table(diff_Prepupae_vs_Pupae_hypo, file="Prepupae_vs_Pupae_Prepupae_10%_sig.txt", quote=F, row.names = T, sep = '\t')

diff_P_vs_P <- getData(diff_Prepupae_vs_Pupae_10)

colnames(diff_P_vs_P)[7]<-"meth_diff"
head(diff_P_vs_P)
head(genome_annotation)

output <- sqldf("SELECT diff.chr,
                diff.start,
                diff.end,
                diff.pvalue,
                diff.qvalue,
                diff.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM diff_P_vs_P AS diff
                LEFT JOIN genome_annotation AS ga
                ON diff.chr = ga.chr
                AND (diff.start >= ga.start AND diff.start <= ga.end)")

write.table(output, file="Prepupae_vs_Pupae_DMS_w10%diff_with_geneID.txt",
          sep="\t", row.names=F, quote=F)

diff_P_vs_P <- getData(diff_Prepupae_vs_Pupae_hyper)
colnames(diff_P_vs_P)[7]<-"meth_diff"

output <- sqldf("SELECT diff.chr,
                diff.start,
                diff.end,
                diff.pvalue,
                diff.qvalue,
                diff.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM diff_P_vs_P AS diff
                LEFT JOIN genome_annotation AS ga
                ON diff.chr = ga.chr
                AND (diff.start >= ga.start AND diff.start <= ga.end)")

write.table(output, file="Prepupae_vs_Pupae_DMS_Hypermethylated_with_geneID.txt",
          sep="\t", row.names=F, quote=F)

diff_P_vs_P <- getData(diff_Prepupae_vs_Pupae_hypo)
colnames(diff_P_vs_P)[7]<-"meth_diff"

output <- sqldf("SELECT diff.chr,
                diff.start,
                diff.end,
                diff.pvalue,
                diff.qvalue,
                diff.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM diff_P_vs_P AS diff
                LEFT JOIN genome_annotation AS ga
                ON diff.chr = ga.chr
                AND (diff.start >= ga.start AND diff.start <= ga.end)")

write.table(output, file="Prepupae_vs_Pupae_DMS_Hypomethylated_with_geneID.txt",
          sep="\t", row.names=F, quote=F)

#------------------------------------------------------------------------------
#Pupae_vs_Adult

diff_Pupae_vs_Adult <- calculateDiffMeth(Pupae_vs_Adult, adjust = c("BH"))

write.table(diff_Pupae_vs_Adult, file="Pupae_vs_Adult_Complete_Result.txt", quote=F, row.names = F, sep = '\t')

diff_Pupae_vs_Adult_10 <- getMethylDiff(diff_Pupae_vs_Adult, difference=10, qvalue=0.0125)
write.table(diff_Pupae_vs_Adult_10, file="Pupae_vs_Adult_10%_sig.txt", quote=F, row.names = T, sep = '\t')

diff_Pupae_vs_Adult_hyper=getMethylDiff(diff_Pupae_vs_Adult,difference=10,qvalue=0.0125,type="hyper")
write.table(diff_Pupae_vs_Adult_hyper, file="Pupae_vs_Adult_Adult_10%_sig.txt", quote=F, row.names = T, sep = '\t')

diff_Pupae_vs_Adult_hypo=getMethylDiff(diff_Pupae_vs_Adult,difference=10,qvalue=0.0125,type="hypo")
write.table(diff_Pupae_vs_Adult_hypo, file="Pupae_vs_Adult_Pupae_10%_sig.txt", quote=F, row.names = T, sep = '\t')

diff_P_vs_A <- getData(diff_Pupae_vs_Adult_10)

colnames(diff_P_vs_A)[7]<-"meth_diff"
head(diff_P_vs_A)
head(genome_annotation)

output <- sqldf("SELECT diff.chr,
                diff.start,
                diff.end,
                diff.pvalue,
                diff.qvalue,
                diff.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM diff_P_vs_A AS diff
                LEFT JOIN genome_annotation AS ga
                ON diff.chr = ga.chr
                AND (diff.start >= ga.start AND diff.start <= ga.end)")

write.table(output, file="Pupae_vs_Adult_DMS_w10%diff_with_geneID.txt",
          sep="\t", row.names=F, quote=F)

diff_P_vs_A <- getData(diff_Pupae_vs_Adult_hyper)
colnames(diff_P_vs_A)[7]<-"meth_diff"

output <- sqldf("SELECT diff.chr,
                diff.start,
                diff.end,
                diff.pvalue,
                diff.qvalue,
                diff.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM diff_P_vs_A AS diff
                LEFT JOIN genome_annotation AS ga
                ON diff.chr = ga.chr
                AND (diff.start >= ga.start AND diff.start <= ga.end)")

write.table(output, file="Pupae_vs_Adult_DMS_Hypermethylated_with_geneID.txt",
          sep="\t", row.names=F, quote=F)

diff_P_vs_A <- getData(diff_Pupae_vs_Adult_hypo)
colnames(diff_P_vs_A)[7]<-"meth_diff"

output <- sqldf("SELECT diff.chr,
                diff.start,
                diff.end,
                diff.pvalue,
                diff.qvalue,
                diff.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM diff_P_vs_A AS diff
                LEFT JOIN genome_annotation AS ga
                ON diff.chr = ga.chr
                AND (diff.start >= ga.start AND diff.start <= ga.end)")

write.table(output, file="Pupae_vs_Adult_DMS_Hypomethylated_with_geneID.txt",
          sep="\t", row.names=F, quote=F)

#------------------------------------------------------------------------------
#Closest TSS

gene_info <-read.table("~/GCF_009193385.2_Nvit_psr_1.1_feature_table.txt",
                                 sep="\t",header = T)
colnames(gene_info) <- c("feature", "class", "assembly","assembly_unit","seq_type","chromosome",
                        "genomic_accession","start","end","strand","product_accession","non-redundant_refseq",
                        "related_accession","name","symbol","GeneID","locus_tag","feature_interval_length",
                        "product_length","attributes")

gene.obj=readTranscriptFeatures("~/GCF_009193385.2_Nvit_psr_1.1_genomic.bed.txt",
up.flank = 100, down.flank = 100, unique.prom = FALSE)

diff_Embryo_vs_Larva_10$chr <- str_replace(diff_Embryo_vs_Larva_10$chr, "_",".")
ann <- genomation::annotateWithGeneParts(as(diff_Embryo_vs_Larva_10,"GRanges"),gene.obj)
distance <- getAssociationWithTSS(ann)

colnames(distance) <- c("row","disttss","transcript","strand")
distance$transcript<-gsub("rna-","",as.character(distance$transcript))
nrow(distance)
nrow(diff_Embryo_vs_Larva_10)

#If distance and diff_Embryo don't match this is likely due to DMS
#on chromosome without gene. Here is how to check which chromosome / row number
diff_Embryo_vs_Larva_10$row <- seq.int(nrow(diff_Embryo_vs_Larva_10))
diff_Embryo_vs_Larva_10[!(diff_Embryo_vs_Larva_10$row %in% distance$row)]

clean_Embryo_vs_Larva_10 <- diff_Embryo_vs_Larva_10[-c(116615), ]
nrow(clean_Embryo_vs_Larva_10)

output <- sqldf("SELECT distance.row,
                distance.disttss,
                distance.transcript,
                distance.strand,
                ga.genomic_accession,
                ga.start,
                ga.end,
                ga.strand,
                ga.product_accession,
                ga.symbol
                FROM distance AS distance
                LEFT JOIN gene_info AS ga
                ON distance.transcript = ga.product_accession")

output$location <- clean_Embryo_vs_Larva_10$start

write.table(output, file="Embryo_vs_Larva_10%_sig_with_distance.txt",quote=F)

diff_Larva_vs_Prepupae_10$chr <- str_replace(diff_Larva_vs_Prepupae_10$chr, "_",".")
ann <- genomation::annotateWithGeneParts(as(diff_Larva_vs_Prepupae_10,"GRanges"),gene.obj)
distance <- getAssociationWithTSS(ann)
colnames(distance) <- c("row","disttss","transcript","strand")
distance$transcript<-gsub("rna-","",as.character(distance$transcript))

output <- sqldf("SELECT distance.row,
                distance.disttss,
                distance.transcript,
                distance.strand,
                ga.genomic_accession,
                ga.start,
                ga.end,
                ga.strand,
                ga.product_accession,
                ga.symbol
                FROM distance AS distance
                LEFT JOIN gene_info AS ga
                ON distance.transcript = ga.product_accession")

output$location <- diff_Larva_vs_Prepupae_10$start

write.table(output, file="Larva_vs_Prepupae_10%_sig_with_distance.txt",quote=F)

diff_Prepupae_vs_Pupae_10$chr <- str_replace(diff_Prepupae_vs_Pupae_10$chr, "_",".")
ann <- genomation::annotateWithGeneParts(as(diff_Prepupae_vs_Pupae_10,"GRanges"),gene.obj)
distance <- getAssociationWithTSS(ann)
colnames(distance) <- c("row","disttss","transcript","strand")
distance$transcript<-gsub("rna-","",as.character(distance$transcript))

output <- sqldf("SELECT distance.row,
                distance.disttss,
                distance.transcript,
                distance.strand,
                ga.genomic_accession,
                ga.start,
                ga.end,
                ga.strand,
                ga.product_accession,
                ga.symbol
                FROM distance AS distance
                LEFT JOIN gene_info AS ga
                ON distance.transcript = ga.product_accession")

output$location <- diff_Prepupae_vs_Pupae_10$start

write.table(output, file="Prepupae_vs_Pupae_10%_sig_with_distance.txt",quote=F)

diff_Pupae_vs_Adult_10$chr <- str_replace(diff_Pupae_vs_Adult_10$chr, "_",".")
ann <- genomation::annotateWithGeneParts(as(diff_Pupae_vs_Adult_10,"GRanges"),gene.obj)
distance <- getAssociationWithTSS(ann)
colnames(distance) <- c("row","disttss","transcript","strand")
distance$transcript<-gsub("rna-","",as.character(distance$transcript))

output <- sqldf("SELECT distance.row,
                distance.disttss,
                distance.transcript,
                distance.strand,
                ga.genomic_accession,
                ga.start,
                ga.end,
                ga.strand,
                ga.product_accession,
                ga.symbol
                FROM distance AS distance
                LEFT JOIN gene_info AS ga
                ON distance.transcript = ga.product_accession")

output$location <- diff_Pupae_vs_Adult_10$start

write.table(output, file="Pupae_vs_Adult_10%_sig_with_distance.txt",quote=F)
