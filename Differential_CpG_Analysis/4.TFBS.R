library(methylKit)
library(grid)
library(readr)
library(ggplot2)
library(sqldf)

TFs <- read_delim("~/Nasonia_100%_TFBS.txt", "\t", escape_double = F, col_names = T, trim_ws = T)

cols <- c(1,4,5,10)
TFs <- TFs[,c(cols)]
colnames(TFs) <- c("chr", "start", "end" ,"TF")
head(TFs)

TFs$chrom <- gsub(" .*", "", TFs$chr)
head(TFs$chrom)
head(TFs)

cols <- c(2,3,4,5)
new <- TFs[,c(cols)]
head(new)

adu <- read.table("Embryo_vs_Larva_10%_sig_with_distance.txt", sep = " ", header = T)
colnames(adu) <- c("row", "disttss", "transcript", "strand", "genomic_accession","Start", "End", "Strand2", "product_accession", "symbol", "location")
head(adu)

atfs <- sqldf("SELECT dms.disttss,
                dms.strand,
                dms.Start,
                dms.End,
                dms.genomic_accession,
                dms.location,
                dms.symbol,
                tf.chrom,
                tf.start,
                tf.end,
                tf.TF
                FROM adu AS dms
                LEFT JOIN new AS tf
                ON tf.chrom = dms.genomic_accession
                AND (dms.location >= tf.start AND dms.location <= tf.end)")

head(atfs)
atfs <- na.omit(atfs[!is.na(atfs$TF),])
head(atfs)
cols <- c(1,2,6,5,7,11)
atfs <- atfs[,c(cols)]


write.table(atfs, file="Embryo_vs_Larva_10%_sig_TFBS_distance.txt",
            sep="\t", row.names=F, quote=F)
