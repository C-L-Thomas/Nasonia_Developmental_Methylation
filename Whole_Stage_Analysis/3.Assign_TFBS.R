library(grid)
library(readr)
library(ggplot2)
library(sqldf)

TFs <- read_delim("~/Nasonia_100%_TFBS.txt", "\t", escape_double = F, col_names = T, trim_ws = T) #Located in Genome Folder

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

adu <- read_delim("~/Adult_Significant_output.txt", "\t", escape_double = F, col_names = T, trim_ws = T)
colnames(adu) <- c("chr2", "start","cov","c_count", "pVal","FDR")
head(adu)

atfs <- sqldf("SELECT dms.chr2,
                dms.start,
                tf.chrom,
                tf.start,
                tf.end,
                tf.TF
                FROM adu AS dms
                LEFT JOIN new AS tf
                ON tf.chrom = dms.chr2
                AND (dms.start >= tf.start AND dms.start <= tf.end)")
head(atfs)
atfs <- na.omit(atfs[!is.na(atfs$TF),])
atfs[,7] <- paste(atfs[,1], atfs[,4], sep = "-")
atfs[,8] <- paste(atfs[,7], atfs[,6], sep = "-")
atfs <- atfs[!duplicated(atfs[,8]), ]
head(atfs)

write.table(atfs, file="Adult_Transcription_Factors_100%.txt",
            sep="\t", row.names=F, quote=F)

atfs$feature <- "TF"
atfs$feature_numb <- "NA"

cols <- c(1,2,3,9,4,5,6,10)
atfs <- atfs[,c(cols)]

colnames(atfs) <- c("chr","start","chr.1","feature","start.1","end","geneID","feature_numb")

write.table(atfs, file="Adult_Transcription_Factors_100%.txt",
            sep="\t", row.names=F, quote=F)
