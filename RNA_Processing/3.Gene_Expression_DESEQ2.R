library('DESeq2')
library('ggplot2')
library('topGO')

basedir = '~/Gene_Expression/'
setwd(basedir)
Counts = "Counts"
sampleDataFilename = "Gene_Metadata.txt"

Embryo = 'Emb'
Larva = 'Lar'
Prepupae = 'Pre'
Pupae = 'Pup'
Adult = 'Adu'

sampleTable = read.table(sampleDataFilename, header=TRUE)
head(sampleTable)

dds <-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory = Counts, design = ~ condition)

keep <- rowSums(counts(dds)) >= 1 
dds <- dds[keep,]
ddsTC <- DESeq(dds, test="LRT", reduced = ~ 1)
rlg<- rlog(ddsTC)
plotPCA(rlg)

#############################################################################################
#To compare gene counts
#Embryo

plot(log2(counts(ddsTC, normalized=TRUE)[,c(1,2)] + 1),
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)

plot(log2(counts(ddsTC, normalized=TRUE)[,c(2,3)] + 1),
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)

plot(log2(counts(ddsTC, normalized=TRUE)[,c(1,3)] + 1),
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)

#Larva
plot(log2(counts(ddsTC, normalized=TRUE)[,c(4,5)] + 1),
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)

plot(log2(counts(ddsTC, normalized=TRUE)[,c(5,6)] + 1),
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)

plot(log2(counts(ddsTC, normalized=TRUE)[,c(4,6)] + 1),
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)

#Prepupae
plot(log2(counts(ddsTC, normalized=TRUE)[,c(7,8)] + 1),
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)

plot(log2(counts(ddsTC, normalized=TRUE)[,c(8,9)] + 1),
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)

plot(log2(counts(ddsTC, normalized=TRUE)[,c(7,9)] + 1),
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)

#Pupae
plot(log2(counts(ddsTC, normalized=TRUE)[,c(10,11)] + 1),
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)

plot(log2(counts(ddsTC, normalized=TRUE)[,c(11,12)] + 1),
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)

plot(log2(counts(ddsTC, normalized=TRUE)[,c(10,12)] + 1),
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)

#Adult
plot(log2(counts(ddsTC, normalized=TRUE)[,c(13,14)] + 1),
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)

plot(log2(counts(ddsTC, normalized=TRUE)[,c(14,15)] + 1),
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)

plot(log2(counts(ddsTC, normalized=TRUE)[,c(13,15)] + 1),
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)

#############################################################################################

resemblar <- results(ddsTC, contrast = c('condition', 'Emb', 'Lar'), alpha = 0.0125)
summary(resemblar, alpha = 0.0125)
resemblar$Gene <- row.names(resemblar)
write.table(resemblar,"Results/EmbryovsLarva_Gene_results.txt",row.names = F, col.names= F,sep = "\t",quote = F)
lfc_emblar <- subset(resemblar ,log2FoldChange >  1.5 | log2FoldChange < -1.5 )
significant_emblar <- subset(lfc_emblar , padj < 0.0125 )
nrow(significant_emblar)
write.table(significant_emblar,"Results/EmbryovsLarva_Significant_Genes.txt",row.names = F, col.names= F,sep = "\t",quote = F)
upregulated_emblar <- subset(resemblar , padj < 0.0125 & log2FoldChange >  1.5)
nrow(upregulated_emblar)
write.table(upregulated_emblar,"Results/EmbryovsLarva_Embryo_Upregulated_Genes.txt",row.names = F, col.names= F,sep = "\t",quote = F)
downregulated_emblar<- subset(resemblar , padj < 0.0125 & log2FoldChange <  -1.5)
nrow(downregulated_emblar)
write.table(downregulated_emblar,"Results/EmbryovsLarva_Larva_Upregulated_Genes.txt",row.names = F, col.names= F,sep = "\t",quote = F)

reslarpre <- results(ddsTC, contrast = c('condition', 'Lar', 'Pre'), alpha = 0.0125)
summary(reslarpre, alpha = 0.0125)
reslarpre$Gene <- row.names(reslarpre)
write.table(reslarpre,"Results/LarvavsPrepupae_Gene_results.txt",row.names = F, col.names= F,sep = "\t",quote = F)
lfc_larpre <- subset(reslarpre ,log2FoldChange >  1.5 | log2FoldChange < -1.5 )
significant_larpre <- subset(lfc_larpre , padj < 0.0125 )
nrow(significant_larpre)
write.table(significant_larpre,"Results/LarvavsPrepupae_Significant_Genes.txt",row.names = F, col.names= F,sep = "\t",quote = F)
upregulated_larpre <- subset(reslarpre , padj < 0.0125 & log2FoldChange >  1.5)
nrow(upregulated_larpre)
write.table(upregulated_larpre,"Results/LarvavsPrepupae_Larva_Upregulated_Genes.txt",row.names = F, col.names= F,sep = "\t",quote = F)
downregulated_larpre<- subset(reslarpre , padj < 0.0125 & log2FoldChange <  -1.5)
nrow(downregulated_larpre)
write.table(downregulated_larpre,"Results/LarvavsPrepupae_Prepupae_Upregulated_Genes.txt",row.names = F, col.names= F,sep = "\t",quote = F)

resprepup <- results(ddsTC, contrast = c('condition', 'Pre', 'Pup'), alpha = 0.0125)
summary(resprepup, alpha = 0.0125)
resprepup$Gene <- row.names(resprepup)
write.table(resprepup,"Results/PrepupaevsPupae_Gene_results.txt",row.names = F, col.names= F,sep = "\t",quote = F)
lfc_prepup <- subset(resprepup ,log2FoldChange >  1.5 | log2FoldChange < -1.5 )
significant_prepup <- subset(lfc_prepup , padj < 0.0125 )
nrow(significant_prepup)
write.table(significant_prepup,"Results/PrepupaevsPupae_Significant_Genes.txt",row.names = F, col.names= F,sep = "\t",quote = F)
upregulated_prepup <- subset(resprepup , padj < 0.0125 & log2FoldChange >  1.5)
nrow(upregulated_prepup)
write.table(upregulated_prepup,"Results/PrepupaevsPupae_Prepupae_Upregulated_Genes.txt",row.names = F, col.names= F,sep = "\t",quote = F)
downregulated_prepup<- subset(resprepup , padj < 0.0125 & log2FoldChange <  -1.5)
nrow(downregulated_prepup)
write.table(downregulated_prepup,"Results/PrepupaevsPupae_Pupae_Upregulated_Genes.txt",row.names = F, col.names= F,sep = "\t",quote = F)

respupadu <- results(ddsTC, contrast = c('condition', 'Pup', 'Adu'), alpha = 0.0125)
summary(respupadu, alpha = 0.0125)
respupadu$Gene <- row.names(respupadu)
write.table(respupadu,"Results/PupaevsAdult_Gene_results.txt",row.names = F, col.names= F,sep = "\t",quote = F)
lfc_pupadu <- subset(respupadu ,log2FoldChange >  1.5 | log2FoldChange < -1.5 )
significant_pupadu <- subset(lfc_pupadu , padj < 0.0125 )
nrow(significant_pupadu)
write.table(significant_pupadu,"Results/PupaevsAdult_Significant_Genes.txt",row.names = F, col.names= F,sep = "\t",quote = F)
upregulated_pupadu <- subset(respupadu , padj < 0.0125 & log2FoldChange >  1.5)
nrow(upregulated_pupadu)
write.table(upregulated_pupadu,"Results/PupaevsAdult_Pupae_Upregulated_Genes.txt",row.names = F, col.names= F,sep = "\t",quote = F)
downregulated_pupadu<- subset(respupadu , padj < 0.0125 & log2FoldChange <  -1.5)
nrow(downregulated_pupadu)
write.table(downregulated_pupadu,"Results/PupaevsAdult_Adult_Upregulated_Genes.txt",row.names = F, col.names= F,sep = "\t",quote = F)

#######################Normalized counts ##############################################
normalized.counts <-fpm(ddsTC, robust = TRUE)
normalized.counts <- as.data.frame(normalized.counts)

Dnmt1a <- normalized.counts["Dnmt1a",]
Dnmt1a <- t(Dnmt1a)
Dnmt1a <- as.data.frame(Dnmt1a)
names(Dnmt1a)[1] <- "Counts"
Dnmt1a$Gene <- c("Dnmt1a", "Dnmt1a", "Dnmt1a", "Dnmt1a", "Dnmt1a" , "Dnmt1a","Dnmt1a", "Dnmt1a", 
                 "Dnmt1a", "Dnmt1a", "Dnmt1a" , "Dnmt1a", "Dnmt1a", "Dnmt1a" , "Dnmt1a")
Dnmt1a$Condition <- c("Embryo","Embryo","Embryo","Larva","Larva","Larva","Prepupa","Prepupa","Prepupa","Pupa","Pupa",
                      "Pupa","Adult","Adult","Adult")

LOC100115455 <- normalized.counts["LOC100115455",]
LOC100115455 <- t(LOC100115455)
LOC100115455 <- as.data.frame(LOC100115455)
names(LOC100115455)[1] <- "Counts"
LOC100115455$Gene <- c("Dnmt1b", "Dnmt1b", "Dnmt1b", "Dnmt1b", "Dnmt1b", "Dnmt1b", "Dnmt1b", "Dnmt1b",
                       "Dnmt1b", "Dnmt1b" , "Dnmt1b", "Dnmt1b", "Dnmt1b" , "Dnmt1b", "Dnmt1b")
LOC100115455$Condition <- c("Embryo","Embryo","Embryo","Larva","Larva","Larva","Prepupa","Prepupa","Prepupa","Pupa","Pupa",
                            "Pupa","Adult","Adult","Adult")

LOC100123657 <- normalized.counts["LOC100123657",]
LOC100123657 <- t(LOC100123657)
LOC100123657 <- as.data.frame(LOC100123657)
names(LOC100123657)[1] <- "Counts"
LOC100123657$Gene <- c("Dnmt1c", "Dnmt1c", "Dnmt1c","Dnmt1c", "Dnmt1c", "Dnmt1c",
                       "Dnmt1c", "Dnmt1c" , "Dnmt1c", "Dnmt1c", "Dnmt1c", "Dnmt1c",
                       "Dnmt1c", "Dnmt1c" , "Dnmt1c")
LOC100123657$Condition <- c("Embryo","Embryo","Embryo","Larva","Larva","Larva","Prepupa","Prepupa","Prepupa","Pupa","Pupa",
                            "Pupa","Adult","Adult","Adult")

LOC100114044 <- normalized.counts["LOC100114044",]
LOC100114044 <- t(LOC100114044)
LOC100114044 <- as.data.frame(LOC100114044)
names(LOC100114044)[1] <- "Counts"
LOC100114044$Gene <- c("Dnmt3", "Dnmt3", "Dnmt3","Dnmt3", "Dnmt3", "Dnmt3",
                       "Dnmt3", "Dnmt3" , "Dnmt3","Dnmt3", "Dnmt3" , "Dnmt3",
                       "Dnmt3", "Dnmt3" , "Dnmt3")
LOC100114044$Condition <- c("Embryo","Embryo","Embryo","Larva","Larva","Larva","Prepupa","Prepupa","Prepupa","Pupa","Pupa",
                            "Pupa","Adult","Adult","Adult")

LOC100114438<- normalized.counts["LOC100114438",]
LOC100114438 <- t(LOC100114438)
LOC100114438<- as.data.frame(LOC100114438)
names(LOC100114438)[1] <- "Counts"
LOC100114438$Gene <- c("Tet", "Tet","Tet", "Tet","Tet", "Tet","Tet", "Tet","Tet", "Tet","Tet", "Tet",
                       "Tet", "Tet","Tet")
LOC100114438$Condition <- c("Embryo","Embryo","Embryo","Larva","Larva","Larva","Prepupa","Prepupa","Prepupa","Pupa","Pupa",
                            "Pupa","Adult","Adult","Adult")


Mbd <- normalized.counts["Mbd",]
Mbd <- t(Mbd)
Mbd <- as.data.frame(Mbd)
names(Mbd)[1] <- "Counts"
Mbd$Gene <- c("Mbd", "Mbd", "Mbd","Mbd", "Mbd", "Mbd",
              "Mbd", "Mbd" , "Mbd","Mbd", "Mbd", "Mbd",
              "Mbd", "Mbd" , "Mbd")
Mbd$Condition <- c("Embryo","Embryo","Embryo","Larva","Larva","Larva","Prepupa","Prepupa","Prepupa","Pupa","Pupa",
                   "Pupa","Adult","Adult","Adult")
summary(Mbd) # 2812 average
mean(Dnmt1a$Counts)/mean(Mbd$Counts) # 0.2481033 Ratio

LOC100122899 <- normalized.counts["LOC100122899",]
LOC100122899 <- t(LOC100122899)
LOC100122899 <- as.data.frame(LOC100122899)
names(LOC100122899)[1] <- "Counts"
LOC100122899$Gene <- c("Tip60", "Tip60", "Tip60","Tip60", "Tip60", "Tip60",
                       "Tip60", "Tip60" , "Tip60","Tip60", "Tip60", "Tip60",
                       "Tip60", "Tip60" , "Tip60")
LOC100122899$Condition <- c("Embryo","Embryo","Embryo","Larva","Larva","Larva","Prepupa","Prepupa","Prepupa","Pupa","Pupa",
                            "Pupa","Adult","Adult","Adult")
summary(LOC100122899) # 1080.8 average

total_dnmts <- rbind(Dnmt1a, LOC100115455, LOC100123657, LOC100114044, LOC100114438,LOC100122899,Mbd)
total_dnmts$Condition <- factor(total_dnmts$Condition, levels=c("Embryo", "Larva", "Prepupa", "Pupa", "Adult"))
total_dnmts$sqrt <- sqrt(total_dnmts$Counts)

ggplot(total_dnmts, aes(x=Condition, y=sqrt, group=Gene)) +
  geom_line(aes(col=Gene), stat = "summary", fun = "mean") +
  xlab("Developmental Stage") +
  ylab("Sqrt FPM") +
  geom_point(aes(color=Gene, shape=Gene)) + 
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6)) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),  # X-axis label size and bold
    axis.title.y = element_text(size = 14, face = "bold")   # Y-axis label size and bold
  )

############################################# GO Terms ############################################################
geneID2GO <- readMappings(file = "~/Nasonia_vitripennis_HGD_go_annotation_condensed.txt", sep = "\t", IDsep = ",")
geneUniverse <- names(geneID2GO)

########################################## Embryo vs Larva  #############################################################

genes <- rownames(significant_emblar)

geneList <- factor(as.integer(geneUniverse %in% genes))
names(geneList) <- geneUniverse

# Biological Process

myGOdata <- new("topGOdata", description="My project", #labelling project
                ontology="BP", #This can be changed to CC or MF
                allGenes=geneList, #Include all genes in your list 
                annot = annFUN.gene2GO, #how to map genes
                gene2GO = geneID2GO) #genes are located in geneID2GO

sg <- sigGenes(myGOdata)
str(sg)
numSigGenes(myGOdata)
options(scipen = 999) # scientific notation messes this up. Turn it off until after topGO
resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
BPRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "weightFisher", ranksOf = "Fisher", topNodes = length(usedGO(myGOdata)))
length(usedGO(myGOdata))
as.numeric(BPRes$classicFisher)
sigBPRes <- BPRes[BPRes$classicFisher < 0.05 ,] 
write.table(sigBPRes, "~Results/GO_Terms/Embryo_vs_Larva_Biological_Process_GO_Terms.txt", quote = F, sep = "\t")

# Cellular Component

myGOdata <- new("topGOdata", description="My project", #labelling project
                ontology="CC", #This can be changed to CC or MF
                allGenes=geneList, #Include all genes in your list 
                annot = annFUN.gene2GO, #how to map genes
                gene2GO = geneID2GO) #genes are located in geneID2GO

sg <- sigGenes(myGOdata)
str(sg)
numSigGenes(myGOdata)
resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
CCRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "weightFisher", ranksOf = "Fisher", topNodes = length(usedGO(myGOdata)))
length(usedGO(myGOdata))
sigCCRes <- CCRes[CCRes$classicFisher < 0.05,] 
write.table(sigCCRes, "~Results/GO_Terms/Embryo_vs_Larva_Cellular_Component_GO_Terms.txt", quote = F, sep = "\t")

# Molecular Function

myGOdata <- new("topGOdata", description="My project", #labelling project
                ontology="MF", #This can be changed to CC or MF
                allGenes=geneList, #Include all genes in your list 
                annot = annFUN.gene2GO, #how to map genes
                gene2GO = geneID2GO) #genes are located in geneID2GO

sg <- sigGenes(myGOdata)
str(sg)
numSigGenes(myGOdata)
resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
MFRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "weightFisher", ranksOf = "Fisher", topNodes = length(usedGO(myGOdata)))
length(usedGO(myGOdata))
sigMFRes <- MFRes[MFRes$classicFisher < 0.05,] 
write.table(sigMFRes, "~Results/GO_Terms/Embryo_vs_Larva_Molecular_Function_GO_Terms.txt", quote = F, sep = "\t")

########################################## Larva vs Prepupa  #############################################################

genes <- rownames(significant_larpre)

geneList <- factor(as.integer(geneUniverse %in% genes))
names(geneList) <- geneUniverse

# Biological Process

myGOdata <- new("topGOdata", description="My project", #labelling project
                ontology="BP", #This can be changed to CC or MF
                allGenes=geneList, #Include all genes in your list 
                annot = annFUN.gene2GO, #how to map genes
                gene2GO = geneID2GO) #genes are located in geneID2GO

sg <- sigGenes(myGOdata)
str(sg)
numSigGenes(myGOdata)
options(scipen = 999) # scientific notation messes this up. Turn it off until after topGO
resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
BPRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "weightFisher", ranksOf = "Fisher", topNodes = length(usedGO(myGOdata)))
length(usedGO(myGOdata))
as.numeric(BPRes$classicFisher)
sigBPRes <- BPRes[BPRes$classicFisher < 0.05 ,] 
write.table(sigBPRes, "~Results/GO_Terms/Larva_vs_Prepupa_Biological_Process_GO_Terms.txt", quote = F, sep = "\t")

# Cellular Component

myGOdata <- new("topGOdata", description="My project", #labelling project
                ontology="CC", #This can be changed to CC or MF
                allGenes=geneList, #Include all genes in your list 
                annot = annFUN.gene2GO, #how to map genes
                gene2GO = geneID2GO) #genes are located in geneID2GO

sg <- sigGenes(myGOdata)
str(sg)
numSigGenes(myGOdata)
resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
CCRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "weightFisher", ranksOf = "Fisher", topNodes = length(usedGO(myGOdata)))
length(usedGO(myGOdata))
sigCCRes <- CCRes[CCRes$classicFisher < 0.05,] 
write.table(sigCCRes, "~Results/GO_Terms/Larva_vs_Prepupa_Cellular_Component_GO_Terms.txt", quote = F, sep = "\t")

# Molecular Function

myGOdata <- new("topGOdata", description="My project", #labelling project
                ontology="MF", #This can be changed to CC or MF
                allGenes=geneList, #Include all genes in your list 
                annot = annFUN.gene2GO, #how to map genes
                gene2GO = geneID2GO) #genes are located in geneID2GO

sg <- sigGenes(myGOdata)
str(sg)
numSigGenes(myGOdata)
resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
MFRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "weightFisher", ranksOf = "Fisher", topNodes = length(usedGO(myGOdata)))
length(usedGO(myGOdata))
sigMFRes <- MFRes[MFRes$classicFisher < 0.05,] 
write.table(sigMFRes, "~Results/GO_Terms/Larva_vs_Prepupa_Molecular_Function_GO_Terms.txt", quote = F, sep = "\t")

##########################################  Prepupa vs Pupa  #############################################################

genes <- rownames(significant_prepup)

geneList <- factor(as.integer(geneUniverse %in% genes))
names(geneList) <- geneUniverse

# Biological Process

myGOdata <- new("topGOdata", description="My project", #labelling project
                ontology="BP", #This can be changed to CC or MF
                allGenes=geneList, #Include all genes in your list 
                annot = annFUN.gene2GO, #how to map genes
                gene2GO = geneID2GO) #genes are located in geneID2GO

sg <- sigGenes(myGOdata)
str(sg)
numSigGenes(myGOdata)
options(scipen = 999) # scientific notation messes this up. Turn it off until after topGO
resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
BPRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "weightFisher", ranksOf = "Fisher", topNodes = length(usedGO(myGOdata)))
length(usedGO(myGOdata))
as.numeric(BPRes$classicFisher)
sigBPRes <- BPRes[BPRes$classicFisher < 0.05 ,] 
write.table(sigBPRes, "~Results/GO_Terms/Prepupa_vs_Pupa_Biological_Process_GO_Terms.txt", quote = F, sep = "\t")

# Cellular Component

myGOdata <- new("topGOdata", description="My project", #labelling project
                ontology="CC", #This can be changed to CC or MF
                allGenes=geneList, #Include all genes in your list 
                annot = annFUN.gene2GO, #how to map genes
                gene2GO = geneID2GO) #genes are located in geneID2GO

sg <- sigGenes(myGOdata)
str(sg)
numSigGenes(myGOdata)
resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
CCRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "weightFisher", ranksOf = "Fisher", topNodes = length(usedGO(myGOdata)))
length(usedGO(myGOdata))
sigCCRes <- CCRes[CCRes$classicFisher < 0.05,] 
write.table(sigCCRes, "~Results/GO_Terms/Prepupa_vs_Pupa_Cellular_Component_GO_Terms.txt", quote = F, sep = "\t")

# Molecular Function

myGOdata <- new("topGOdata", description="My project", #labelling project
                ontology="MF", #This can be changed to CC or MF
                allGenes=geneList, #Include all genes in your list 
                annot = annFUN.gene2GO, #how to map genes
                gene2GO = geneID2GO) #genes are located in geneID2GO

sg <- sigGenes(myGOdata)
str(sg)
numSigGenes(myGOdata)
resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
MFRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "weightFisher", ranksOf = "Fisher", topNodes = length(usedGO(myGOdata)))
length(usedGO(myGOdata))
sigMFRes <- MFRes[MFRes$classicFisher < 0.05,] 
write.table(sigMFRes, "~Results/GO_Terms/Prepupa_vs_Pupa_Molecular_Function_GO_Terms.txt", quote = F, sep = "\t")

##########################################  Pupa vs Adult  #############################################################

genes <- rownames(significant_pupadu)

geneList <- factor(as.integer(geneUniverse %in% genes))
names(geneList) <- geneUniverse

# Biological Process

myGOdata <- new("topGOdata", description="My project", #labelling project
                ontology="BP", #This can be changed to CC or MF
                allGenes=geneList, #Include all genes in your list 
                annot = annFUN.gene2GO, #how to map genes
                gene2GO = geneID2GO) #genes are located in geneID2GO

sg <- sigGenes(myGOdata)
str(sg)
numSigGenes(myGOdata)
options(scipen = 999) # scientific notation messes this up. Turn it off until after topGO
resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
BPRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "weightFisher", ranksOf = "Fisher", topNodes = length(usedGO(myGOdata)))
length(usedGO(myGOdata))
as.numeric(BPRes$classicFisher)
sigBPRes <- BPRes[BPRes$classicFisher < 0.05 ,] 
write.table(sigBPRes, "~Results/GO_Terms/Pupa_vs_Adult_Biological_Process_GO_Terms.txt", quote = F, sep = "\t")

# Cellular Component

myGOdata <- new("topGOdata", description="My project", #labelling project
                ontology="CC", #This can be changed to CC or MF
                allGenes=geneList, #Include all genes in your list 
                annot = annFUN.gene2GO, #how to map genes
                gene2GO = geneID2GO) #genes are located in geneID2GO

sg <- sigGenes(myGOdata)
str(sg)
numSigGenes(myGOdata)
resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
CCRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "weightFisher", ranksOf = "Fisher", topNodes = length(usedGO(myGOdata)))
length(usedGO(myGOdata))
sigCCRes <- CCRes[CCRes$classicFisher < 0.05,] 
write.table(sigCCRes, "~Results/GO_Terms/Pupa_vs_Adult_Cellular_Component_GO_Terms.txt", quote = F, sep = "\t")

# Molecular Function

myGOdata <- new("topGOdata", description="My project", #labelling project
                ontology="MF", #This can be changed to CC or MF
                allGenes=geneList, #Include all genes in your list 
                annot = annFUN.gene2GO, #how to map genes
                gene2GO = geneID2GO) #genes are located in geneID2GO

sg <- sigGenes(myGOdata)
str(sg)
numSigGenes(myGOdata)
resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
MFRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "weightFisher", ranksOf = "Fisher", topNodes = length(usedGO(myGOdata)))
length(usedGO(myGOdata))
sigMFRes <- MFRes[MFRes$classicFisher < 0.05,] 
write.table(sigMFRes, "~Results/GO_Terms/Pupa_vs_Adult_Molecular_Function_GO_Terms.txt", quote = F, sep = "\t")
