library(GenomicFeatures)
library(Rsamtools)
library(GenomicAlignments)
library(DEXSeq)

txdb = makeTxDbFromGFF("GCF_009193385.2_Nvit_psr_1.1_genomic.gtf.gz")
flattenedAnnotation = exonicParts( txdb, linked.to.single.gene.only=TRUE )
names(flattenedAnnotation) =
    sprintf("%s:E%0.3d", flattenedAnnotation$gene_id, flattenedAnnotation$exonic_part)

bamFiles = c("Embryo1.bam","Embryo2.bam","Larva1.bam","Larva2.bam","Larva3.bam","Prepupa1","Prepupa2","Prepupa3","Pupa1","Pupa2","Pupa3","Adult1","Adult2","Adult3")
bamFiles = BamFileList( bamFiles )
#seqlevelsStyle(flattenedAnnotation) = "NBCI"
se = summarizeOverlaps(
    flattenedAnnotation, BamFileList(bamFiles), singleEnd=FALSE,
    fragments=TRUE, ignore.strand=TRUE )

save.image(file = "Initial_Load.RData")

colData(se)$condition = factor(c("Adult","Adult","Adult","Embryo","Embryo","Embryo","Larva","Larva","Larva",
"Prepupae","Prepupae","Prepupae","Pupae","Pupae","Pupae"))


dxd = DEXSeqDataSetFromSE( se, design= ~ sample + exon + condition:exon )
colData(dxd)
head( counts(dxd), 5 )
sampleAnnotation( dxd )

#DEXSeq
dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd )
normalized.counts <- counts(dxd, normalized=TRUE)
write.table(as.data.frame(normalized.counts),file="Normalised_Exons.txt", quote = FALSE, sep = "\t")
