#Count Reads using HTSeq

for file in $(ls ~/*.bam)
do
    base=$(basename ${file} ".bam")
    htseq-count --idattr=Name --additional-attr=gene --format=bam --type=gene \
    --mode intersection-strict \
    ${file} \
~/GCF_009193385.2_Nvit_psr_1.1_genomic.gff \
-s reverse \
> counts/${base}.counts.txt
done
