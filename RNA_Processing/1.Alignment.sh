#Align using STAR package (Note you will need to create a star index first)

REF=~/star_index/

for file in $(ls *1.fq.gz)
do
    base=$(basename $file "1.fq.gz")
    STAR \
    --runThreadN 8 \
    --twopassMode Basic \
    --outSAMtype BAM SortedByCoordinate \
    --limitBAMsortRAM 9203070763 \
    --genomeDir ${REF} \
    --readFilesCommand zcat \
    --readFilesIn ${base}1.fq.gz ${base}2.fq.gz \
    --outFileNamePrefix ~/${base}
done
