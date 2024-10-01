#When looking at the each developmental stage on the whole, I combined each of the fastq files for each stage together using:

Embryo1_1.fq.gz Embryo2_1.fq > Embryo_1.fq
Embryo1_2.fq.gz Embryo2_2.fq > Embryo_2.fq

#Align to reference genome (Note you will need to run bismark genome prepare first)

REF_FA=~/Nvit_psr_1.1/

for file in $(ls *1.fq.gz)
do
  	base=$(basename $file "1.fq.gz")
        ~/Bismark-0.24.2/bismark ${REF_FA} -1 ${base}1.fq.gz -2 ${base}2.fq.gz
done

#Align to lambda genome

REF_FA=~/Lambda/

for file in $(ls *1.fastq.gz)
do
  	base=$(basename $file "1.fastq.gz")
        ~/Bismark-0.24.2/bismark ${REF_FA} -1 ${base}1.fq.gz -2 ${base}2.fq.gz
done

#Deduplicate the aligned bams

for file in $(ls *bam)
do
  	base=$(basename $file ".bam")
       ~/Bismark-0.24.2/deduplicate_bismark -p --bam ${file}
done

#Extract the methylation status of CpGs

for file in $(ls *.bam)
do
  	base=$(basename $file ".bam")
        ~/Bismark-0.24.2/bismark_methylation_extractor -p --no_overlap  --comprehensive --bedgraph \
        --genome_folder ~/Nvit_psr_1.1/ \
        --report --cytosine_report --buffer_size 50% \
${file}  -o ${base}_Extract
done

#Destrand the samples

for file in $(ls *cov.gz)
do
  	base=$(basename $file ".cov.gz")
        ~/Bismark-0.24.2/coverage2cytosine -o ${base}_mergeCpGs --merge_CpGs \
        --genome_folder ~/Nvit_psr_1.1/ \
${file}
done


