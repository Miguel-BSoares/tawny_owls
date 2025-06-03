#####
# Methodology to trim, align and count mRNA reads
# Using <reference> for the genome or cds used as maping reference

module load trimmomatic
mkdir -p test_trim2
for R1 in *_1.fq.gz; do \
        R2="${R1%_1.fq.gz}_2.fq.gz" \
        trimmomatic PE -phred33 $R1 $R2 "${R1%.*}_PE.fq" "${R1%.*}_SR.fq" "${R2%.*}_PE.fq" "${R2%.*}_SR.fq" \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 ILLUMINACLIP:TruSeq2-PE.fa:2:30:10 SLIDINGWINDOW:5:20 MINLEN:36  \  
        2> test_trim2/test.log; \
done

## Construction of genome index ##
hisat2-build <reference>.fasta genome_index

## extracting names from files and mapping loop ##
ls *2.fq_PE.fq | paste -sd'\n' | sed 's%_2.fq_PE.fq% %g'>files

## add files content to bash script and map against genome/cds ##
files=""
for sample in $files
do
		 hisat2 -p 10 -x hisat_index_cds/cds_index -1 ../test_trim2/${f}_1.fq_PE.fq -2 ../test_trim2/${f}_2.fq_PE.fq |
         samtools view -b |
         samtools sort --threads 10 > aligned/${sample}.bam
done

## counts ##
# require an annotation file #

module load htseq
module load samtools
module load igv

mkdir -p ht_seq_counts_nostrand

for i in *.bam

        do htseq-count -f bam ${i} -r name v5_annotation.gff --type gene --idattr ID --stranded no >> ht_seq_counts_nostrand/htseq_count${i}.txt;

done
