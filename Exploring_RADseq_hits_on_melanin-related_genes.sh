######################################
###Sequence of steps to ######
######################################
#populations.snps.vcf is the final snp panel 
#
#loop to inspect annotated file
list="g12671.t1
g11537.t1
g14279.t1
g17775.t1
g11868.t1
g5619.t1
g7780.t1
g815.t1
g11540.t1
g6385.t1
g7206.t1
g2824.t1
g13007.t1
g16980.t1
g16384.t1"

for item in $list
do vcftools --vcf populations.snps.vcf --chr $item --out ${item}_name --recode;done

##loop to inspect in which contig annotated melanin regions lie
list="g12671.t1
g11537.t1
g14279.t1
g17775.t1
g11868.t1
g5619.t1
g7780.t1
g815.t1
g11540.t1
g6385.t1
g7206.t1
g2824.t1
g13007.t1
g16980.t1"

for item in $list
do grep $item v5_cds_fasta.fasta >> test; done

##loop to inspect how many radseq hits existed in the region 

list="contig_48|arrow
contig_398|arrow
contig_602|arrow
contig_876|arrow
contig_398|arrow
contig_1954|arrow
contig_2965|arrow
contig_1022|arrow
contig_398|arrow
contig_2237|arrow
contig_2850|arrow
contig_1326|arrow
contig_513|arrow
contig_793|arrow"

for item in $list
do vcftools --vcf populations.snps.vcf --chr $item --out melanin/${item}_name --recode;done

##retrieve contig sizes

list="contig_48|arrow
contig_398|arrow
contig_602|arrow
contig_876|arrow
contig_398|arrow
contig_1954|arrow
contig_2965|arrow
contig_1022|arrow
contig_398|arrow
contig_2237|arrow
contig_2850|arrow
contig_1326|arrow
contig_513|arrow
contig_793|arrow"

for item in $list
do grep $item populations.snps.vcf >> contig_sizes; done

##extract contigs hits from GEMMA
list="chr
contig_48|arrow
contig_398|arrow
contig_602|arrow
contig_876|arrow
contig_398|arrow
contig_1954|arrow
contig_2965|arrow
contig_1022|arrow
contig_398|arrow
contig_2237|arrow
contig_2850|arrow
contig_1326|arrow
contig_513|arrow
contig_793|arrow"

for item in $list
do grep $item R08.assoc.txt >> melanin_contigs; done

###extract mapped statistics per contig####
##first index bamfiles##
samtools idxstats input_alignments_sorted.bam
##extract statistics from bamfiles - it take from all samples##
for i in *.bam
do samtools idxstats ${i} >>melanin_contigs; done
##extract contigs from bamfiles##

for item in $list
do grep $item *.bam >> melanin_contigs; done

list="chr
contig_48|arrow
contig_398|arrow
contig_602|arrow
contig_876|arrow
contig_398|arrow
contig_1954|arrow
contig_2965|arrow
contig_1022|arrow
contig_398|arrow
contig_2237|arrow
contig_2850|arrow
contig_1326|arrow
contig_513|arrow
contig_793|arrow"

for item in $list
do grep $item melanin_contigs >> melanin_contigs_sorted; done


### average mapped to each melanin contig###
for i in *.name
do awk '{ total += $3 } END { print total/NR }' ${i} > ${i}.avg;
done

## regions within contigs
 
contig_48|arrow 22070833:22098051
contig_398|arrow 48035584:48036561
contig_602|arrow 32715607:32729374
contig_876|arrow 6062406:6063383
contig_398|arrow 79816315:79817310
contig_1954|arrow 1382950:1383894
contig_2965|arrow 9775309:9776649
contig_1022|arrow 5236575:5237963
contig_398|arrow 48085765:48086724
contig_2237|arrow 2416945:2418078
contig_2850|arrow 5475376:5476356
contig_1326|arrow 11282764:11291334
contig_513|arrow 283407:308956
contig_793|arrow 3815682:3825591

for i in *.bam
do samtools view -b ${i} "contig_48|arrow:22070833-22098051" | samtools view -c  > g12671; done
#20kb
for i in *.bam
do samtools view -b ${i} "contig_48|arrow:22050833-22228051" | samtools view -c  > g12671_20kb; done
##coverage
samtools coverage -b test -r 'contig_48|arrow:22070833-22098051' -o cov_g12671
samtools coverage -b test -r 'contig_48|arrow:22050833-22228051' -o cov_g12671_20kb
#########################
#########################
#########################
for i in *.bam
do samtools view -b ${i} "contig_398|arrow:48035584-48036561" | samtools view -c > g11537; done
#20kb
for i in *.bam
do samtools view -b ${i} "contig_398|arrow:48015584-48056561" | samtools view -c > g11537_20kb; done
##coverage
samtools coverage -b test -r 'contig_398|arrow:48035584-48036561' -o cov_g11537
samtools coverage -b test -r 'contig_398|arrow:48015584-48056561' -o cov_g11537_20kb
#########################
#########################
#########################
for i in *.bam
do samtools view -b ${i} "contig_602|arrow:32715607-32729374" | samtools view -c  > g14279; done
#20kb
for i in *.bam
do samtools view -b ${i} "contig_602|arrow:32695607-32929374" | samtools view -c  > g14279_20kb; done
##coverage
samtools coverage -b test -r 'contig_602|arrow:32715607-32729374' -o cov_g14279
samtools coverage -b test -r 'contig_602|arrow:32695607-32929374' -o cov_g14279_20kb
#########################
#########################
#########################
for i in *.bam
do samtools view -b ${i} "contig_876|arrow:6062406-6063383" | samtools view -c > g17775; done
#20kb
for i in *.bam
do samtools view -b ${i} "contig_876|arrow:6042406-6083383" | samtools view -c > g17775_20kb; done
##coverage
samtools coverage -b test -r 'contig_876|arrow:6062406-6063383' -o cov_g17775
samtools coverage -b test -r 'contig_876|arrow:6042406-6083383' -o cov_g17775_20kb
#########################
#########################
#########################
for i in *.bam
do samtools view -b ${i} "contig_398|arrow:79816315-79817310" | samtools view -c > g11868; done
#20kb
for i in *.bam
do samtools view -b ${i} "contig_398|arrow:79796315-79837310" | samtools view -c > g11868_20kb; done
##coverage
samtools coverage -b test -r 'contig_398|arrow:79816315-79817310' -o cov_g11868
samtools coverage -b test -r 'contig_398|arrow:79796315-79837310' -o cov_g11868_20kb
#########################
#########################
#########################
for i in *.bam
do samtools view -b ${i} "contig_1954|arrow:1382950-1383894" | samtools view -c > g5619; done
#20kb
for i in *.bam
do samtools view -b ${i} "contig_1954|arrow:1362950-1403894" | samtools view -c > g5619_20kb; done
##coverage
samtools coverage -b test -r 'contig_1954|arrow:1382950-1383894' -o cov_g5619
samtools coverage -b test -r 'contig_1954|arrow:1362950-1403894' -o cov_g5619_20kb

#########################
#########################
#########################
for i in *.bam
do samtools view -b ${i} "contig_2965|arrow:9775309-9776649" | samtools view -c > g7780; done
#20kb
for i in *.bam
do samtools view -b ${i} "contig_2965|arrow:9755309-9796649" | samtools view -c > g7780_20kb; done
##coverage
samtools coverage -b test -r 'contig_2965|arrow:9775309-9776649' -o cov_g7780
samtools coverage -b test -r 'contig_2965|arrow:9755309-9796649' -o cov_g7780_20kb
#########################
#########################
#########################
for i in *.bam
do samtools view -b ${i} "contig_1022|arrow:5236575-5237963" | samtools view -c > g815; done
#20kb
for i in *.bam
do samtools view -b ${i} "contig_1022|arrow:5216575-5257963" | samtools view -c > g815_20kb; done
##coverage
samtools coverage -b test -r 'contig_1022|arrow:5236575-5237963' -o cov_g815
samtools coverage -b test -r 'contig_1022|arrow:5216575-5257963' -o cov_g815_20kb
#########################
#########################
#########################
for i in *.bam
do samtools view -b ${i} "contig_398|arrow:48085765-48086724" | samtools view -c > g11540; done
#20kb
for i in *.bam
do samtools view -b ${i} "contig_398|arrow:48065765-48106724" | samtools view -c > g11540_20kb; done
#coverage
samtools coverage -b test -r 'contig_398|arrow:48085765-48086724' -o cov_g11540
samtools coverage -b test -r 'contig_398|arrow:48065765-48106724' -o cov_g11540_20kb
#########################
#########################
#########################
for i in *.bam
do samtools view -b ${i} "contig_2237|arrow:2416945-2418078" | samtools view -c > g6385; done
#20kb
for i in *.bam
do samtools view -b ${i} "contig_2237|arrow:2396945-2438078" | samtools view -c > g6385_20kb; done
##coverage
samtools coverage -b test -r 'contig_2237|arrow:2416945-2418078' -o cov_g6385
samtools coverage -b test -r 'contig_2237|arrow:2396945-2438078' -o cov_g6385_20kb
#########################
#########################
#########################
for i in *.bam
do samtools view -b ${i} "contig_2850|arrow:5475376-5476356" | samtools view -c > g7206; done
#20kb
for i in *.bam
do samtools view -b ${i} "contig_2850|arrow:5455376-5496356" | samtools view -c > g7206_20kb; done
#coverage
samtools coverage -b test -r 'contig_2850|arrow:5475376-5476356' -o cov_g7206
samtools coverage -b test -r 'contig_2850|arrow:5455376-5496356' -o cov_g7206_20kb
#########################
#########################
#########################
for i in *.bam
do samtools view -b ${i} "contig_1326|arrow:11282764-11291334" | samtools view -c > g2824; done
#20kb
for i in *.bam
do samtools view -b ${i} "contig_1326|arrow:11262764-11311334" | samtools view -c > g2824_20kb; done
#coverage
samtools coverage -b test -r 'contig_1326|arrow:11282764-11291334' -o cov_g2824
samtools coverage -b test -r 'contig_1326|arrow:11262764-11311334' -o cov_g2824_20kb
#########################
#########################
#########################
for i in *.bam
do samtools view -b ${i} "contig_513|arrow:283407-308956" | samtools view -c > g13007; done
#20kb
for i in *.bam
do samtools view -b ${i} "contig_513|arrow:263407-328956" | samtools view -c > g13007_20kb; done
#coverage
samtools coverage -b test -r 'contig_513|arrow:283407-308956' -o cov_g13007
samtools coverage -b test -r 'contig_513|arrow:263407-328956' -o cov_g13007_20kb
#########################
#########################
#########################
for i in *.bam
do samtools view -b ${i} "contig_793|arrow:3815682-3825591" | samtools view -c > g16980; done
#20kb
for i in *.bam
do samtools view -b ${i} "contig_793|arrow:3795682-3845591" | samtools view -c > g16980_20kb; done
#coverage
samtools coverage -b test -r 'contig_793|arrow:3815682-3825591' -o cov_g16980
samtools coverage -b test -r 'contig_793|arrow:3795682-3845591' -o cov_g16980_20kb

######################################
###COMPARE V5 AND V6 FOR MELANIN######
######################################

#USE BLASR TO FIND CDS OF V6 IN V5
blasr v6_cds.fasta candidate_index.fasta -m 1 --hitPolicy all > candidate_v6
awk '{print$2}' candidate_v6 > v6_list.txt
#USE SEQTK TO SAMPLE OUT FROM THE WHOLE V6 CDS# grep also works
PATH=$PATH:/users/baltazar/software/seqtk
seqtk subseq v6_cds.fasta v6_list.txt > v6_candidates.fasta
