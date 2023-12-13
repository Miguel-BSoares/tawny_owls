############################################################################
###----------------Aligment of RADseq reads to tawny owl genome----------###
############################################################################
####aligning of RADseq reads to tawny owl genome
####software utilized:
####gcc/9.1.0
####bwa-mem2/2.2
####samtools/1.7
####stacks/2.65
#alignment
files="JUST PREFIX FILE NAME"
for sample in $files
do
        bwa-mem2 mem -t 10 genome/v6.fasta ${sample}_R1_001.fastq.gz ${sample}_R2_001.fastq.gz |
         samtools view -b |
         samtools sort --threads 10 > aligned/${sample}.bam
done

#Stacks
gstacks -I ./  -M popmap_final -O 01_colour -t 20
populations -P ./  -M popmap_final -O 02_populations -r 0.8 --min-maf 0.05 --ordered-export --write-random-snp -t 20

##----Replicates of Fst with 100 randomly picked loci----####

for (( i = 1; i<= 100;i++ ));do
        grep -v "^#" populations.sumstats.tsv | cut -f 2,3 | sort | uniq | shuf | head -n 2 | sort -n > whitelist/w$i.tsv
        done
##2 extract from original vcf file and create new ones to be passed to populations
for file in whitelist/*.tsv; do vcftools --vcf populations.snps.vcf --positions whitelist/$file --out whitelist/${file} --recode;done
##3 perform populations across all vcf files
for file in whitelist/*.vcf; do populations -V $file -M popmap_final_years -O 02_resampling/ --fstats -t 10; done

####----Plotting Fst histogram of 100 randomly picked loci per pairwise year comparison----####
library(ggplot2)
p2009_2011<-read.csv("2009_2011.csv", sep = ";", header = TRUE)
m<-mean(p2009_2011$fst)
sd<-sd(p2009_2011$fst)
a<-p2009_2011$fst
ci<-m+(1.96*(sd/sqrt(100)))
p<-ggplot(p2009_2011, aes(x=fst)) + geom_histogram(color="black", fill="darkgrey")+
geom_vline(aes(xintercept=ci),color="red", linetype="dashed", size=1)+
geom_vline(aes(xintercept=0.02337),color="black", linetype="dashed", size=1)+
geom_vline(aes(xintercept=0.015316),color="green", linetype="dashed", size=1)+
theme_classic()+  
geom_density(alpha=.1, fill="yellow") 
p

############################################################################
###----------------Extracting zigosity from vcf----------###
############################################################################
#extract specific positions
vcftools --vcf populations.snps.vcf --positions gwas --out gwas_snps --recode
##extract specific chromossomes --chr can be used multiple times 
vcftools --vcf populations.snps.vcf --chr contig_602\|arrow --out contig_602 --recode
##normalize
bcftools norm -m-any populations.snps.vcf -Ov > pop.norm.vcf
##remove IDs with few sites and sites uncalled
vcftools --vcf contigs_3657_602.recode.vcf --keep final_individuals --max-missing 1 --recode --out contigs_filtered
##make genotype matrix
vcftools --vcf contigs_3657_602.recode.vcf --012 --out contigs_matrix

##count heterozigotes/homozygotes from vcf file
paste <(bcftools view contigs_3657_602.norm.vcf |\
awk -F"\t" 'BEGIN {print "CHR\tPOS\tID\tREF\tALT"} \
      !/^#/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5}') \
    \
  <(bcftools query -f '[\t%SAMPLE=%GT]\n' contigs_3657_602.norm.vcf |\
    awk 'BEGIN {print "nHet"} {print gsub(/0\|1|1\|0|0\/1|1\/0/, "")}') \
    \
  <(bcftools query -f '[\t%SAMPLE=%GT]\n' contigs_3657_602.norm.vcf |\
    awk 'BEGIN {print "nHomAlt"} {print gsub(/1\|1|1\/1/, "")}') \
    \
  <(bcftools query -f '[\t%SAMPLE=%GT]\n' contigs_3657_602.norm.vcf |\
    awk 'BEGIN {print "nHomRef"} {print gsub(/0\|0|0\/0/, "")}') \
    \
  <(bcftools view contigs_3657_602.norm.vcf | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HetSamples"} \
    !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\|1|1\|0|0\/1|1\/0/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
    \
  <(bcftools view contigs_3657_602.norm.vcf | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HomSamplesAlt"} \
    !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/1\|1|1\/1/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
    \
  <(bcftools view contigs_3657_602.norm.vcf | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HomSamplesRef"} \
    !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\|0|0\/0/,"", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
    \
  | sed 's/,\t/\t/g' | sed 's/,$//g'

#coding
#heterozygotes=0
#homozygotes_ref=1
#homozygotes_alt=2

####################################
##------------adegenet ----------###
####################################
library(adegenet)
library(RColorBrewer)
library(ggplot2)
colour2<-c("gainsboro","darkgoldenrod3")
n_color<-2
pallete20<-rainbow(2)
colour3<-c("gray78","sienna","dodgerblue")
color<-read.genetix("color.gtx")
color$pop
#find K
color_clusters<-find.clusters(color,max.n.clust = 20,choose=FALSE, n.iter = 2000 )
plot(color_clusters$Kstat, type="o", xlab="number of clusters (K)", ylab="BIC",
     col="blue", main="Detection based on BIC")
points(2, color_clusters$Kstat[2], pch="x", cex=3)
mtext(3, tex="'X' indicates the actual number of clusters")

#choose K based on previous step
color_clusters<-find.clusters(color,max.n.clust = 20)
color_dapc<-dapc.genind(color,color_clusters$grp)###defined by find clusters
pops<-color_clusters$grp##pops defined by clusters
write.csv(color_clusters$grp,file = "groups")
pops_origin<-color$pop##original pops
plot2<-cbind(pops,color_dapc$ind.coord)
write.csv(plot2,file = "dapc_coords")
#dapc_coords are the covariates for PLINK#
color_dapc

#plot
plot_ne<-ggplot(plot2,aes(LD1,LD2, colour=color$pop))+
  geom_point(colour = "gray50", size = 4)+geom_point(size=3)+
  scale_colour_manual(values = colour3,name="Original groups",labels=c("Gray", "Brown", "Unknown"))+
  theme_bw()+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot_ne

## plot for loadings =1
d<-density(color_dapc$ind.coord)
plot(d, main =" Distribution LD 1 for K = 2", xlab="coordinates of LD1")
polygon(d, col="gainsboro", border="black") 

################################
##------------SNPrelate----------###
################################
#plotting and exploring PCA
library(gdsfmt)
library(SNPRelate)
library(scatterplot3d)
library(RColorBrewer)
library(snpStats)
library(rtracklayer)
library(biomaRt)
vcf.fn <- "populations.snps.vcf"
snpgdsVCF2GDS(vcf.fn, "final2.gds", method="biallelic.only")
snpgdsSummary("final2.gds")
genofile <- openfn.gds("final2.gds")
head(genofile)
####LD prunning - find threshold of linkage desiquilibiurm##
ld_02 <- snpgdsLDpruning(genofile, ld.threshold=0.2,autosome.only=FALSE,missing.rate=0.2)
ld_05 <- snpgdsLDpruning(genofile, ld.threshold=0.5,autosome.only=FALSE,missing.rate=0.20)
snpset.id <- unlist(unname(ld_05))
####pcas and trees####
bubu_pca <- snpgdsPCA(genofile,autosome.only=FALSE,snp.id=snpset.id)
pc.percent <- bubu_pca$varprop*100
pcas<-head(round(pc.percent,2))
plot(pcas,cex=1,col="blue",pch=20, main="PC contribution", xlab="Principal Component",ylab="% variation")
tab<-data.frame(sample.id=bubu_pca$sample.id,
                EV1=bubu_pca$eigenvect[,1],EV2=bubu_pca$eigenvect[,2],EV3=bubu_pca$eigenvect[,3],
                strigsAsFactors=FALSE) 
#draw
par(mar=c(5,5,4,4))
plot(tab$EV2,tab$EV1,pch=20)
#3d plot
library(scatterplot3d)
# shading relation to y axis
scatterplot3d(tab$EV2,tab$EV1,tab$EV3,pch=20,grid=FALSE,tick.marks=TRUE,type="p",highlight.3d=TRUE,angle=45)
ibs.hc <- snpgdsHCluster(snpgdsIBS(genofile, num.thread=2,autosome.only=FALSE),snp.id=snpset.id)
rv2 <- snpgdsCutTree(ibs.hc)
rv <- snpgdsCutTree(ibs.hc)
plot(rv$dendrogram, leaflab="none")
##########GROUPS#######################
sample.id <- scan("id.txt", what=character())
col_code <- scan("col.txt", what=character())
fam_code <- scan("fam.txt", what=character())
year_code <- scan("year.txt", what=character())
tail(cbind(sample.id,col_code,fam_code,year_code))
full_file<-cbind(sample.id,col_code,fam_code,year_code)
##colors#
col_years <- brewer.pal(11, "BrBG")
col_fams <- 
#morphs
tab<-data.frame(sample.id=bubu_pca$sample.id,
                col=factor(col_code)[match(bubu_pca$sample.id,sample.id)],
                EV1=bubu_pca$eigenvect[,1],EV2=bubu_pca$eigenvect[,2],EV3=bubu_pca$eigenvect[,3],
                strigsAsFactors=FALSE) 
par(mar=c(5,5,4,4))
plot(tab$EV1,tab$EV2,col=as.integer(tab$col),pch=20,ylab="PC2(2.83%)",xlab="PC1(3.30%)", main="22093 SNPs, 370 owls")
legend("bottomleft", inset=0.05, bty="n", cex=0.75, title="Morphs", legend=levels(tab$col), col=1:nlevels(tab$col),fill = 1:nlevels(tab$col))
#family##plotting can be adjusted
tab<-data.frame(sample.id=bubu_pca$sample.id,
                fam=factor(fam_code)[match(bubu_pca$sample.id,sample.id)],
                EV1=bubu_pca$eigenvect[,1],EV2=bubu_pca$eigenvect[,2],EV3=bubu_pca$eigenvect[,3],
                strigsAsFactors=FALSE) 
par(mar=c(5,5,4,4))
plot(tab$EV1,tab$EV2,col=as.integer(tab$fam),pch=20,ylab="PC2(2.83%)",xlab="PC1(3.30%)", main="22093 SNPs, 370 owls")
legend("bottomleft", inset=0.05, bty="n", cex=0.75, title="Families", legend=levels(tab$fam), 
       col=1:nlevels(tab$fam),fill = 1:nlevels(tab$fam))
#year
tab<-data.frame(sample.id=bubu_pca$sample.id,
                year=factor(year_code)[match(bubu_pca$sample.id,
                sample.id)],EV1=bubu_pca$eigenvect[,1],EV2=bubu_pca$eigenvect[,2],EV3=bubu_pca$eigenvect[,3],
                strigsAsFactors=FALSE) 
par(mar=c(5,5,4,4))
plot(tab$EV1,tab$EV2,col=col_years,pch=20,ylab="PC2(2.83%)",xlab="PC1(3.30%)", main="22093 SNPs, 370 owls")
legend("bottomleft", inset=0.05, bty="n", cex=0.75, title="Year", legend=levels(tab$year), col=col_years,fill = col_years)

##
##########################################################################################################
####################---------------Genome wide association analyses----------#############################
##########################################################################################################

################################
##------------Plink----------###
################################

plink2 -bfile genome_maf005_noR --hwe 1e-6 --geno 0.2 --maf 0.05 --make-bed --out genome_maf005_variant_filtered --allow-extra-chr
plink2 --bfile genome_maf005_variant_filtered --indep-pairwise 50 10 0.1 --allow-extra-chr --out genome_maf005_filtered_hwe_kingoff_ld
plink2 --bfile genome_maf005_variant_filtered --extract genome_maf005_filtered_hwe_kingoff_ld.prune.in --out genome_maf005_fullfiltered --allow-extra-chr --make-bed
plink2 --bfile genome_maf005_fullfiltered --out report_genome_maf005_fullfiltered --allow-extra-chr --nonfounders --glm --covar dapc_coords -parameters 1-2 --pfilter 0.001

################################
##------------GEMMA----------###
################################

Gemma –bfile [plink file] –gk [1 or 2] –o  [prefix for the output]
Gemma –bfile [plink file] –k [relatedness file] –lmm [num]–o  [prefix for the output] –c [covariates] –n [phenotype column number]
                                                               
#explore GEMMA results
assoc2<-read.csv("R08_dapc2.assoc.txt", sep="\t",dec=".", header = TRUE)
assoc2
n<-12673 ##Number snps
sig_level=-log10(0.05/n)#signficant level to be applied as a function of the number of SNPs
sug_level=-log10(1/n)#suggestive level to be applied as a function of the number of SNPs
assoc2$log_p_lrt<-(-log10(assoc2$p_lrt))##log transforms the p
assoc_lrt_sugg<-subset(assoc2,log_p_lrt>=4)
assoc_lrt_sugg
write.csv(assoc_lrt_sugg, file = "results_dapc2")

################################
##------------Wtest----------###
################################

library(wtest)
color_allele<-read.csv("wtest_allele.csv",header = TRUE,sep = ";")
color_pheno<-read.csv("wtest_pheno.csv",header = FALSE, sep = ";")
###main effect calculation
hf1 <- hf(data = color_test, w.order = 1, B = 100)
w1<-wtest(data = color_test, y = color_pheno, w.order = 1, hf1 = hf1 )
w1$results
w.diagnosis(data = color_test, w.order = 1, n.rep = 10,hf1 = hf1,hf2 = hf2)
w.qqplot(data = color_test, y = color_pheno,w.order = 1, input.poolsize = 200,hf1 = hf1,hf2 = hf2)

###main effect calculation alleles
hf1 <- hf(data = color_allele, w.order = 1, B = 100)
w1<-wtest(data = color_allele, y = color_pheno, w.order = 1, hf1 = hf1 )
w.diagnosis(data = color_allele, w.order = 1, n.rep = 1000,hf1 = hf1)
w.qqplot(data = color_allele, y = color_pheno, w.order = 1, input.poolsize = 200,hf1 = hf1)
w1$results
write.csv(w1$results, file = "1st_order_alleles")
##pairwise interaction (order 2)
hf2<-hf(data = color_allele, w.order = 2, B = 100)
w2<-wtest(data = color_allele, y = color_pheno, w.order = 2, input.pval = 1, input.poolsize = 10, output.pval = 1, hf1=hf1, hf2=hf2)
w.diagnosis(data = color_allele, w.order = 2, n.rep = 1000,hf1 = hf1, hf2=hf2)
w.qqplot(data = color_allele, y = color_pheno, w.order = 2, input.poolsize = 200,hf1 = hf1,hf2=hf2)
w2$results
write.csv(w2$results, file = "2nd_order_alleles")

#####all genotypes in contigs###
color_allcontigs<-read.csv("wtest_fullcontig.csv",header = TRUE,sep = ";")
hf1 <- hf(data = color_allcontigs, w.order = 1, B = 400)
w1<-wtest(data = color_allcontigs, y = color_pheno, w.order = 1,output.pval = 1, hf1 = hf1)
w1$results
write.csv(w1$results, file = "1st_order_genotypes")
w.diagnosis(data = color_allcontigs, w.order = 1, n.rep = 1000,hf1 = hf1)
w.qqplot(data = color_allcontigs, y = color_pheno,w.order = 1, input.poolsize = 200,hf1 = hf1,hf2=hf2)
w1$results

w1$results["p-value"]<- -log10(w1$results["p-value"])
plot_genotypes<-cbind(w1$results["marker"],w1$results["p-value"])
write.csv(plot_genotypes, file ="genotypes")

##2nd order interactions
hf2<-hf(data = color_allcontigs, w.order = 2, B = 400)
w2<-wtest(data = color_allcontigs, y = color_pheno, w.order = 2, input.pval = 1, 
          input.poolsize = 10, output.pval = 1, hf1=hf1, hf2=hf2)
w.diagnosis(data = color_allcontigs, w.order = 2, n.rep = 1000,hf1 = hf1, hf2=hf2)
w.qqplot(data = color_allcontigs, y = color_pheno, w.order = 2, input.poolsize = 200,hf1 = hf1,hf2=hf2)
w2$results
write.csv(w2$results, file = "2nd_order_contigs")

