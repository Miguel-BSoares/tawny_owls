
##
#


## 1- RevertSam ##

list="593
594"

mkdir -p unmapped_v6_bam
mkdir -p tmp

for sample in $list
        do picard32 RevertSam \
          -I aligned_v6/${sample}.bam \
          -O unmapped_v6_bam/${sample}.bam \
          SANITIZE=TRUE \
          ATTRIBUTE_TO_CLEAR=XT \
          ATTRIBUTE_TO_CLEAR=XN \
          ATTRIBUTE_TO_CLEAR=AS \
          ATTRIBUTE_TO_CLEAR=OC \
          ATTRIBUTE_TO_CLEAR=OP \
          TMP_DIR=tmp \
done

## 2 - MergeBamAlignment ##

list="593
594"

mkdir -p merged_alignments_v5

for sample in $list
        do picard MergeBamAlignment \
        --ALIGNED_BAM aligned_v5/${sample}.bam \
        -UNMAPPED_BAM unmapped_v5_bam/${sample}.bam \
        -O merged_alignments_v5/${sample}.merged.bam \
        -VALIDATION_STRINGENCY SILENT \
        -R v5/<reference> \
        TMP_DIR=tmp \
done

## 3 - SortSam ##

mkdir -p sorted_v5 

for sample in $list

        do picard SortSam -I merged_alignments_v5/${sample}.bam -O sorted_v5/${sample}.bam -SO coordinate

done

## 4 - CollectAlignmentSummaryMetrics

mkdir -p statistics_v5

for sample in $list

		do picard CollectAlignmentSummaryMetrics -I sorted_v5/${sample}.bam -R v5_genome_copy2.fasta -O statistics_v5/

done

## 5 - AddReadGroup
	
picard AddOrReplaceReadGroups -I merged_alignments_lq/${sample}.merged.bam -O replaced_readgroup/${sample}.bam  --RGID 4 --RGLB lib1 --RGPL ILLUMINA --RGPU unit1 --RGSM 20


## 6 - MarkDuplicates

picard MarkDuplicates \
  --INPUT replaced_readgroup/${sample}.bam \
  --OUTPUT removed_duplicates_lq/${sample}.removed.bam \
  --METRICS_FILE dup_metrics.txt \
  --MAX_RECORDS_IN_RAM 99999  \
  --REMOVE_DUPLICATES true \
  --TMP_DIR tmp/

## 7 - SplitNCigarReads 

gatk SplitNCigarReads \
  -R v5_genome_copy2.fasta \
  -I removed_duplicates_lq/${sample}.removed.bam \
  -O splitcigar_lq/${sample}.split.bam \
  --tmp-dir tmp/

### 8 - HaplotypeCaller ## CREATE INDEXED FILES

gatk HaplotypeCaller \
  -R v5_genome_copy2.fasta \
  -I splitcigar_lq/${sample}.split.bam \
  -O gvfc_aftermarkduplicates/${sample}.g.vcf.gz \
  -ERC GVCF \
  --min-base-quality-score 30 \
  --tmp-dir tmp/

## 10 - Rename VCFs ## Because readgroups add to be added for haplotype caller ##

for sample in $list

        do picard RenameSampleInVcf -I ${sample}.g.vcf.gz -O renamed_vcf/${sample}.g.vcf.gz --NEW_SAMPLE_NAME ${sample}

done

## 11 - CombineVCF - make a list ##

gatk CombineGVCFs -R v5_genome_copy2.fasta \

	--variant gvfc_aftermarkduplicates/12A.g.vcf.gz \
	--variant gvfc_aftermarkduplicates/13A.g.vcf.gz \
	--variant gvfc_aftermarkduplicates/14A.g.vcf.gz \
	--variant gvfc_aftermarkduplicates/15A.g.vcf.gz \
	--variant gvfc_aftermarkduplicates/16A.g.vcf.gz \
	--variant gvfc_aftermarkduplicates/17A.g.vcf.gz \
	--variant gvfc_aftermarkduplicates/18A.g.vcf.gz \
	--variant gvfc_aftermarkduplicates/18BA.g.vcf.gz \
	--variant gvfc_aftermarkduplicates/19A.g.vcf.gz \
	--variant gvfc_aftermarkduplicates/1A.g.vcf.gz \
	--variant gvfc_aftermarkduplicates/20A.g.vcf.gz \
	--variant gvfc_aftermarkduplicates/20BA.g.vcf.gz \
	--variant gvfc_aftermarkduplicates/26A.g.vcf.gz \
	--variant gvfc_aftermarkduplicates/28A.g.vcf.gz \
	--variant gvfc_aftermarkduplicates/29A.g.vcf.gz \
	--variant gvfc_aftermarkduplicates/2A.g.vcf.gz \
	--variant gvfc_aftermarkduplicates/30A.g.vcf.gz \
	--variant gvfc_aftermarkduplicates/31A.g.vcf.gz \
	--variant gvfc_aftermarkduplicates/32A.g.vcf.gz \
	--variant gvfc_aftermarkduplicates/33A.g.vcf.gz \
	--variant gvfc_aftermarkduplicates/34A.g.vcf.gz \
	--variant gvfc_aftermarkduplicates/35A.g.vcf.gz \
	--variant gvfc_aftermarkduplicates/36A.g.vcf.gz \
	--variant gvfc_aftermarkduplicates/37A.g.vcf.gz \
	--variant gvfc_aftermarkduplicates/41A.g.vcf.gz \
	--variant gvfc_aftermarkduplicates/7A.g.vcf.gz \
	--variant gvfc_aftermarkduplicates/9A.g.vcf.gz \
	-O combined/rna_joint.g.vcf.gz --tmp-dir tmp/

## 11 - JointCall

gatk GenotypeGVCFs -R v5/v5_genome_copy2.fasta -V combined_v5/joint.g.vcf.gz -O genotype_gvcfs_final_v5/joint_v5.vcf.gz -stand-call-conf 30 --tmp-dir tmp/

## 12 - Subset vcf in SNPs and Indels to be filtered independently (default/suggested parameters)

gatk SelectVariants -V joint_v5.vcf.gz -O snps_v5.vcf.gz -select-type SNP --tmp-dir tmp/

## 13 - VariantFiltration 

gatk VariantFiltration \
                -R ../../v5_genome_copy2.fasta \
                -V rna_snps.vcf.gz \
                --window 35 \
                --cluster 3 \
                -filter "QD < 2.0" --filter-name "QD2" \
                -filter "QUAL < 30.0" --filter-name "QUAL30" \
                -filter "SOR > 2.0" --filter-name "SOR2" \
                -filter "FS > 60.0" --filter-name "FS60" \
                -filter "MQ < 60.0" --filter-name "MQ60" \
                -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
                -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
                -O rna_snps_filtered.vcf.gz

##Utilize bedtools to collect SNPs located in gff annotated file

intersectBed -a rna_snps_filtered.vcf.gz -b v5_curated_names.gff3 -header > rna_cds.vcf

##Filter just those who pass the previous filters

bcftools view -f 'PASS' rna_cds.vcf > rna_cds_unique_final.vcf

bcftools view --max-alleles 2 -i 'F_MISSING<0.3' -q 0.05:minor rna_cds_unique_final.vcf > rna_cds_unique_biallelic_g03_maf005_final.vcf


##		
