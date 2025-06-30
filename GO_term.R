#GOterm anaylsis sex DEGs owls
#DEGS from sex effect ~ sex + morph + sex:morph

#set working directory:
setwd("C:/Users/melan/OneDrive/Dokumente/Results/Tawny owl/htseq/hisat_htseq_counts_final/GOterms")


#load packages:
library(ape)
library(ggplot2)
library(ggfortify)
library(gplots)
library(annotate)
library(GOstats)
library(biomaRt)
library(AnnotationForge)
library(goEnrichment) # download failed
library(AnnotationDbi)
library(GSEABase)
library(stringr)
library(data.table)
library(tidyr)
library(dplyr)
library(reshape2)
library(cowplot)
library(RColorBrewer)
theme_set(theme_cowplot())




#functions:
getAttributeField <- function (x, field, attrsep = ";") {
  s = strsplit(x, split = attrsep, fixed = TRUE)
  sapply(s, function(atts) {
    a = strsplit(atts, split = "=", fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      rv = a[[m]][2]
    }
    else {
      rv = as.character(NA)
    }
    return(rv)
  })
}







#------------------load data-------------------------------------------------------------------------------------------
DEGs.sex <- read.csv2("C:/Users/melan/OneDrive/Dokumente/Results/Tawny owl/htseq/hisat_htseq_counts_final/htseq_counts_padj0.1_sex_lq_final.csv",header = T) %>% unique()
DEGs.sex$X # gene ids that are used for GOterm enrichment analysis

#load data and adjust counts
pheno_Data <- read.csv2("C:/Users/melan/OneDrive/Dokumente/Results/Tawny owl/sample_info.csv",row.names = 1)
pheno_Data_LQ <- subset(pheno_Data,library_type=="low_Q")


countData <- read.csv2("C:/Users/melan/OneDrive/Dokumente/Results/Tawny owl/htseq/hisat_htseq_counts_final/read_counts/htseq_counts_final.csv",row.names=1, header=T, check.names = FALSE)
countData_LQ <- countData[,rownames(pheno_Data_LQ)]

#ensures correct order
rownames(pheno_Data_LQ) == colnames(countData_LQ)#needs to be TRUE

#filter count data
hts_dds_lq <- DESeqDataSetFromMatrix(countData_LQ, pheno_Data_LQ, ~ Sex + Morph + Sex:Morph) #cov always first
hts_dds_lq <- hts_dds_lq[ rowSums(counts(hts_dds_lq)) > 1, ]#15846

rownames(counts(hts_dds_lq)) #gene universe to compare our DEGs to.




#------------------extract GOterms-------------------------------------------------------------------------------------------


#load annotation file and transfer GOterms from default version of the gff:
Ann.gff <- read.gff("C:/Users/melan/OneDrive/Dokumente/Results/Tawny owl/v5_test_gff3.gff3")
Default.gff <- read.gff("C:/Users/melan/OneDrive/Dokumente/Results/Tawny owl/grey_flye_1k_polished.evm.decorated.gff")

#create common identifyer that allows to reliably merge information:
mRNA.ann <- subset(Ann.gff,type=="mRNA")
mRNA.ann$ID <- paste(mRNA.ann$seqid,mRNA.ann$start,mRNA.ann$end,sep="_")
mRNA.ann$gene_id <- getAttributeField(mRNA.ann$attributes, "Parent")
mRNA.def <- subset(Default.gff,type=="mRNA")
mRNA.def$ID <- paste(mRNA.def$seqid,mRNA.def$start,mRNA.def$end,sep="_")
mRNA.def$GO <- getAttributeField(mRNA.def$attributes, "em_GOs")

mRNA.all <- merge(mRNA.ann[,c("ID","gene_id")],mRNA.def[,c("ID","GO")],by="ID",all.x=T)




#---------------create universe--------------------------------------------------------------------

# in the GOterm universe we only want genes that were sequenced, survived any filtering and have a GOterm attributed:
universe <- subset(mRNA.all,gene_id %in% rownames(counts(hts_dds_lq)))
universe <- na.omit(universe[,c("gene_id","GO")]) # check that this is gene_id and GO!!!
universe <- universe %>%
  mutate(go_id = strsplit(as.character(GO), ",")) %>%
  unnest(go_id)
universe$go_linkage_type <- "IEA"
universe <- as.data.frame(universe[,c("go_id","go_linkage_type","gene_id")])

#IMPORTANT: I just put IEA for now as I dont know what to put here,
#sticklebacks (biomart) have 82701 IEA and 63 ISS.

#create gene set collection
goFrame=GOFrame(universe,organism="Strix aluco")
goAllFrame=GOAllFrame(goFrame)
gsc_universe <- GeneSetCollection(goAllFrame, setType = GOCollection())


###### condtional!!!
#The GO ontology is set up as a directed acyclic graph, where a parent term is comprised of all its child terms. If you do a standard
#hypergeometric, you might e.g., find 'positive regulation of kinase activity' to be significant.
#If you then test 'positive regulation of catalytic activity', which is a parent term, then it might be significant as well, but only because of
#the terms coming from positive regulation of kinase activity.

#The conditional hypergeometric takes this into account, and only uses
#those terms that were not already significant when testing a higher
#order (parent) term.




#--------------------------------- GOterm analysis for sex DEGs -------------------------------------------------

params_sequenced <- GSEAGOHyperGParams(name="GO_over_MF_conditional_DEG_sex",
                                       geneSetCollection=gsc_universe,
                                       geneIds = as.vector(unique(DEGs.sex$X)),
                                       universeGeneIds = as.vector(universe[,3]),
                                       ontology = "MF",
                                       pvalueCutoff = 0.05,
                                       conditional = TRUE,
                                       testDirection = "over")
Over_sequenced_MF <- hyperGTest(params_sequenced)

Over_sequenced_MF_df <- as.data.frame(summary(Over_sequenced_MF)); head(Over_sequenced_MF_df)
htmlReport(Over_sequenced_MF, file="GO_over_MF_conditional_DEG_sex.html")
write.table(as.data.frame(summary(Over_sequenced_MF)),"GO_over_MF_conditional_DEG_sex.txt", sep="\t")


params_sequenced <- GSEAGOHyperGParams(name="GO_over_BP_conditional_DEG_sex",
                                       geneSetCollection=gsc_universe,
                                       geneIds = as.vector(unique(DEGs.sex$X)),
                                       universeGeneIds = as.vector(universe[,3]),
                                       ontology = "BP",
                                       pvalueCutoff = 0.05,
                                       conditional = TRUE,
                                       testDirection = "over")
Over_sequenced_BP <- hyperGTest(params_sequenced)

Over_sequenced_BP_df <- as.data.frame(summary(Over_sequenced_BP)); head(Over_sequenced_BP_df)
htmlReport(Over_sequenced_BP, file="GO_over_BP_conditional_DEG_sex.html")
write.table(as.data.frame(summary(Over_sequenced_BP)),"GO_over_BP_conditional_DEG_sex.txt", sep="\t")


params_sequenced <- GSEAGOHyperGParams(name="GO_over_CC_conditional_DEG_sex",
                                       geneSetCollection=gsc_universe,
                                       geneIds = as.vector(unique(DEGs.sex$X)),
                                       universeGeneIds = as.vector(universe[,3]),
                                       ontology = "CC",
                                       pvalueCutoff = 0.05,
                                       conditional = TRUE,
                                       testDirection = "over")
Over_sequenced_CC <- hyperGTest(params_sequenced)

Over_sequenced_CC_df <- as.data.frame(summary(Over_sequenced_CC)); head(Over_sequenced_CC_df)
htmlReport(Over_sequenced_CC, file="GO_over_CC_conditional_DEG_sex.html")
write.table(as.data.frame(summary(Over_sequenced_CC)),"GO_over_CC_conditional_DEG_sex.txt", sep="\t")









################################# FDR with goEnrichment:


params_fdr <- joinGOEnrichResults(goEnrichTest(gsc=gsc_universe,
                                               gene.ids = as.vector(unique(DEGs.sex$X)),
                                               univ.gene.ids = as.vector(universe[,3]),
                                               ontologies = "MF",
                                               pvalue.cutoff = 0.05,
                                               cond = TRUE,
                                               test.dir = "over"),
                                  p.adjust.method = "fdr")

head(params_fdr)
write.table(params_fdr,"FDR_goEnr_GO_over_MF_conditional_DEG_sex.txt", sep="\t")


params_fdr <- joinGOEnrichResults(goEnrichTest(gsc=gsc_universe,
                                               gene.ids = as.vector(unique(DEGs.sex$X)),
                                               univ.gene.ids = as.vector(universe[,3]),
                                               ontologies = "BP",
                                               pvalue.cutoff = 0.05,
                                               cond = TRUE,
                                               test.dir = "over"),
                                  p.adjust.method = "fdr")

head(params_fdr)
write.table(params_fdr,"FDR_goEnr_GO_over_BP_conditional_DEG_sex.txt", sep="\t")

params_fdr <- joinGOEnrichResults(goEnrichTest(gsc=gsc_universe,
                                               gene.ids = as.vector(unique(DEGs.sex$X)),
                                               univ.gene.ids = as.vector(universe[,3]),
                                               ontologies = "CC",
                                               pvalue.cutoff = 0.05,
                                               cond = TRUE,
                                               test.dir = "over"),
                                  p.adjust.method = "fdr")

head(params_fdr)
write.table(params_fdr,"FDR_goEnr_GO_over_CC_conditional_DEG_sex.txt", sep="\t")

