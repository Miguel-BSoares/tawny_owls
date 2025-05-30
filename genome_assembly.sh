####All analyses were run in a HPC environment @CSC – IT Center for Science, Finland####
####Different assemblies were produce by switching -m in flye
####Comparisons were performed with Quast that also computes BUSCO.

####################################
##alignment of raw reads (sub-reads)
flye --pacbio-raw owl_data/pacbio.V5A.subreads.fastq.gz --out-dir v5_flye_2 --resume -m 1000 --genome-size 1g --asm-coverage 20 --threads 4

#polishing assembly
pbmm2 index flye_assemblies/flye_v6_assembly_2_1000.fasta pbmm2_v6_flye_2_1000/v6_flye_2_1000_assembly.mmi
pbmm2 align --sort --log-level DEBUG pbmm2_v6_flye_2_1000/v6_flye_2_1000_assembly.mmi owl_data/m64048_210523_061850.subreads.bam pbmm2_v6_flye_2_1000/v6_flye_2_1000.bam
gcpp pbmm2_v5_flye/v5_flye.bam -r flye_assemblies/v5_flye_assembly.fasta -o v5_flye_polished/v5_flye_polished.fasta --algorithm=arrow -j 20

#comparison and BUSCO
quast.py flye_assemblies/polished/*fasta --output-dir quast_flye_polished -k -b -e -f --eukaryote --large --est-ref-size 1 --threads 4

#plot BUSCO
library(tidyverse)
library(wesanderson)
setwd("E:/NTFS/TURKU - CURRENT JOB/denovo_assembly/assemblies/genome/busco")
aves_busco<-read.csv("E:/NTFS/TURKU - CURRENT JOB/denovo_assembly/assemblies/genome/busco/aves.csv", sep = ";" , header = TRUE)
p <- ggplot(aves_busco, aes(fill= factor(Type,levels=c("Complete and single copy","Complete and duplicated","Fragmented","Missing")), y=BUSCOS, x=ï..Genome)) +
  geom_bar(position="stack", stat="identity", color = "black") +
  labs(y="Total",x="Genomes",fill="")+
  scale_fill_manual(values = wes_palette("Zissou1")) +
  ggtitle("BUSCO results for Aves database")
p + theme_bw()+coord_flip()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text=element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14))

         
setwd("E:/NTFS/TURKU - CURRENT JOB/denovo_assembly/assemblies/genome/busco")
aves_busco<-read.csv("E:/NTFS/TURKU - CURRENT JOB/denovo_assembly/assemblies/genome/busco/owls.csv", sep = ";" , header = TRUE)
p <- ggplot(aves_busco, aes(fill= factor(Type,levels=c("Complete and single copy","Complete and duplicated","Fragmented","Missing")), y=BUSCOS, x=ï..Genome)) +
  geom_bar(position="stack", stat="identity", color = "black") +
  labs(y="Total",x="Genomes",fill="")+
  scale_fill_manual(values = wes_palette("Zissou1")) +
  ggtitle("BUSCO results for Eukaryote database")
p + theme_bw()+coord_flip()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text=element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14))

