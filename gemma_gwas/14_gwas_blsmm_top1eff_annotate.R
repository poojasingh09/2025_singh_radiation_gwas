

# This script is used after trait_name_gemma_mean_parameters.py
# This script is used to sum GEMMA parameters across 20kb windows
#------------------------------------------------------------------------------------------
# set working directory / trait names / .gff file directory / scaffold sizes file directory
#------------------------------------------------------------------------------------------


library(plyranges)
library(dplyr)
wk_dir <- "./"

# Reads input parameters
args<-commandArgs(trailingOnly=T)
trait=args[1]


sizes <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2023_bcftools_pnye3/pnye_v3.fa.fai.2col", header=FALSE, stringsAsFactors = F)
#genes <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2023_bcftools_pnye3/Bouillabase_Tilapia_NCBI_nyerereitrack.pundcross.gapsEstimated.genes2pnye2.bed", stringsAsFactors = FALSE)
genes <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/pnye_v3_assembly/punda.v3.annotation_v2.gene.bed.txt", stringsAsFactors = FALSE)

### overlap with genes ### top EFF 1%

colnames(genes) <- c("seqnames", "start","end","gene")
head(genes)
genes <- genes %>% as_granges()


  snps <- read.table(paste(wk_dir,"all.SNPs.filt.dentition.imp.bslmm.dentition.",trait,".top0.1eff.dsv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  colnames(snps) <- c("seqnames", "pos","alpha","beta","gamma", "eff")
  snps$id <- paste(snps$seqnames, snps$pos, sep="-")
  snps$start <- snps$pos
  snps$end <- snps$pos
  snps$seqnames <- paste("chr", snps$seqnames, sep="")
  snps1 <- snps[, c(1,7,8,6,3,4,5)]

  snps <- snps %>% as_granges()
  
  snps <- join_overlap_left(snps, genes)
  snps <- as(snps, "data.frame")
  snps <- snps[order(snps$gamma, decreasing = TRUE),]
  write.table(snps, paste(wk_dir,"all.SNPs.filt.dentition.imp.bslmm.dentition.",trait,".top1eff.anno.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = '\t')



## top PIP 1 %

  snps <- read.table(paste(wk_dir,"all.SNPs.filt.dentition.imp.bslmm.dentition.",trait,".pip01.dsv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  colnames(snps) <- c("seqnames", "pos","alpha","beta","gamma", "eff")
  snps$id <- paste(snps$seqnames, snps$pos, sep="-")
  snps$start <- snps$pos
  snps$end <- snps$pos
  snps$seqnames <- paste("chr", snps$seqnames, sep="")
  snps1 <- snps[, c(1,7,8,6,3,4,5)]

  snps <- snps %>% as_granges()

  snps <- join_overlap_left(snps, genes)
  snps <- as(snps, "data.frame")
  snps <- snps[order(snps$gamma, decreasing = TRUE),]
  write.table(snps, paste(wk_dir,"all.SNPs.filt.dentition.imp.bslmm.dentition.",trait,".pip01.anno.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = '\t')



## top PIP 5 %

  snps <- read.table(paste(wk_dir,"all.SNPs.filt.dentition.imp.bslmm.dentition.",trait,".pip05.dsv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  colnames(snps) <- c("seqnames", "pos","alpha","beta","gamma", "eff")
  snps$id <- paste(snps$seqnames, snps$pos, sep="-")
  snps$start <- snps$pos
  snps$end <- snps$pos
  snps$seqnames <- paste("chr", snps$seqnames, sep="")
  snps1 <- snps[, c(1,7,8,6,3,4,5)]

  snps <- snps %>% as_granges()

  snps <- join_overlap_left(snps, genes)
  snps <- as(snps, "data.frame")
  snps <- snps[order(snps$gamma, decreasing = TRUE),]
  write.table(snps, paste(wk_dir,"all.SNPs.filt.dentition.imp.bslmm.dentition.",trait,".pip05.anno.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = '\t')

