

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


#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("plyranges")



colnames(sizes) <- c("seqnames", "end")
sizes$start <- 1

roundUp <- function(x,m) m*ceiling(x / m)
sizes$end <- roundUp(sizes$end, 20000)
bins <- sizes %>% as_granges() %>% tile_ranges(width = 20000L)

  print(trait)
  snp_params <- read.table(paste(wk_dir,"all.SNPs.filt.dentition.imp.bslmm.dentition.",trait,".mean.params.txt", sep = ""), header = T, row.names = NULL, stringsAsFactors = FALSE)

  colnames(snp_params) <- c("seqnames", "pos","alpha","beta","gamma")
  snp_params$id <- paste(snp_params$seqnames, snp_params$pos, sep="-")
  snp_params$start <- snp_params$pos
  snp_params$end <- snp_params$pos
  snp_params$seqnames <- 	paste("chr", snp_params$seqnames, sep="")
  snp_params1 <- snp_params[, c(1,7,8,6,3,4,5)]
  	
  write.table(snp_params1, paste(wk_dir,"all.SNPs.filt.dentition.imp.bslmm.dentition.", trait,".mean.params.bed", sep = ""), quote = FALSE, row.names = FALSE, sep = '\t')
  print(trait)
  print('top hits snps')
  print(head(snp_params1[order(snp_params1$gamma,decreasing = TRUE),]))
  print('  ')

  snps <- snp_params1 %>% as_granges()  %>% filter_by_overlaps(bins)
  sums <- bins %>%
    join_overlap_inner(snps) %>%
    disjoin_ranges(sum_alpha_20kb = sum(alpha), 
                   sum_beta_20kb = sum(beta), 
                   sum_gamma_20kb = sum(gamma))
  
  sums <-  as(sums, "data.frame")
  write.table(sums, paste(wk_dir,"all.SNPs.filt.dentition.imp.bslmm.dentition.", trait,".mean.params.20kb.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = '\t')
  print('top hits windows')
  print(head(sums[order(sums$sum_gamma_20kb,decreasing = TRUE),]))








### overlap with genes ###

colnames(genes) <- c("seqnames", "start","end","gene")
head(genes)
dim(genes)
genes <- genes %>% as_granges()

  snps <- read.table(paste(wk_dir,"all.SNPs.filt.dentition.imp.bslmm.dentition.",trait,".mean.params.txt", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  colnames(snps) <- c("seqnames", "pos","alpha","beta","gamma")
  snps$id <- paste(snps$seqnames, snps$pos, sep="-")
  snps$start <- snps$pos
  snps$end <- snps$pos
  snps$seqnames <- paste("chr", snps$seqnames, sep="")
  snps1 <- snps[, c(1,7,8,6,3,4,5)]

  sums <- read.table(paste(wk_dir,"all.SNPs.filt.dentition.imp.bslmm.dentition.", trait,".mean.params.20kb.txt", sep = ""), header = TRUE, stringsAsFactors = FALSE)  
  snps <- snps %>% as_granges()
  sums <- sums %>% as_granges()
  
  snps <- join_overlap_left(snps, genes)
  snps <- as(snps, "data.frame")
  snps <- snps[order(snps$gamma, decreasing = TRUE),]
  write.table(snps, paste(wk_dir,"all.SNPs.filt.dentition.imp.bslmm.dentition.",trait,".mean.params.anno.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = '\t')
  sums <- join_overlap_left(sums, genes)
  sums <- as(sums, "data.frame")
  sums <- sums[order(sums$sum_gamma_20kb, decreasing = TRUE),]
  write.table(sums, paste(wk_dir,"all.SNPs.filt.dentition.imp.bslmm.dentition.", trait,".mean.params.20kb.anno.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = '\t')
  
  #snp_thresh <- quantile(snps$gamma,.99)[[1]]
  #snps <- snps[which(snps$gamma >= snp_thresh),]
  #write.table(snps, paste(wk_dir,"all.SNPs.filt.dentition.imp.bslmm.dentition.",trait,".mean.params.99Q.anno.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = '\t')
  #sum_thresh <- quantile(sums$sum_gamma_20kb,.99)[[1]]
  #sums <- sums[which(sums$sum_gamma_20kb >= sum_thresh),]
  #write.table(sums, paste(wk_dir,"all.SNPs.filt.dentition.imp.bslmm.dentition.", trait,".mean.params.20kb.99Q.anno.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = '\t')
