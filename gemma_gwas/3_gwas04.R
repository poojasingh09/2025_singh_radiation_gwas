## pooja.singh09@gmail.com
## this scripts annotates SNPs with genes using a gff (keep only gene annotations in gff)


library(plyranges)
library(dplyr)
wk_dir <- "./"

sizes <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2023_bcftools_pnye3/pnye_v3.fa.fai.2col", header=FALSE, stringsAsFactors = F)
genes <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2023_bcftools_pnye3/Bouillabase_Tilapia_NCBI_nyerereitrack.pundcross.gapsEstimated.genes2pnye2.bed", stringsAsFactors = FALSE)

args<-commandArgs(trailingOnly=T)
trait=args[1]


### overlap with genes ### top p_score 1%

colnames(genes) <- c("seqnames", "start","end","gene")
head(genes)
genes <- genes %>% as_granges()

  snps1 <- read.table(paste0("all.SNPs.filt.jaw.imp.lmm.jaw.", trait, ".assoc.txt"), header = TRUE, stringsAsFactors = FALSE)
  snps1$rs <- paste0(snps1$chr, "_",snps1$ps)
  snps1$chr <- paste0("chr", snps1$chr)
  snps <- snps1[snps1$p_score < 0.01,]
  snps <- snps[,c(1,2,3,7,12,14)]
  colnames(snps) <- c("seqnames", "id", "pos","af","p_wald","p_score")
  snps$start <- snps$pos
  snps$end <- snps$pos

  snps <- snps %>% as_granges()
  
  snps <- join_overlap_left(snps, genes)
  snps <- as(snps, "data.frame")
  snps <- snps[order(snps$p_score, decreasing = FALSE),]
  write.table(snps, paste0("all.SNPs.filt.jaw.imp.lmm.jaw.", trait,".assoc.top1percent.anno.txt"), quote = FALSE, row.names = FALSE, sep = '\t')
