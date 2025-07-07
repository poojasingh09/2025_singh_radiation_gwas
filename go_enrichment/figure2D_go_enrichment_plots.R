

## pooja.singh09@gmail.com
## july 2024

library(ggplot2)
library(cowplot)
library(ggpubr)
library(ComplexHeatmap)
if (!require("devtools")) {
  install.packages("devtools", dependencies = TRUE)
  library(devtools)
}
install_github("raivokolde/pheatmap")
library(pheatmap)
library(RColorBrewer)

#make figure 2 for paper

setwd("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/go_enrichment/")
fb <- read.table("FB_topgo_fisher_weight_tests_BP_q0.05.txt", header=T, sep="\t")
fb$trait <- "FB"
yf <- read.table("YF_topgo_fisher_weight_tests_BP_q0.05.txt", header=T, sep="\t")
yf$trait <- "YF"
bb <- read.table("BB_topgo_fisher_weight_tests_BP_q0.05.txt", header=T, sep="\t")
bb$trait <- "BB"
mb <- read.table("MB_topgo_fisher_weight_tests_BP_q0.05.txt", header=T, sep="\t")
mb$trait <- "MB"
rc <- read.table("RC_topgo_fisher_weight_tests_BP_q0.05.txt", header=T, sep="\t")
rc$trait <- "RC"
vb <- read.table("VB_topgo_fisher_weight_tests_BP_q0.05.txt", header=T, sep="\t")
vb$trait <- "VB"

mlb <- read.table("MLB_topgo_fisher_weight_tests_BP_q0.05.txt", header=T, sep="\t")
mlb$trait <- "MLB"

#dlb <- read.table("DLB_topgo_fisher_weight_tests_BP_q0.05.txt", header=T, sep="\t")
#dlb$trait <- "DLB"

pc1 <- read.table("PC1_topgo_fisher_weight_tests_BP_q0.05.txt", header=T, sep="\t")
pc1$trait <- "PC1"

#pc2 <- read.table("PC2_topgo_fisher_weight_tests_BP_q0.05.txt", header=T, sep="\t")
#pc2$trait <- "PC2"

setwd("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/go_enrichment/")


#kd <- read.table("kd_topgo_fisher_weight_tests_BP_q0.05.txt", header=T, sep="\t")
#kd$trait <- "KD"

cs <- read.table("cs_topgo_fisher_weight_tests_BP_q0.05.txt", header=T, sep="\t")
cs$trait <- "CS"

ir <- read.table("ir_topgo_fisher_weight_tests_BP_q0.05.txt", header=T, sep="\t")
ir$trait <- "IR"

tm <- read.table("tm_topgo_fisher_weight_tests_BP_q0.05.txt", header=T, sep="\t")
tm$trait <- "TM"


## reformat data for heatmap
traits <- list(cs, ir, tm, pc1, vb,mlb, fb, yf, rc, mb, bb )
traits_new <- list()
counter=0
for (trait in traits){
  counter = counter + 1
  trait1 <- trait[,c(2,11)]
  trait1$topgoFisher.padjust <- -log10(trait1$topgoFisher.padjust)
  colnames(trait1) <- c("term", trait[1,12])
  traits_new[[counter]] <- trait1
  
}

input <- Reduce(function(x, y) merge(x, y, all=T, by=c("term")), traits_new, accumulate=F)

input1 <- input[,c(2:12)]
rownames(input1) <- input$term


svg("figure_2d_go_enrichment_heatmap_v2.svg", height=10, width=8)
pheatmap(as.matrix(input1), legend=T, cluster_rows=FALSE, 
         cluster_cols=FALSE, na_col="white", cellwidth = 12, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="YlGn")))(10))

dev.off()


## full dataset


setwd("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/go_enrichment/")
fb <- read.table("FB_topgo_fisher_weight_tests_BP_full.txt", header=T, sep="\t")
fb$trait <- "FB"
yf <- read.table("YF_topgo_fisher_weight_tests_BP_full.txt", header=T, sep="\t")
yf$trait <- "YF"
bb <- read.table("BB_topgo_fisher_weight_tests_BP_full.txt", header=T, sep="\t")
bb$trait <- "BB"
mb <- read.table("MB_topgo_fisher_weight_tests_BP_full.txt", header=T, sep="\t")
mb$trait <- "MB"
rc <- read.table("RC_topgo_fisher_weight_tests_BP_full.txt", header=T, sep="\t")
rc$trait <- "RC"
vb <- read.table("VB_topgo_fisher_weight_tests_BP_full.txt", header=T, sep="\t")
vb$trait <- "VB"

mlb <- read.table("MLB_topgo_fisher_weight_tests_BP_full.txt", header=T, sep="\t")
mlb$trait <- "MLB"

dlb <- read.table("DLB_topgo_fisher_weight_tests_BP_full.txt", header=T, sep="\t")
dlb$trait <- "DLB"

pc1 <- read.table("PC1_topgo_fisher_weight_tests_BP_full.txt", header=T, sep="\t")
pc1$trait <- "PC1"

pc2 <- read.table("PC2_topgo_fisher_weight_tests_BP_full.txt", header=T, sep="\t")
pc2$trait <- "PC2"

setwd("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/go_enrichment/")


kd <- read.table("kd_topgo_fisher_weight_tests_BP_full.txt", header=T, sep="\t")
kd$trait <- "KD"

cs <- read.table("cs_topgo_fisher_weight_tests_BP_full.txt", header=T, sep="\t")
cs$trait <- "CS"

ir <- read.table("ir_topgo_fisher_weight_tests_BP_full.txt", header=T, sep="\t")
ir$trait <- "IR"

tm <- read.table("tm_topgo_fisher_weight_tests_BP_full.txt", header=T, sep="\t")
tm$trait <- "TM"


## reformat data for heatmap
traits <- list(cs, ir, kd, tm, pc1,pc2, vb,mlb, dlb, fb, yf, rc, mb, bb )
traits_new <- list()
counter=0
for (trait in traits){
  counter = counter + 1
  trait1 <- trait[,c(2,11)]
  trait1$topgoFisher.padjust <- -log10(trait1$topgoFisher.padjust)
  colnames(trait1) <- c("term", trait[1,12])
  traits_new[[counter]] <- trait1
  
}

input <- Reduce(function(x, y) merge(x, y, all=T, by=c("term")), traits_new, accumulate=F)

input1 <- input[,c(2:15)]
rownames(input1) <- input$term


svg("figure_2d_go_enrichment_heatmap_v2_full.svg", height=14, width=8)
pheatmap(as.matrix(input1), legend=T, cluster_rows=FALSE, 
         cluster_cols=FALSE, na_col="white", cellwidth = 12, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="YlGn")))(10))

dev.off()





