## making manhttan plots with GWAS candidates highlighted
## pooja.singh09@gmail.com


#install.packages('qqman')
library(qqman)
library(scales)
library(data.table)
library(dplyr)
library(GenomicRanges)
library(plyranges)
library(IRanges)


setwd('/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/fst/')
file_dir <- '/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_dentition/fst'
## read in gene annotations

genes <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/pnye_v3_assembly/punda.v3.annotation_v2.gene.bed.txt", stringsAsFactors = FALSE)
colnames(genes) <- c("seqnames", "start","end","anno")
#genes <- genes %>% as_granges()


## get outliers to highlight

setwd('/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/snp_effect_overlaps/')
mb_outliers <- fread("mb_top_99.9_and_pip01.txt", data.table=F)
mb_outliers$SNP <- paste0(mb_outliers$seqnames, "_",mb_outliers$start)
mb_outliers$SNP <- gsub("chr", "", mb_outliers$SNP)

fb_outliers <- fread("fb_top_99.9_and_pip01.txt",header=T,sep="\t", data.table=F)
fb_outliers$SNP <- paste0(fb_outliers$seqnames, "_",fb_outliers$start)
fb_outliers$SNP <- gsub("chr", "", fb_outliers$SNP)


yf_outliers <- fread("yf_top_99.9_and_pip01.txt",header=T,sep="\t", data.table=F)
yf_outliers$SNP <- paste0(yf_outliers$seqnames, "_",yf_outliers$start)
yf_outliers$SNP <- gsub("chr", "", yf_outliers$SNP)


bb_outliers  <- fread("bb_top_99.9_and_pip01.txt",header=T,sep="\t", data.table=F)
bb_outliers$SNP <- paste0(bb_outliers$seqnames, "_",bb_outliers$start)
bb_outliers$SNP <- gsub("chr", "", bb_outliers$SNP)

rc_outliers  <- fread("rc_top_99.9_and_pip01.txt",header=T,sep="\t", data.table=F)
rc_outliers$SNP <- paste0(rc_outliers$seqnames, "_",rc_outliers$start)
rc_outliers$SNP <- gsub("chr", "", rc_outliers$SNP)


dlb_outliers  <- fread("dlb_top_99.9_and_pip01.txt",header=T,sep="\t", data.table=F)
dlb_outliers$SNP <- paste0(dlb_outliers$seqnames, "_",dlb_outliers$start)
dlb_outliers$SNP <- gsub("chr", "", dlb_outliers$SNP)


mlb_outliers  <- fread("mlb_top_99.9_and_pip01.txt",header=T,sep="\t", data.table=F)
mlb_outliers$SNP <- paste0(mlb_outliers$seqnames, "_",mlb_outliers$start)
mlb_outliers$SNP <- gsub("chr", "", mlb_outliers$SNP)

vb_outliers  <- fread("vb_top_99.9_and_pip01.txt",header=T,sep="\t", data.table=F)
vb_outliers$SNP <- paste0(vb_outliers$seqnames, "_",vb_outliers$start)
vb_outliers$SNP <- gsub("chr", "", vb_outliers$SNP)


pc1_outliers  <- fread("pc1_top_99.9_and_pip01.txt",header=T,sep="\t", data.table=F)
pc1_outliers$SNP <- paste0(pc1_outliers$seqnames, "_",pc1_outliers$start)
pc1_outliers$SNP <- gsub("chr", "", pc1_outliers$SNP)

pc2_outliers  <- fread("pc2_top_99.9_and_pip01.txt",header=T,sep="\t", data.table=F)
pc2_outliers$SNP <- paste0(pc2_outliers$seqnames, "_",pc2_outliers$start)
pc2_outliers$SNP <- gsub("chr", "", pc2_outliers$SNP)

cs_outliers  <- fread("cs_top_99.9_and_pip01.txt",header=T,sep="\t", data.table=F)
cs_outliers$SNP <- paste0(cs_outliers$seqnames, "_",cs_outliers$start)
cs_outliers$SNP <- gsub("chr", "", cs_outliers$SNP)

ir_outliers  <- fread("ir_top_99.9_and_pip01.txt",header=T,sep="\t", data.table=F)
ir_outliers$SNP <- paste0(ir_outliers$seqnames, "_",ir_outliers$start)
ir_outliers$SNP <- gsub("chr", "", ir_outliers$SNP)

kd_outliers  <- fread("kd_top_99.9_and_pip01.txt",header=T,sep="\t", data.table=F)
kd_outliers$SNP <- paste0(kd_outliers$seqnames, "_",kd_outliers$start)
kd_outliers$SNP <- gsub("chr", "", kd_outliers$SNP)

tm_outliers  <- fread("tm_top_99.9_and_pip01.txt",header=T,sep="\t", data.table=F)
tm_outliers$SNP <- paste0(tm_outliers$seqnames, "_",tm_outliers$start)
tm_outliers$SNP <- gsub("chr", "", tm_outliers$SNP)


## modified function for manhattan plots to highlight outliers from gwas

manhattan1<-function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c(alpha("gray60",0.5),
                                                                               alpha("gray80", 0.5)), chrlabs = NULL, suggestiveline = -log10(1e-05),
                      genomewideline = -log10(5e-08), highlight1 = NULL, highlight2 = NULL, highlight3 = NULL, highlight4 = NULL, logp = TRUE,
                      highlight5 = NULL, highlight6 = NULL, highlight7 = NULL, highlight8 = NULL, highlight9 = NULL, highlight10 = NULL, highlight11 = NULL, highlight12 = NULL, highlight13 = NULL, highlight14 = NULL,
                      ...)
{
  CHR = BP = P = index = NULL
  if (!(chr %in% names(x)))
    stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x)))
    stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x)))
    stop(paste("Column", p, "not found!"))
  if (!(snp %in% names(x)))
    warning(paste("No SNP column found. OK unless you're trying to highlight."))
  if (!is.numeric(x[[chr]]))
    stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]]))
    stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]]))
    stop(paste(p, "column should be numeric."))
  d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
  if (!is.null(x[[snp]]))
    d = transform(d, SNP = x[[snp]])
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d <- d[order(d$CHR, d$BP), ]
  if (logp) {
    d$logp <- -log10(d$P)
  }
  else {
    d$logp <- d$P
  }
  d$pos = NA
  d$index = NA
  ind = 0
  for (i in unique(d$CHR)) {
    ind = ind + 1
    d[d$CHR == i, ]$index = ind
  }
  nchr = length(unique(d$CHR))
  if (nchr == 1) {
    options(scipen = 999)
    d$pos = d$BP/1e+06
    ticks = floor(length(d$pos))/2 + 1
    xlabel = paste("Chromosome", unique(d$CHR), "position(Mb)")
    labs = ticks
  }
  else {
    lastbase = 0
    ticks = NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos = d[d$index == i, ]$BP
      }
      else {
        lastbase = lastbase + tail(subset(d, index ==
                                            i - 1)$BP, 1)
        d[d$index == i, ]$pos = d[d$index == i, ]$BP +
          lastbase
      }
      ticks = c(ticks, (min(d[d$CHR == i, ]$pos) + max(d[d$CHR ==
                                                           i, ]$pos))/2 + 1)
    }
    xlabel = "Chromosome"
    labs <- unique(d$CHR)
  }
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i",
                   las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0,
                                                                     ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log10))
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in%
                                            names(dotargs)]))
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      }
      else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    }
    else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  if (nchr == 1) {
    axis(1, ...)
  }
  else {
    axis(1, at = ticks, labels = labs, ...)
  }
  col = rep(col, max(d$CHR))
  if (nchr == 1) {
    with(d, points(pos, logp, pch = 20, col = col[1], ...))
  }
  else {
    icol = 1
    for (i in unique(d$index)) {
      with(d[d$index == unique(d$index)[i], ], points(pos,
                                                      logp, col = col[icol], pch = 20, ...))
      icol = icol + 1
    }
  }
  if (suggestiveline)
    abline(h = suggestiveline, col = "blue")
  if (genomewideline)
    abline(h = genomewideline, col = "red")
  if (!is.null(highlight1)) {
    if (any(!(highlight1 %in% d$SNP)))
      warning("You're trying to highlight1 SNPs that don't exist in your results.")
    d.highlight1 = d[which(d$SNP %in% highlight1), ]
    with(d.highlight1, points(pos, logp, col = alpha("#E64B35", 1), pch = 20,
                              ...))
  }
  if (!is.null(highlight2)) {
    if (any(!(highlight2 %in% d$SNP)))
      warning("You're trying to highlight2 SNPs that don't exist in your results.")
    d.highlight2 = d[which(d$SNP %in% highlight2), ]
    with(d.highlight2, points(pos, logp, col = alpha("#E64B35", 1.0), pch = 20,
                              ...))
  }
  if (!is.null(highlight3)) {
    if (any(!(highlight3 %in% d$SNP)))
      warning("You're trying to highlight3 SNPs that don't exist in your results.")
    d.highlight3 = d[which(d$SNP %in% highlight3), ]
    with(d.highlight3, points(pos, logp, col = alpha("#E64B35", 1.0), pch = 20,
                              ...))
  }
  
  if (!is.null(highlight4)) {
    if (any(!(highlight4 %in% d$SNP)))
      warning("You're trying to highlight4 SNPs that don't exist in your results.")
    d.highlight4 = d[which(d$SNP %in% highlight4), ]
    with(d.highlight4, points(pos, logp, col = alpha("#E64B35", 1.0), pch = 20,
                              ...))
  }
  
  if (!is.null(highlight5)) {
    if (any(!(highlight5 %in% d$SNP)))
      warning("You're trying to highlight5 SNPs that don't exist in your results.")
    d.highlight5 = d[which(d$SNP %in% highlight5), ]
    with(d.highlight5, points(pos, logp, col = alpha("#E64B35", 1.0), pch = 20,
                              ...))
  }
  
  if (!is.null(highlight6)) {
    if (any(!(highlight6 %in% d$SNP)))
      warning("You're trying to highlight6 SNPs that don't exist in your results.")
    d.highlight6 = d[which(d$SNP %in% highlight6), ]
    with(d.highlight6, points(pos, logp, col = alpha("#009E73", 1.0), pch = 20,  
                              ...))
  }
  
  if (!is.null(highlight7)) {
    if (any(!(highlight7 %in% d$SNP)))
      warning("You're trying to highlight7 SNPs that don't exist in your results.")
    d.highlight7 = d[which(d$SNP %in% highlight7), ]
    with(d.highlight7, points(pos, logp, col = alpha("#009E73", 1.0), pch = 20,  
                              ...))
  }
  
  if (!is.null(highlight8)) {
    if (any(!(highlight8 %in% d$SNP)))
      warning("You're trying to highlight8 SNPs that don't exist in your results.")
    d.highlight8 = d[which(d$SNP %in% highlight8), ]
    with(d.highlight8, points(pos, logp, col = alpha("#009E73", 1.0), pch = 20,  
                              ...))
  }
  
  if (!is.null(highlight9)) {
    if (any(!(highlight9 %in% d$SNP)))
      warning("You're trying to highlight9 SNPs that don't exist in your results.")
    d.highlight9 = d[which(d$SNP %in% highlight9), ]
    with(d.highlight9, points(pos, logp, col = alpha("#4DBBD5", 1.0), pch = 20,  
                              ...))
  }
  
  if (!is.null(highlight10)) {
    if (any(!(highlight10 %in% d$SNP)))
      warning("You're trying to highlight10 SNPs that don't exist in your results.")
    d.highlight10 = d[which(d$SNP %in% highlight10), ]
    with(d.highlight10, points(pos, logp, col = alpha("#4DBBD5", 1), pch = 20,  
                               ...))
  }
  
  if (!is.null(highlight11)) {
    if (any(!(highlight11 %in% d$SNP)))
      warning("You're trying to highlight10 SNPs that don't exist in your results.")
    d.highlight11 = d[which(d$SNP %in% highlight11), ]
    with(d.highlight11, points(pos, logp, col = alpha("#4DBBD5", 1), pch = 20,  
                               ...))
  }
  
  if (!is.null(highlight12)) {
    if (any(!(highlight12 %in% d$SNP)))
      warning("You're trying to highlight10 SNPs that don't exist in your results.")
    d.highlight12 = d[which(d$SNP %in% highlight12), ]
    with(d.highlight12, points(pos, logp, col = alpha("#4DBBD5", 1), pch = 20,  
                               ...))
  }
  if (!is.null(highlight13)) {
    if (any(!(highlight13 %in% d$SNP)))
      warning("You're trying to highlight10 SNPs that don't exist in your results.")
    d.highlight13 = d[which(d$SNP %in% highlight13), ]
    with(d.highlight13, points(pos, logp, col = alpha("#4DBBD5", 1), pch = 20,  
                               ...))
  }
  if (!is.null(highlight14)) {
    if (any(!(highlight14 %in% d$SNP)))
      warning("You're trying to highlight10 SNPs that don't exist in your results.")
    d.highlight14 = d[which(d$SNP %in% highlight14), ]
    with(d.highlight14, points(pos, logp, col = alpha("#4DBBD5", 1), pch = 20,  
                               ...))
  }
  
  
}

### make manhattan plots with gwas outliers
setwd('/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_dentition/fst')

counter <- 0
all <- list()
files <- list.files(pattern = "\\.fst$")
for (file in files) {
  counter = counter +1
  name <- paste(strsplit(file, "\\_20kb")[[1]][1])
  print(name)
  fst <- read.table(file, header=T)
  head(fst)
  fstsubset<-fst[complete.cases(fst),]
  fstsubset$WEIGHTED_FST <- replace(fstsubset$WEIGHTED_FST, fstsubset$WEIGHTED_FST < 0, 0)
  #cutoff <- quantile(fstsubset$WEIGHTED_FST, prob=0.95)
  fstsubset_without_chr10 <- fstsubset[fstsubset$CHROM != "chr10",]
  cutoff <- quantile(fstsubset_without_chr10$WEIGHTED_FST, prob=0.95)

  
  SNP<-c(1:(nrow(fstsubset)))
  mydf<-data.frame(SNP,fstsubset)
  mydf$CHROM <- gsub('chr', '', mydf$CHROM)
  mydf$CHROM <- as.numeric(mydf$CHROM)
  mydf$SNP <- paste0(mydf$CHROM, "_", mydf$BIN_START)
  mydf$seqnames <- paste0("chr", mydf$CHROM)
  colnames(mydf) <- c("SNP","CHROM","start","end","N_VARIANTS","WEIGHTED_FST", "MEAN_FST","seqnames")
  mydf <- mydf[, c(8,3,4,1,2,5,6,7)]
  
  mydf <- data.table(mydf)
  setkeyv(mydf, names(mydf))
  
  yf_outliers <- data.table(yf_outliers)
  setkeyv(yf_outliers, names(yf_outliers))
  ans <- foverlaps(mydf, yf_outliers, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames","start","end"))
  yf_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  yf_outliers1$trait <- "YF"
  
  
  fb_outliers <- data.table(fb_outliers)
  setkeyv(fb_outliers, names(fb_outliers))
  ans <- foverlaps(mydf, fb_outliers, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames","start","end"))
  fb_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  fb_outliers1$trait <- "FB"
  
  bb_outliers <- data.table(bb_outliers)
  setkeyv(bb_outliers, names(bb_outliers))
  ans <- foverlaps(mydf, bb_outliers, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames","start","end"))
  bb_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  bb_outliers1$trait <- "BB"
  
  mb_outliers <- data.table(mb_outliers)
  setkeyv(mb_outliers, names(mb_outliers))
  ans <- foverlaps(mydf, mb_outliers, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames","start","end"))
  mb_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  mb_outliers1$trait <- "MB"
  
  rc_outliers <- data.table(rc_outliers)
  setkeyv(rc_outliers, names(rc_outliers))
  ans <- foverlaps(mydf, rc_outliers, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames","start","end"))
  rc_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  rc_outliers1$trait <- "RC"
  
  mlb_outliers <- data.table(mlb_outliers)
  setkeyv(mlb_outliers, names(mlb_outliers))
  ans <- foverlaps(mydf, mlb_outliers, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames","start","end"))
  mlb_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  mlb_outliers1$trait <- "MLB"
  
  dlb_outliers <- data.table(dlb_outliers)
  setkeyv(dlb_outliers, names(dlb_outliers))
  ans <- foverlaps(mydf, dlb_outliers, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames","start","end"))
  dlb_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  dlb_outliers1$trait <- "DLB"
  
  vb_outliers <- data.table(vb_outliers)
  setkeyv(vb_outliers, names(vb_outliers))
  ans <- foverlaps(mydf, vb_outliers, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames","start","end"))
  vb_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  vb_outliers1$trait <- "VB"
  
  pc1_outliers <- data.table(pc1_outliers)
  setkeyv(pc1_outliers, names(pc1_outliers))
  ans <- foverlaps(mydf, pc1_outliers, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames","start","end"))
  pc1_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  pc1_outliers1$trait <- "PC1"
  
  pc2_outliers <- data.table(pc2_outliers)
  setkeyv(pc2_outliers, names(pc2_outliers))
  ans <- foverlaps(mydf, pc2_outliers, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames","start","end"))
  pc2_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  pc2_outliers1$trait <- "PC2"
  
  cs_outliers <- data.table(cs_outliers)
  setkeyv(cs_outliers, names(cs_outliers))
  ans <- foverlaps(mydf, cs_outliers, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames","start","end"))
  cs_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  cs_outliers1$trait <- "CS"
  
  ir_outliers <- data.table(ir_outliers)
  setkeyv(ir_outliers, names(ir_outliers))
  ans <- foverlaps(mydf, ir_outliers, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames","start","end"))
  ir_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  ir_outliers1$trait <- "IR"
  
  kd_outliers <- data.table(kd_outliers)
  setkeyv(kd_outliers, names(kd_outliers))
  ans <- foverlaps(mydf, kd_outliers, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames","start","end"))
  kd_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  kd_outliers1$trait <- "KD"
  
  tm_outliers <- data.table(tm_outliers)
  setkeyv(tm_outliers, names(tm_outliers))
  ans <- foverlaps(mydf, tm_outliers, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames","start","end"))
  tm_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  tm_outliers1$trait <- "TM"
  
  
  #mydf1 <- sample_n(mydf, 100000)
  mydf1 <- unique(rbind(yf_outliers1, fb_outliers1, bb_outliers1,mb_outliers1,rc_outliers1,vb_outliers1,mlb_outliers1,dlb_outliers1,pc1_outliers1,pc2_outliers1,cs_outliers1,ir_outliers1,kd_outliers1,tm_outliers1))
  print(dim(mydf1))
  all[[file]] <- mydf1

  
  #pdf(paste0(file, ".pip01_GWAS_SNPoutliers_0.95_FST_20kb_outliers.pdf"), h=4, w=10)
  #svg(paste0(file, "mahanttan_gwasoutliers_0.999_gwas_0.95fst_topgenes.svg"), h=4, w=12)
  png(paste0("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/fst/",name, ".99.9_pip01_GWAS_SNPoutliers_0.95_FST_20kb_outliers.png"), width = 2800, height = 1000, res = 300)
  
  manhattan1(mydf,chr="CHROM",bp="start",p="WEIGHTED_FST",snp="SNP",logp=FALSE,ylab="Weir and Cockerham Fst 20Kb windows",
             main=paste0(name, " PIP > 0.01 & top 5% FST"), highlight1 = mb_outliers1$i.SNP,
             highlight2=fb_outliers1$i.SNP, highlight3=yf_outliers1$i.SNP, highlight4=bb_outliers1$i.SNP, highlight5=rc_outliers1$i.SNP, highlight6=dlb_outliers1$i.SNP, 
             highlight7=mlb_outliers1$i.SNP, highlight8=vb_outliers1$i.SNP, highlight9=pc1_outliers1$i.SNP, highlight10=pc2_outliers1$i.SNP,
             highlight11=cs_outliers1$i.SNP,highlight12=ir_outliers1$i.SNP,highlight13=kd_outliers1$i.SNP,highlight14=tm_outliers1$i.SNP, cex=1)
  dev.off()

}

all_anno1 <- rbindlist(all, fill=TRUE, idcol = "comparison")
setwd('/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/fst/')
write.table(all_anno1, "allcomparisons_top_gwas_pip01_SNP_0.95fst_tophits_SNPoverlaponly.txt", quote=F, row.names=F, sep="\t")
write.table(data.frame(table(all_anno1$comparison)), "allcomparisons_top_gwas_pip01_SNP_0.95fst_tophits_SNPoverlaponly_outliersummary.txt", quote=F, row.names=F, sep="\t")




###################################
### make manhattan plots with gwas outlier GENES 
setwd('/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_dentition/fst')

counter <- 0
all <- list()
files <- list.files(pattern = "\\.fst$")


for (file in files) {
  counter = counter +1
  fst <- read.table(file, header=T)
  head(fst)
  fstsubset<-fst[complete.cases(fst),]
  fstsubset$WEIGHTED_FST <- replace(fstsubset$WEIGHTED_FST, fstsubset$WEIGHTED_FST < 0,0)
  fstsubset_without_chr10 <- fstsubset[fstsubset$CHROM != "chr10",]
  cutoff <- quantile(fstsubset_without_chr10$WEIGHTED_FST, prob=0.95)
  #cutoff <- quantile(fstsubset$WEIGHTED_FST, prob=0.95)
  SNP<-c(1:(nrow(fstsubset)))
  mydf<-data.frame(SNP,fstsubset)
  mydf$CHROM <- gsub('chr', '', mydf$CHROM)
  mydf$CHROM <- as.numeric(mydf$CHROM)
  mydf$SNP <- paste0(mydf$CHROM, "_", mydf$BIN_START)
  mydf$seqnames <- paste0("chr", mydf$CHROM)
  colnames(mydf) <- c("SNP","CHROM","start","end","N_VARIANTS","WEIGHTED_FST", "MEAN_FST","seqnames")
  mydf <- mydf[, c(8,3,4,1,2,5,6,7)]
  
  mydf <- data.table(mydf)
  setkeyv(mydf, names(mydf))
  
  yf_outliers <- data.table(yf_outliers)
  yf_genes <- merge(genes, yf_outliers, by.y="gene", by.x="anno")
  yf_genes1 <- data.table(yf_genes[, c(2,3,4,15)]) # chr, pos, 
  yf_intergenic <- yf_outliers[!(complete.cases(yf_outliers$gene)),]
  yf_intergenic1 <- data.table(yf_intergenic[, c(1,2,3,12)])
  yf_all <- rbind(yf_genes1, yf_intergenic1, use.names=FALSE)
  setkeyv(yf_all, names(yf_all))
  ans <- foverlaps(mydf, yf_all, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames.x","start.x","end.x"))
  yf_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  yf_outliers1$trait <- "YF"
  
  
  fb_outliers <- data.table(fb_outliers)
  fb_genes <- merge(genes, fb_outliers, by.y="gene", by.x="anno")
  fb_genes1 <- data.table(fb_genes[, c(2,3,4,15)])
  fb_intergenic <- fb_outliers[!(complete.cases(fb_outliers$gene)),]
  fb_intergenic1 <- data.table(fb_intergenic[, c(1,2,3,12)])
  fb_all <- rbind(fb_genes1, fb_intergenic1, use.names=FALSE)
  setkeyv(fb_all, names(fb_all))
  ans <- foverlaps(mydf, fb_all, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames.x","start.x","end.x"))
  fb_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  fb_outliers1$trait <- "fb"
  
  mb_outliers <- data.table(mb_outliers)
  mb_genes <- merge(genes, mb_outliers, by.y="gene", by.x="anno")
  mb_genes1 <- data.table(mb_genes[, c(2,3,4,15)])
  mb_intergenic <- mb_outliers[!(complete.cases(mb_outliers$gene)),]
  mb_intergenic1 <- data.table(mb_intergenic[, c(1,2,3,12)])
  mb_all <- rbind(mb_genes1, mb_intergenic1, use.names=FALSE)
  setkeyv(mb_all, names(mb_all))
  ans <- foverlaps(mydf, mb_all, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames.x","start.x","end.x"))
  mb_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  mb_outliers1$trait <- "mb"
  
  rc_outliers <- data.table(rc_outliers)
  rc_genes <- merge(genes, rc_outliers, by.y="gene", by.x="anno")
  rc_genes1 <- data.table(rc_genes[, c(2,3,4,15)])
  rc_intergenic <- rc_outliers[!(complete.cases(rc_outliers$gene)),]
  rc_intergenic1 <- data.table(rc_intergenic[, c(1,2,3,12)])
  rc_all <- rbind(rc_genes1, rc_intergenic1, use.names=FALSE)
  setkeyv(rc_all, names(rc_all))
  ans <- foverlaps(mydf, rc_all, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames.x","start.x","end.x"))
  rc_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  rc_outliers1$trait <- "rc"
  
  mlb_outliers <- data.table(mlb_outliers)
  mlb_genes <- merge(genes, mlb_outliers, by.y="gene", by.x="anno")
  mlb_genes1 <- data.table(mlb_genes[, c(2,3,4,15)])
  mlb_intergenic <- mlb_outliers[!(complete.cases(mlb_outliers$gene)),]
  mlb_intergenic1 <- data.table(mlb_intergenic[, c(1,2,3,12)])
  mlb_all <- rbind(mlb_genes1, mlb_intergenic1, use.names=FALSE)
  setkeyv(mlb_all, names(mlb_all))
  ans <- foverlaps(mydf, mlb_all, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames.x","start.x","end.x"))
  mlb_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  mlb_outliers1$trait <- "mlb"
  
  dlb_outliers <- data.table(dlb_outliers)
  dlb_genes <- merge(genes, dlb_outliers, by.y="gene", by.x="anno")
  dlb_genes1 <- data.table(dlb_genes[, c(2,3,4,15)])
  dlb_intergenic <- dlb_outliers[!(complete.cases(dlb_outliers$gene)),]
  dlb_intergenic1 <- data.table(dlb_intergenic[, c(1,2,3,12)])
  dlb_all <- rbind(dlb_genes1, dlb_intergenic1, use.names=FALSE)
  setkeyv(dlb_all, names(dlb_all))
  ans <- foverlaps(mydf, dlb_all, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames.x","start.x","end.x"))
  dlb_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  dlb_outliers1$trait <- "dlb"
  
  vb_outliers <- data.table(vb_outliers)
  vb_genes <- merge(genes, vb_outliers, by.y="gene", by.x="anno")
  vb_genes1 <- data.table(vb_genes[, c(2,3,4,15)])
  vb_intergenic <- vb_outliers[!(complete.cases(vb_outliers$gene)),]
  vb_intergenic1 <- data.table(vb_intergenic[, c(1,2,3,12)])
  vb_all <- rbind(vb_genes1, vb_intergenic1, use.names=FALSE)
  setkeyv(vb_all, names(vb_all))
  ans <- foverlaps(mydf, vb_all, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames.x","start.x","end.x"))
  vb_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  vb_outliers1$trait <- "vb"
  
  bb_outliers <- data.table(bb_outliers)
  bb_genes <- merge(genes, bb_outliers, by.y="gene", by.x="anno")
  bb_genes1 <- data.table(bb_genes[, c(2,3,4,15)])
  bb_intergenic <- bb_outliers[!(complete.cases(bb_outliers$gene)),]
  bb_intergenic1 <- data.table(bb_intergenic[, c(1,2,3,12)])
  bb_all <- rbind(bb_genes1, bb_intergenic1, use.names=FALSE)
  setkeyv(bb_all, names(bb_all))
  ans <- foverlaps(mydf, bb_all, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames.x","start.x","end.x"))
  bb_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  bb_outliers1$trait <- "bb"
  
  pc1_outliers <- data.table(pc1_outliers)
  pc1_genes <- merge(genes, pc1_outliers, by.y="gene", by.x="anno")
  pc1_genes1 <- data.table(pc1_genes[, c(2,3,4,15)])
  pc1_intergenic <- pc1_outliers[!(complete.cases(pc1_outliers$gene)),]
  pc1_intergenic1 <- data.table(pc1_intergenic[, c(1,2,3,12)])
  pc1_all <- rbind(pc1_genes1, pc1_intergenic1, use.names=FALSE)
  setkeyv(pc1_all, names(pc1_all))
  ans <- foverlaps(mydf, pc1_all, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames.x","start.x","end.x"))
  pc1_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  pc1_outliers1$trait <- "pc1"
  
  pc2_outliers <- data.table(pc2_outliers)
  pc2_genes <- merge(genes, pc2_outliers, by.y="gene", by.x="anno")
  pc2_genes1 <- data.table(pc2_genes[, c(2,3,4,15)])
  pc2_intergenic <- pc2_outliers[!(complete.cases(pc2_outliers$gene)),]
  pc2_intergenic1 <- data.table(pc2_intergenic[, c(1,2,3,12)])
  pc2_all <- rbind(pc2_genes1, pc2_intergenic1, use.names=FALSE)
  setkeyv(pc2_all, names(pc2_all))
  ans <- foverlaps(mydf, pc2_all, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames.x","start.x","end.x"))
  pc2_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  pc2_outliers1$trait <- "pc2"
  
  cs_outliers <- data.table(cs_outliers)
  cs_genes <- merge(genes, cs_outliers, by.y="gene", by.x="anno")
  cs_genes1 <- data.table(cs_genes[, c(2,3,4,15)])
  cs_intergenic <- cs_outliers[!(complete.cases(cs_outliers$gene)),]
  cs_intergenic1 <- data.table(cs_intergenic[, c(1,2,3,12)])
  cs_all <- rbind(cs_genes1, cs_intergenic1, use.names=FALSE)
  setkeyv(cs_all, names(cs_all))
  ans <- foverlaps(mydf, cs_all, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames.x","start.x","end.x"))
  cs_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  cs_outliers1$trait <- "cs"
  
  ir_outliers <- data.table(ir_outliers)
  ir_genes <- merge(genes, ir_outliers, by.y="gene", by.x="anno")
  ir_genes1 <- data.table(ir_genes[, c(2,3,4,15)])
  ir_intergenic <- ir_outliers[!(complete.cases(ir_outliers$gene)),]
  ir_intergenic1 <- data.table(ir_intergenic[, c(1,2,3,12)])
  ir_all <- rbind(ir_genes1, ir_intergenic1, use.names=FALSE)
  setkeyv(ir_all, names(ir_all))
  ans <- foverlaps(mydf, ir_all, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames.x","start.x","end.x"))
  ir_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  ir_outliers1$trait <- "ir"
  
  kd_outliers <- data.table(kd_outliers)
  kd_genes <- merge(genes, kd_outliers, by.y="gene", by.x="anno")
  kd_genes1 <- data.table(kd_genes[, c(2,3,4,15)])
  kd_intergenic <- kd_outliers[!(complete.cases(kd_outliers$gene)),]
  kd_intergenic1 <- data.table(kd_intergenic[, c(1,2,3,12)])
  kd_all <- rbind(kd_genes1, kd_intergenic1, use.names=FALSE)
  setkeyv(kd_all, names(kd_all))
  ans <- foverlaps(mydf, kd_all, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames.x","start.x","end.x"))
  kd_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  kd_outliers1$trait <- "kd"
  
  tm_outliers <- data.table(tm_outliers)
  tm_genes <- merge(genes, tm_outliers, by.y="gene", by.x="anno")
  tm_genes1 <- data.table(tm_genes[, c(2,3,4,15)])
  tm_intergenic <- tm_outliers[!(complete.cases(tm_outliers$gene)),]
  tm_intergenic1 <- data.table(tm_intergenic[, c(1,2,3,12)])
  tm_all <- rbind(tm_genes1, tm_intergenic1, use.names=FALSE)
  setkeyv(tm_all, names(tm_all))
  ans <- foverlaps(mydf, tm_all, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames.x","start.x","end.x"))
  tm_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  tm_outliers1$trait <- "tm"
  
  
  mydf1 <- rbind(yf_outliers1, fb_outliers1, bb_outliers1,mb_outliers1,rc_outliers1,vb_outliers1,mlb_outliers1,dlb_outliers1,pc1_outliers1,pc2_outliers1,cs_outliers1,ir_outliers1,kd_outliers1,tm_outliers1)
  print(dim(mydf1))
  
  all[[file]] <- mydf1
  
  
  #pdf(paste0(file, ".mahanttan_gwasoutliers_0.999_gwasGENEINETRGENICoutliers_0.95fst_20kb_withoutchr10.pdf"), h=5, w=12)
  #svg(paste0(file, "mahanttan_gwasoutliers_0.999_gwas_0.95fst_topgenes.svg"), h=4, w=12)
  png(paste0("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/fst/",file, ".99.9_pip01_GWAS_SNPoutliers_GENES_0.95_FST_20kb_outliers.png"), width = 2800, height = 1000, res = 300)
  #par(mar = c(1, 1, 1, 1))
  manhattan1(mydf,chr="CHROM",bp="start",p="WEIGHTED_FST",snp="SNP",logp=FALSE,ylab="",
             xlab="", main="", highlight1 = mb_outliers1$i.SNP,
             highlight2=fb_outliers1$i.SNP, highlight3=yf_outliers1$i.SNP, highlight4=bb_outliers1$i.SNP, highlight5=rc_outliers1$i.SNP, highlight6=mlb_outliers1$i.SNP, 
             highlight7=dlb_outliers1$i.SNP, highlight8=vb_outliers1$i.SNP, highlight9=pc1_outliers1$i.SNP, highlight10=pc2_outliers1$i.SNP, highlight11=cs_outliers1$i.SNP,highlight12=ir_outliers1$i.SNP,highlight13=kd_outliers1$i.SNP,highlight14=tm_outliers1$i.SNP, cex=1)
  
  dev.off()
  
  
}

setwd('/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/fst/')
all_anno1 <- rbindlist(all, fill=TRUE, idcol = "comparison")
write.table(all_anno1, "allcomparisons_top_99.9_gwasGENEINTERGENIC_0.95fst_tophits_withoutchr10.txt", quote=F, row.names=F, sep="\t")




##### manhattan plots for pundamilia

manhattan1<-function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c(alpha("gray60",0.5),
                                                                               alpha("gray80", 0.5)), chrlabs = NULL, suggestiveline = -log10(1e-05),
                      genomewideline = -log10(5e-08), highlight1 = NULL, highlight2 = NULL, highlight3 = NULL, highlight4 = NULL, logp = TRUE,
                      highlight5 = NULL, highlight6 = NULL, highlight7 = NULL, highlight8 = NULL, highlight9 = NULL, highlight10 = NULL, highlight11 = NULL, highlight12 = NULL, highlight13 = NULL, highlight14 = NULL,
                      ...)
{
  CHR = BP = P = index = NULL
  if (!(chr %in% names(x)))
    stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x)))
    stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x)))
    stop(paste("Column", p, "not found!"))
  if (!(snp %in% names(x)))
    warning(paste("No SNP column found. OK unless you're trying to highlight."))
  if (!is.numeric(x[[chr]]))
    stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]]))
    stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]]))
    stop(paste(p, "column should be numeric."))
  d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
  if (!is.null(x[[snp]]))
    d = transform(d, SNP = x[[snp]])
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d <- d[order(d$CHR, d$BP), ]
  if (logp) {
    d$logp <- -log10(d$P)
  }
  else {
    d$logp <- d$P
  }
  d$pos = NA
  d$index = NA
  ind = 0
  for (i in unique(d$CHR)) {
    ind = ind + 1
    d[d$CHR == i, ]$index = ind
  }
  nchr = length(unique(d$CHR))
  if (nchr == 1) {
    options(scipen = 999)
    d$pos = d$BP/1e+06
    ticks = floor(length(d$pos))/2 + 1
    xlabel = paste("Chromosome", unique(d$CHR), "position(Mb)")
    labs = ticks
  }
  else {
    lastbase = 0
    ticks = NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos = d[d$index == i, ]$BP
      }
      else {
        lastbase = lastbase + tail(subset(d, index ==
                                            i - 1)$BP, 1)
        d[d$index == i, ]$pos = d[d$index == i, ]$BP +
          lastbase
      }
      ticks = c(ticks, (min(d[d$CHR == i, ]$pos) + max(d[d$CHR ==
                                                           i, ]$pos))/2 + 1)
    }
    xlabel = "Chromosome"
    labs <- unique(d$CHR)
  }
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i",
                   las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0,
                                                                     ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log10))
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in%
                                            names(dotargs)]))
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      }
      else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    }
    else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  if (nchr == 1) {
    axis(1, ...)
  }
  else {
    axis(1, at = ticks, labels = labs, ...)
  }
  col = rep(col, max(d$CHR))
  if (nchr == 1) {
    with(d, points(pos, logp, pch = 20, col = col[1], ...))
  }
  else {
    icol = 1
    for (i in unique(d$index)) {
      with(d[d$index == unique(d$index)[i], ], points(pos,
                                                      logp, col = col[icol], pch = 20, ...))
      icol = icol + 1
    }
  }
  if (suggestiveline)
    abline(h = suggestiveline, col = "blue")
  if (genomewideline)
    abline(h = genomewideline, col = "red")
  if (!is.null(highlight1)) {
    if (any(!(highlight1 %in% d$SNP)))
      warning("You're trying to highlight1 SNPs that don't exist in your results.")
    d.highlight1 = d[which(d$SNP %in% highlight1), ]
    with(d.highlight1, points(pos, logp, col = alpha("#E64B35", 1), pch = 20,
                              ...))
  }
  if (!is.null(highlight2)) {
    if (any(!(highlight2 %in% d$SNP)))
      warning("You're trying to highlight2 SNPs that don't exist in your results.")
    d.highlight2 = d[which(d$SNP %in% highlight2), ]
    with(d.highlight2, points(pos, logp, col = alpha("#E64B35", 1.0), pch = 20,
                              ...))
  }
  if (!is.null(highlight3)) {
    if (any(!(highlight3 %in% d$SNP)))
      warning("You're trying to highlight3 SNPs that don't exist in your results.")
    d.highlight3 = d[which(d$SNP %in% highlight3), ]
    with(d.highlight3, points(pos, logp, col = alpha("#E64B35", 1.0), pch = 20,
                              ...))
  }
  
  if (!is.null(highlight4)) {
    if (any(!(highlight4 %in% d$SNP)))
      warning("You're trying to highlight4 SNPs that don't exist in your results.")
    d.highlight4 = d[which(d$SNP %in% highlight4), ]
    with(d.highlight4, points(pos, logp, col = alpha("#E64B35", 1.0), pch = 20,
                              ...))
  }
  
  if (!is.null(highlight5)) {
    if (any(!(highlight5 %in% d$SNP)))
      warning("You're trying to highlight5 SNPs that don't exist in your results.")
    d.highlight5 = d[which(d$SNP %in% highlight5), ]
    with(d.highlight5, points(pos, logp, col = alpha("#E64B35", 1.0), pch = 20,
                              ...))
  }
  
  if (!is.null(highlight6)) {
    if (any(!(highlight6 %in% d$SNP)))
      warning("You're trying to highlight6 SNPs that don't exist in your results.")
    d.highlight6 = d[which(d$SNP %in% highlight6), ]
    with(d.highlight6, points(pos, logp, col = alpha("#009E73", 1.0), pch = 20,  
                              ...))
  }
  
  if (!is.null(highlight7)) {
    if (any(!(highlight7 %in% d$SNP)))
      warning("You're trying to highlight7 SNPs that don't exist in your results.")
    d.highlight7 = d[which(d$SNP %in% highlight7), ]
    with(d.highlight7, points(pos, logp, col = alpha("#009E73", 1.0), pch = 20,  
                              ...))
  }
  
  if (!is.null(highlight8)) {
    if (any(!(highlight8 %in% d$SNP)))
      warning("You're trying to highlight8 SNPs that don't exist in your results.")
    d.highlight8 = d[which(d$SNP %in% highlight8), ]
    with(d.highlight8, points(pos, logp, col = alpha("#009E73", 1.0), pch = 20,  
                              ...))
  }
  
  if (!is.null(highlight9)) {
    if (any(!(highlight9 %in% d$SNP)))
      warning("You're trying to highlight9 SNPs that don't exist in your results.")
    d.highlight9 = d[which(d$SNP %in% highlight9), ]
    with(d.highlight9, points(pos, logp, col = alpha("#4DBBD5", 1.0), pch = 20,  
                              ...))
  }
  
  if (!is.null(highlight10)) {
    if (any(!(highlight10 %in% d$SNP)))
      warning("You're trying to highlight10 SNPs that don't exist in your results.")
    d.highlight10 = d[which(d$SNP %in% highlight10), ]
    with(d.highlight10, points(pos, logp, col = alpha("#4DBBD5", 1), pch = 20,  
                               ...))
  }
  
  if (!is.null(highlight11)) {
    if (any(!(highlight11 %in% d$SNP)))
      warning("You're trying to highlight10 SNPs that don't exist in your results.")
    d.highlight11 = d[which(d$SNP %in% highlight11), ]
    with(d.highlight11, points(pos, logp, col = alpha("#4DBBD5", 1), pch = 20,  
                               ...))
  }
  
  if (!is.null(highlight12)) {
    if (any(!(highlight12 %in% d$SNP)))
      warning("You're trying to highlight10 SNPs that don't exist in your results.")
    d.highlight12 = d[which(d$SNP %in% highlight12), ]
    with(d.highlight12, points(pos, logp, col = alpha("#4DBBD5", 1), pch = 20,  
                               ...))
  }
  if (!is.null(highlight13)) {
    if (any(!(highlight13 %in% d$SNP)))
      warning("You're trying to highlight10 SNPs that don't exist in your results.")
    d.highlight13 = d[which(d$SNP %in% highlight13), ]
    with(d.highlight13, points(pos, logp, col = alpha("#4DBBD5", 1), pch = 20,  
                               ...))
  }
  if (!is.null(highlight14)) {
    if (any(!(highlight14 %in% d$SNP)))
      warning("You're trying to highlight10 SNPs that don't exist in your results.")
    d.highlight14 = d[which(d$SNP %in% highlight14), ]
    with(d.highlight14, points(pos, logp, col = alpha("#4DBBD5", 1), pch = 20,  
                               ...))
  }
  
  
}
library(rlist)
setwd('/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_repeatcolour/fst')

counter <- 0
all <- list()
files <- c("pnye_makobe_vs_ppun_makobe_20kb.windowed.weir.fst","pnye_ruti_vs_ppun_ruti_20kb.windowed.weir.fst","pnye_kiss_vs_ppun_kiss_20kb.windowed.weir.fst","pnye_python_vs_ppun_python_20kb.windowed.weir.fst")

pdf("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/fst/2025_figure_4_mahanttan_99.9_GWAS_SNPoutliers_0.95_FST_20kb_outliers_v2.pdf")
par(mfrow=c(4,1), oma = c(0, 0, 0, 0), mai = c(0.3, 0.5, 0.1, 0.1))

for (file in files) {
  counter = counter +1
  fst <- read.table(file, header=T)
  head(fst)
  fstsubset<-fst[complete.cases(fst),]
  fstsubset$WEIGHTED_FST <- replace(fstsubset$WEIGHTED_FST, fstsubset$WEIGHTED_FST < 0,0)
  fstsubset_without_chr10 <- fstsubset[fstsubset$CHROM != "chr10",]
  cutoff <- quantile(fstsubset_without_chr10$WEIGHTED_FST, prob=0.95)
  #cutoff <- quantile(fstsubset$WEIGHTED_FST, prob=0.95)
  SNP<-c(1:(nrow(fstsubset)))
  mydf<-data.frame(SNP,fstsubset)
  mydf$CHROM <- gsub('chr', '', mydf$CHROM)
  mydf$CHROM <- as.numeric(mydf$CHROM)
  mydf$SNP <- paste0(mydf$CHROM, "_", mydf$BIN_START)
  mydf$seqnames <- paste0("chr", mydf$CHROM)
  colnames(mydf) <- c("SNP","CHROM","start","end","N_VARIANTS","WEIGHTED_FST", "MEAN_FST","seqnames")
  mydf <- mydf[, c(8,3,4,1,2,5,6,7)]
  
  mydf <- data.table(mydf)
  setkeyv(mydf, names(mydf))
  
  yf_outliers <- data.table(yf_outliers)
  yf_genes <- merge(genes, yf_outliers, by.y="gene", by.x="anno")
  yf_genes1 <- data.table(yf_genes[, c(2,3,4,15)])
  yf_intergenic <- yf_outliers[!(complete.cases(yf_outliers$gene)),]
  yf_intergenic1 <- data.table(yf_intergenic[, c(1,2,3,12)])
  yf_all <- rbind(yf_genes1, yf_intergenic1, use.names=FALSE)
  setkeyv(yf_all, names(yf_all))
  ans <- foverlaps(mydf, yf_all, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames.x","start.x","end.x"))
  yf_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  yf_outliers1$trait <- "YF"
  
  
  fb_outliers <- data.table(fb_outliers)
  fb_genes <- merge(genes, fb_outliers, by.y="gene", by.x="anno")
  fb_genes1 <- data.table(fb_genes[, c(2,3,4,15)])
  fb_intergenic <- fb_outliers[!(complete.cases(fb_outliers$gene)),]
  fb_intergenic1 <- data.table(fb_intergenic[, c(1,2,3,12)])
  fb_all <- rbind(fb_genes1, fb_intergenic1, use.names=FALSE)
  setkeyv(fb_all, names(fb_all))
  ans <- foverlaps(mydf, fb_all, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames.x","start.x","end.x"))
  fb_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  fb_outliers1$trait <- "fb"
  
  mb_outliers <- data.table(mb_outliers)
  mb_genes <- merge(genes, mb_outliers, by.y="gene", by.x="anno")
  mb_genes1 <- data.table(mb_genes[, c(2,3,4,15)])
  mb_intergenic <- mb_outliers[!(complete.cases(mb_outliers$gene)),]
  mb_intergenic1 <- data.table(mb_intergenic[, c(1,2,3,12)])
  mb_all <- rbind(mb_genes1, mb_intergenic1, use.names=FALSE)
  setkeyv(mb_all, names(mb_all))
  ans <- foverlaps(mydf, mb_all, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames.x","start.x","end.x"))
  mb_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  mb_outliers1$trait <- "mb"
  
  rc_outliers <- data.table(rc_outliers)
  rc_genes <- merge(genes, rc_outliers, by.y="gene", by.x="anno")
  rc_genes1 <- data.table(rc_genes[, c(2,3,4,15)])
  rc_intergenic <- rc_outliers[!(complete.cases(rc_outliers$gene)),]
  rc_intergenic1 <- data.table(rc_intergenic[, c(1,2,3,12)])
  rc_all <- rbind(rc_genes1, rc_intergenic1, use.names=FALSE)
  setkeyv(rc_all, names(rc_all))
  ans <- foverlaps(mydf, rc_all, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames.x","start.x","end.x"))
  rc_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  rc_outliers1$trait <- "rc"
  
  mlb_outliers <- data.table(mlb_outliers)
  mlb_genes <- merge(genes, mlb_outliers, by.y="gene", by.x="anno")
  mlb_genes1 <- data.table(mlb_genes[, c(2,3,4,15)])
  mlb_intergenic <- mlb_outliers[!(complete.cases(mlb_outliers$gene)),]
  mlb_intergenic1 <- data.table(mlb_intergenic[, c(1,2,3,12)])
  mlb_all <- rbind(mlb_genes1, mlb_intergenic1, use.names=FALSE)
  setkeyv(mlb_all, names(mlb_all))
  ans <- foverlaps(mydf, mlb_all, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames.x","start.x","end.x"))
  mlb_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  mlb_outliers1$trait <- "mlb"
  
  dlb_outliers <- data.table(dlb_outliers)
  dlb_genes <- merge(genes, dlb_outliers, by.y="gene", by.x="anno")
  dlb_genes1 <- data.table(dlb_genes[, c(2,3,4,15)])
  dlb_intergenic <- dlb_outliers[!(complete.cases(dlb_outliers$gene)),]
  dlb_intergenic1 <- data.table(dlb_intergenic[, c(1,2,3,12)])
  dlb_all <- rbind(dlb_genes1, dlb_intergenic1, use.names=FALSE)
  setkeyv(dlb_all, names(dlb_all))
  ans <- foverlaps(mydf, dlb_all, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames.x","start.x","end.x"))
  dlb_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  dlb_outliers1$trait <- "dlb"
  
  vb_outliers <- data.table(vb_outliers)
  vb_genes <- merge(genes, vb_outliers, by.y="gene", by.x="anno")
  vb_genes1 <- data.table(vb_genes[, c(2,3,4,15)])
  vb_intergenic <- vb_outliers[!(complete.cases(vb_outliers$gene)),]
  vb_intergenic1 <- data.table(vb_intergenic[, c(1,2,3,12)])
  vb_all <- rbind(vb_genes1, vb_intergenic1, use.names=FALSE)
  setkeyv(vb_all, names(vb_all))
  ans <- foverlaps(mydf, vb_all, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames.x","start.x","end.x"))
  vb_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  vb_outliers1$trait <- "vb"
  
  bb_outliers <- data.table(bb_outliers)
  bb_genes <- merge(genes, bb_outliers, by.y="gene", by.x="anno")
  bb_genes1 <- data.table(bb_genes[, c(2,3,4,15)])
  bb_intergenic <- bb_outliers[!(complete.cases(bb_outliers$gene)),]
  bb_intergenic1 <- data.table(bb_intergenic[, c(1,2,3,12)])
  bb_all <- rbind(bb_genes1, bb_intergenic1, use.names=FALSE)
  setkeyv(bb_all, names(bb_all))
  ans <- foverlaps(mydf, bb_all, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames.x","start.x","end.x"))
  bb_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  bb_outliers1$trait <- "bb"
  
  pc1_outliers <- data.table(pc1_outliers)
  pc1_genes <- merge(genes, pc1_outliers, by.y="gene", by.x="anno")
  pc1_genes1 <- data.table(pc1_genes[, c(2,3,4,15)])
  pc1_intergenic <- pc1_outliers[!(complete.cases(pc1_outliers$gene)),]
  pc1_intergenic1 <- data.table(pc1_intergenic[, c(1,2,3,12)])
  pc1_all <- rbind(pc1_genes1, pc1_intergenic1, use.names=FALSE)
  setkeyv(pc1_all, names(pc1_all))
  ans <- foverlaps(mydf, pc1_all, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames.x","start.x","end.x"))
  pc1_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  pc1_outliers1$trait <- "pc1"
  
  pc2_outliers <- data.table(pc2_outliers)
  pc2_genes <- merge(genes, pc2_outliers, by.y="gene", by.x="anno")
  pc2_genes1 <- data.table(pc2_genes[, c(2,3,4,15)])
  pc2_intergenic <- pc2_outliers[!(complete.cases(pc2_outliers$gene)),]
  pc2_intergenic1 <- data.table(pc2_intergenic[, c(1,2,3,12)])
  pc2_all <- rbind(pc2_genes1, pc2_intergenic1, use.names=FALSE)
  setkeyv(pc2_all, names(pc2_all))
  ans <- foverlaps(mydf, pc2_all, nomatch=0L, by.x=c("seqnames","start","end"), by.y=c("seqnames.x","start.x","end.x"))
  pc2_outliers1 <- ans[ans$MEAN_FST >= cutoff,]
  pc2_outliers1$trait <- "pc2"
  
  mydf1 <- rbind(yf_outliers1, fb_outliers1, bb_outliers1,mb_outliers1,rc_outliers1,vb_outliers1,mlb_outliers1,dlb_outliers1,pc1_outliers1,pc2_outliers1)
  print(dim(mydf1))
  
  all[[file]] <- mydf1
  
  
  #png(paste0(file, ".mahanttan_gwasoutliers_0.999_gwasGENEINETRGENICoutliers_0.95fst_20kb_withoutchr10.png"), h=4, w=10)
  #svg(paste0(file, "mahanttan_gwasoutliers_0.999_gwas_0.95fst_topgenes.svg"), h=5, w=14)
  
  
  manhattan1(mydf,chr="CHROM",bp="start",p="WEIGHTED_FST",snp="SNP",logp=FALSE,ylab="",
             xlab="", main="", highlight1 = mb_outliers1$i.SNP,
             highlight2=fb_outliers1$i.SNP, highlight3=yf_outliers1$i.SNP, highlight4=bb_outliers1$i.SNP, highlight5=rc_outliers1$i.SNP, highlight6=mlb_outliers1$i.SNP, 
             highlight7=dlb_outliers1$i.SNP, highlight8=vb_outliers1$i.SNP,highlight9=pc1_outliers1$i.SNP, highlight10=pc2_outliers1$i.SNP)
  
  #dev.off()
  
  
}
dev.off()

#all_anno1 <- rbindlist(all, fill=TRUE, idcol = "comparison")
#write.table(all_anno1, "allcomparisons_top_99.9_gwas_0.95fst_tophits_withoutchr10_genelevel.txt", quote=F, row.names=F, sep="\t")



