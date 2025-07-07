
#pooja.singh09@gmail.com
#this plots top hits that are in the top 99.9 and pip > 0.01 for each trait

while (dev.cur()>1) dev.off()

options(echo=TRUE)
library(zoo) 
library(basicPlotteR)
library(dplyr)
library(svglite)

## main figure for paper Figure 3 pattern main

#png(file="/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/manhattanplots/figure3_main_colpattern_manhattan.png", width=8,height=8, units="in",res=300)
svglite(file="/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/manhattanplots/figure3_main_colpattern_manhattan.svg", width=8,height=5)

par( mai=c(0.3, 0.6, 0.1, 0.1), mfrow=c(5,1))


#### mlb
# Reads input parameters
args<-commandArgs(trailingOnly=T)
#trait=args[1]
trait=c("MidLateralBand")

setwd('/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_repeatcolour/bslmm/')

# library to speed up loading of big tables
library(data.table)

# Load parameters output

params<-fread(paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".mean.params.anno.txt", sep=""),header=T,sep="\t", data.table=F)


params$pos <- params$start
params$chr <- as.numeric(gsub("chr","", params$seqnames))


# get list of linkage groups/chromosomes
chrs<-sort(as.numeric(unique(params$chr)))

# add linkage group column (chr)
chr<-gsub("lg|_.+","",params$chr)
params["chr"]<-chr


# sort by linkage group and position
params.sort<-params[order(as.numeric(params$chr), params$pos),]
params.sort$lg <- paste("LG",params.sort$chr, sep="")
params.sort$rs <- paste(params.sort$lg, params.sort$pos, sep="-")
params.sort$gene_name <- sapply(strsplit(params.sort$gene,split = ';'), `[`, 2)
params.sort$gene_name <- gsub('Name=','',params.sort$gene_name)
params.sort$eff <- abs(params.sort$beta*params.sort$gamma)


params.sort1 <- params.sort[abs(params.sort$beta) > 0,] # only variants with some effect
params.sort <- params.sort1
######################### plot of sum_gamma



# set up empty plot
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.ylim.png", sep=""), width=11.7,height=8.3,units="in",res=200)
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.png", sep=""), width=12,height=6)


# set up empty plot
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0,max(params.sort1$gamma)),ylab=paste0("Top SNPs ", trait, " (PIP)"),xlab="linkage group", xaxt="n")


# plot grey bands for chromosome/linkage groups

start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(params.sort[params.sort$chr==ch,])
  cat ("CH: ", ch, "\n")
  colour<-adjustcolor("lightgrey", alpha.f=0.1)
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=colour, border=colour)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  
  start<-start+size
}

# Add x axis labels
#axis(side=1,at=lab.pos,labels=chrs,tick=F)


# plot PIP for insignificant variants

#x<-match(params.sort$gamma[params.sort$gamma<threshold],params.sort$gamma)
top <- params.sort[quantile(params.sort$gamma,0.9999)[[1]] & params.sort$gamma > 0.01,]
x<- match(params.sort$gamma[!params.sort$id %in% top$id],params.sort$gamma)

# PIP
y<-params.sort$gamma[!params.sort$id %in% top$id]

# sparse effect size, used for dot size
z<-(params.sort[!params.sort$id %in% top$id,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))

new <- data.frame(x=x,y=y,z=z)
set.seed(123)  # Set a seed for reproducibility
new1 <- new[sample(nrow(new), 20), ]
symbols(new1$x,new1$y,circles=new1$z, bg="grey90",inches=1/26,fg=NULL,add=T)


# plot threshold line
#abline(h=threshold,lty=3,col="black")

# rank of high PIP variants across linkage groups with gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
z<-(params.sort[params.sort$id %in% top$id & !is.na(params.sort$gene_name),])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  


symbols(x,y, circles=z, bg=adjustcolor("#009E73",alpha.f = 0.8), inches=1/26,fg="black",add=T)

# rank of high PIP variants across linkage groups without gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
#z<-(params.sort[params.sort$gamma>=threshold,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#009E73",alpha.f = 0.8), inches=1/26,fg=NA,add=T)

# add label high PIP variants
#text(x,y,labels=params.sort$gene[params.sort$gamma>=threshold], adj=c(0,0), cex=0.4)
top_20_rows <- head(top[!is.na(top$gene_name), ][order(-top$gamma[!is.na(top$gene_name)]), ], 20)


labels <- top_20_rows$gene_name
x1<-match(top_20_rows$gamma,params.sort$gamma)
y1<-top_20_rows$gamma

addTextLabels(x1,y1,tolower(labels), col.label="black", cex.label = 1.1,  col.background = NULL, lty=1)


#### dlb

# Reads input parameters
args<-commandArgs(trailingOnly=T)
#trait=args[1]
trait=c("DorsalLateralBand")

# library to speed up loading of big tables
library(data.table)

# Load parameters output

params<-fread(paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".mean.params.anno.txt", sep=""),header=T,sep="\t", data.table=F)


params$pos <- params$start
params$chr <- as.numeric(gsub("chr","", params$seqnames))


# get list of linkage groups/chromosomes
chrs<-sort(as.numeric(unique(params$chr)))

# add linkage group column (chr)
chr<-gsub("lg|_.+","",params$chr)
params["chr"]<-chr


# sort by linkage group and position
params.sort<-params[order(as.numeric(params$chr), params$pos),]
params.sort$lg <- paste("LG",params.sort$chr, sep="")
params.sort$rs <- paste(params.sort$lg, params.sort$pos, sep="-")
params.sort$gene_name <- sapply(strsplit(params.sort$gene,split = ';'), `[`, 2)
params.sort$gene_name <- gsub('Name=','',params.sort$gene_name)
params.sort$eff <- abs(params.sort$beta*params.sort$gamma)

params.sort1 <- params.sort[abs(params.sort$beta) > 0,] # only variants with some effect
params.sort <- params.sort1
######################### plot of sum_gamma



# set up empty plot
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.ylim.png", sep=""), width=11.7,height=8.3,units="in",res=200)
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.png", sep=""), width=12,height=6)


# set up empty plot
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0,max(params.sort$gamma)),ylab=paste0("Top SNPs ", trait, " (PIP)"),xlab="linkage group", xaxt="n")


# plot grey bands for chromosome/linkage groups

start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(params.sort[params.sort$chr==ch,])
  cat ("CH: ", ch, "\n")
  colour<-adjustcolor("lightgrey", alpha.f=0.1)
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=colour, border=colour)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  
  start<-start+size
}

# Add x axis labels
axis(side=1,at=lab.pos,labels=chrs,tick=F)


# plot PIP for insignificant variants

#x<-match(params.sort$gamma[params.sort$gamma<threshold],params.sort$gamma)
top <- params.sort[quantile(params.sort$gamma,0.9999)[[1]] & params.sort$gamma > 0.01,]
x<- match(params.sort$gamma[!params.sort$id %in% top$id],params.sort$gamma)

# PIP
y<-params.sort$gamma[!params.sort$id %in% top$id]

# sparse effect size, used for dot size
z<-(params.sort[!params.sort$id %in% top$id,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))
new <- data.frame(x=x,y=y,z=z)
set.seed(123)  # Set a seed for reproducibility
new1 <- new[sample(nrow(new), 20), ]
symbols(new1$x,new1$y,circles=new1$z, bg="grey90",inches=1/26,fg=NULL,add=T)

# plot threshold line
#abline(h=threshold,lty=3,col="black")

# rank of high PIP variants across linkage groups with gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
z<-(params.sort[params.sort$id %in% top$id & !is.na(params.sort$gene_name),])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#009E73",alpha.f = 0.8), inches=1/26,fg="black",add=T)

# rank of high PIP variants across linkage groups without gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
#z<-(params.sort[params.sort$gamma>=threshold,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#009E73",alpha.f = 0.8), inches=1/26,fg=NA,add=T)

# add label high PIP variants
#text(x,y,labels=params.sort$gene[params.sort$gamma>=threshold], adj=c(0,0), cex=0.4)
top_20_rows <- head(top[!is.na(top$gene_name), ][order(-top$gamma[!is.na(top$gene_name)]), ], 20)


labels <- top_20_rows$gene_name
x1<-match(top_20_rows$gamma,params.sort$gamma)
y1<-top_20_rows$gamma

addTextLabels(x1,y1,tolower(labels), col.label="black", cex.label = 1.1,  col.background = NULL, lty=1)


#####fb
# Reads input parameters
args<-commandArgs(trailingOnly=T)
#trait=args[1]
trait=c("Flameback")

# library to speed up loading of big tables
library(data.table)

# Load parameters output

params<-fread(paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".mean.params.anno.txt", sep=""),header=T,sep="\t", data.table=F)


params$pos <- params$start
params$chr <- as.numeric(gsub("chr","", params$seqnames))


# get list of linkage groups/chromosomes
chrs<-sort(as.numeric(unique(params$chr)))

# add linkage group column (chr)
chr<-gsub("lg|_.+","",params$chr)
params["chr"]<-chr


# sort by linkage group and position
params.sort<-params[order(as.numeric(params$chr), params$pos),]
params.sort$lg <- paste("LG",params.sort$chr, sep="")
params.sort$rs <- paste(params.sort$lg, params.sort$pos, sep="-")
params.sort$gene_name <- sapply(strsplit(params.sort$gene,split = ';'), `[`, 2)
params.sort$gene_name <- gsub('Name=','',params.sort$gene_name)
params.sort$eff <- abs(params.sort$beta*params.sort$gamma)

params.sort1 <- params.sort[abs(params.sort$beta) > 0,] # only variants with some effect
params.sort <- params.sort1
######################### plot of sum_gamma



# set up empty plot
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.ylim.png", sep=""), width=11.7,height=8.3,units="in",res=200)
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.png", sep=""), width=12,height=6)


# set up empty plot
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0,max(params.sort$gamma)),ylab=paste0("Top SNPs ", trait, " (PIP)"),xlab="linkage group", xaxt="n")


# plot grey bands for chromosome/linkage groups

start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(params.sort[params.sort$chr==ch,])
  cat ("CH: ", ch, "\n")
  colour<-adjustcolor("lightgrey", alpha.f=0.1)
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=colour, border=colour)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  
  start<-start+size
}

# Add x axis labels
#axis(side=1,at=lab.pos,labels=chrs,tick=F)


# plot PIP for insignificant variants

#x<-match(params.sort$gamma[params.sort$gamma<threshold],params.sort$gamma)
top <- params.sort[quantile(params.sort$gamma,0.9999)[[1]] & params.sort$gamma > 0.01,]

x<- match(params.sort$gamma[!params.sort$id %in% top$id],params.sort$gamma)

# PIP
y<-params.sort$gamma[!params.sort$id %in% top$id]

# sparse effect size, used for dot size
z<-(params.sort[!params.sort$id %in% top$id,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))
new <- data.frame(x=x,y=y,z=z)
set.seed(123)  # Set a seed for reproducibility
new1 <- new[sample(nrow(new), 20), ]
symbols(new1$x,new1$y,circles=new1$z, bg="grey90",inches=1/26,fg=NULL,add=T)

# plot threshold line
#abline(h=threshold,lty=3,col="black")

# rank of high PIP variants across linkage groups with gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
z<-(params.sort[params.sort$id %in% top$id & !is.na(params.sort$gene_name),])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#E64B35",alpha.f = 0.8), inches=1/26,fg="black",add=T)

# rank of high PIP variants across linkage groups without gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
#z<-(params.sort[params.sort$gamma>=threshold,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#E64B35",alpha.f = 0.8), inches=1/26,fg=NA,add=T)

# add label high PIP variants
#text(x,y,labels=params.sort$gene[params.sort$gamma>=threshold], adj=c(0,0), cex=0.4)
top_20_rows <- head(top[!is.na(top$gene_name), ][order(-top$gamma[!is.na(top$gene_name)]), ], 20)


labels <- top_20_rows$gene_name
x1<-match(top_20_rows$gamma,params.sort$gamma)
y1<-top_20_rows$gamma

addTextLabels(x1,y1,tolower(labels), col.label="black", cex.label = 1.1,  col.background = NULL, lty=1)


#dev.off()

#tophits <- params.sort[params.sort$gamma>=threshold,]
#write.table(tophits, paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, "99.9.PIPoutlier.SNP.txt", sep=""), quote=F, row.names=F, sep="\t")

#### 2d


# Reads input parameters
args<-commandArgs(trailingOnly=T)
#trait=args[1]
trait=c("YellowFlank")

# library to speed up loading of big tables
library(data.table)

# Load parameters output

params<-fread(paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".mean.params.anno.txt", sep=""),header=T,sep="\t", data.table=F)


params$pos <- params$start
params$chr <- as.numeric(gsub("chr","", params$seqnames))


# get list of linkage groups/chromosomes
chrs<-sort(as.numeric(unique(params$chr)))

# add linkage group column (chr)
chr<-gsub("lg|_.+","",params$chr)
params["chr"]<-chr


# sort by linkage group and position
params.sort<-params[order(as.numeric(params$chr), params$pos),]
params.sort$lg <- paste("LG",params.sort$chr, sep="")
params.sort$rs <- paste(params.sort$lg, params.sort$pos, sep="-")
params.sort$gene_name <- sapply(strsplit(params.sort$gene,split = ';'), `[`, 2)
params.sort$gene_name <- gsub('Name=','',params.sort$gene_name)
params.sort$eff <- abs(params.sort$beta*params.sort$gamma)

params.sort1 <- params.sort[abs(params.sort$beta) > 0,] # only variants with some effect
params.sort <- params.sort1
######################### plot of sum_gamma



# set up empty plot
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.ylim.png", sep=""), width=11.7,height=8.3,units="in",res=200)
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.png", sep=""), width=12,height=6)


# set up empty plot
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0,max(params.sort$gamma)),ylab=paste0("Top SNPs ", trait, " (PIP)"),xlab="linkage group", xaxt="n")


# plot grey bands for chromosome/linkage groups

start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(params.sort[params.sort$chr==ch,])
  cat ("CH: ", ch, "\n")
  colour<-adjustcolor("lightgrey", alpha.f=0.1)
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=colour, border=colour)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  
  start<-start+size
}

# Add x axis labels
#axis(side=1,at=lab.pos,labels=chrs,tick=F)


# plot PIP for insignificant variants

#x<-match(params.sort$gamma[params.sort$gamma<threshold],params.sort$gamma)
top <- params.sort[quantile(params.sort$gamma,0.9999)[[1]] & params.sort$gamma > 0.01,]
dim(params.sort)

x<- match(params.sort$gamma[!params.sort$id %in% top$id],params.sort$gamma)

# PIP
y<-params.sort$gamma[!params.sort$id %in% top$id]

# sparse effect size, used for dot size
z<-(params.sort[!params.sort$id %in% top$id,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))
new <- data.frame(x=x,y=y,z=z)
set.seed(123)  # Set a seed for reproducibility
new1 <- new[sample(nrow(new), 20), ]
symbols(new1$x,new1$y,circles=new1$z, bg="grey90",inches=1/26,fg=NULL,add=T)

# plot threshold line
#abline(h=threshold,lty=3,col="black")

# rank of high PIP variants across linkage groups with gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
z<-(params.sort[params.sort$id %in% top$id & !is.na(params.sort$gene_name),])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#E64B35",alpha.f = 0.8), inches=1/26,fg="black",add=T)

# rank of high PIP variants across linkage groups without gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
#z<-(params.sort[params.sort$gamma>=threshold,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#E64B35",alpha.f = 0.8), inches=1/26,fg=NA,add=T)

# add label high PIP variants
#text(x,y,labels=params.sort$gene[params.sort$gamma>=threshold], adj=c(0,0), cex=0.4)
top_20_rows <- head(top[!is.na(top$gene_name), ][order(-top$gamma[!is.na(top$gene_name)]), ], 20)


labels <- top_20_rows$gene_name
x1<-match(top_20_rows$gamma,params.sort$gamma)
y1<-top_20_rows$gamma

addTextLabels(x1,y1,tolower(labels), col.label="black", cex.label = 1.1,  col.background = NULL, lty=1)

#dev.off()

#tophits <- params.sort[params.sort$gamma>=threshold,]
#write.table(tophits, paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, "99.9.PIPoutlier.SNP.txt", sep=""), quote=F, row.names=F, sep="\t")



dev.off()

################################################### Figure 3 shape main


#png(file="/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/manhattanplots/figure3_main_shape_manhattan.png", width=8,height=8, units="in",res=300)
svglite(file="/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/manhattanplots/figure3_main_shape_manhattan.svg", width=8,height=5)

par( mai=c(0.3, 0.6, 0.1, 0.1), mfrow=c(4,1))


########### cs


setwd('/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_bcftools_pnye3_dentition_repeat/bslmm/') # to read in files

### multiplot 1a
# Reads input parameters
args<-commandArgs(trailingOnly=T)
#trait=args[1]
trait=c("CuspShape")


# library to speed up loading of big tables
library(data.table)

# Load parameters output

params<-fread(paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.", trait, ".mean.params.anno.txt", sep=""),header=T,sep="\t", data.table=F)


params$pos <- params$start
params$chr <- as.numeric(gsub("chr","", params$seqnames))


# get list of linkage groups/chromosomes
chrs<-sort(as.numeric(unique(params$chr)))

# add linkage group column (chr)
chr<-gsub("lg|_.+","",params$chr)
params["chr"]<-chr


# sort by linkage group and position
params.sort<-params[order(as.numeric(params$chr), params$pos),]
params.sort$lg <- paste("LG",params.sort$chr, sep="")
params.sort$rs <- paste(params.sort$lg, params.sort$pos, sep="-")
params.sort$gene_name <- sapply(strsplit(params.sort$gene,split = ';'), `[`, 2)
params.sort$gene_name <- gsub('Name=','',params.sort$gene_name)
params.sort$eff <- abs(params.sort$beta*params.sort$gamma)

params.sort1 <- params.sort[abs(params.sort$beta) > 0,] # only variants with some effect
params.sort <- params.sort1
######################### plot of sum_gamma



# set up empty plot
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.ylim.png", sep=""), width=11.7,height=8.3,units="in",res=200)
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.png", sep=""), width=12,height=6)


# set up empty plot
#plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0,max(params.sort$gamma)),ylab=paste0("Top SNPs ", trait, " (PIP)"), xaxt="n")
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0,1),ylab=paste0("Top SNPs ", trait, " (PIP)"), xaxt="n")


# plot grey bands for chromosome/linkage groups

start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(params.sort[params.sort$chr==ch,])
  cat ("CH: ", ch, "\n")
  colour<-adjustcolor("lightgrey", alpha.f=0.1)
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=colour, border=colour)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  
  start<-start+size
}

# Add x axis labels
#axis(side=1,at=lab.pos,labels=chrs,tick=F)


# highligh high PIP variants

#threshold <- 0.01
#threshold <- 0.01 
#abline(h=threshold,lty=3,col="black")

# plot PIP for insignificant variants

#x<-match(params.sort$gamma[params.sort$gamma<threshold],params.sort$gamma)
top <- params.sort[quantile(params.sort$gamma,0.9999)[[1]] & params.sort$gamma > 0.01,]
x<- match(params.sort$gamma[!params.sort$id %in% top$id],params.sort$gamma)

# PIP
y<-params.sort$gamma[!params.sort$id %in% top$id]

# sparse effect size, used for dot size
z<-(params.sort[!params.sort$id %in% top$id,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))
new <- data.frame(x=x,y=y,z=z)
set.seed(123)  # Set a seed for reproducibility
new1 <- new[sample(nrow(new), 20), ]
symbols(new1$x,new1$y,circles=new1$z, bg="grey90",inches=1/26,fg=NULL,add=T)


# plot threshold line
#abline(h=threshold,lty=3,col="black")

# rank of high PIP variants across linkage groups with gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
z<-(params.sort[params.sort$id %in% top$id & !is.na(params.sort$gene_name),])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#4DBBD5",alpha.f = 0.8), inches=1/22,fg="black",add=T)

# rank of high PIP variants across linkage groups without gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
#z<-(params.sort[params.sort$gamma>=threshold,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#4DBBD5",alpha.f = 0.8), inches=1/22,fg=NA,add=T)

# add label high PIP variants
#text(x,y,labels=params.sort$gene[params.sort$gamma>=threshold], adj=c(0,0), cex=0.4)
top_20_rows <- head(top[!is.na(top$gene_name), ][order(-top$gamma[!is.na(top$gene_name)]), ], 20)


labels <- top_20_rows$gene_name
x1<-match(top_20_rows$gamma,params.sort$gamma)
y1<-top_20_rows$gamma

addTextLabels(x1,y1,tolower(labels), col.label="black", cex.label = 1.1,  col.background = NULL, lty=1)

#dev.off()


#tophits <- params.sort[params.sort$gamma>=threshold,]
#write.table(tophits, paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, "99.9.PIPoutlier.SNP.txt", sep=""), quote=F, row.names=F, sep="\t")



###### kd


# Reads input parameters
args<-commandArgs(trailingOnly=T)
#trait=args[1]
trait=c("LPJkeelDepth")


# library to speed up loading of big tables
library(data.table)

# Load parameters output

params<-fread(paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.", trait, ".mean.params.anno.txt", sep=""),header=T,sep="\t", data.table=F)


params$pos <- params$start
params$chr <- as.numeric(gsub("chr","", params$seqnames))


# get list of linkage groups/chromosomes
chrs<-sort(as.numeric(unique(params$chr)))

# add linkage group column (chr)
chr<-gsub("lg|_.+","",params$chr)
params["chr"]<-chr


# sort by linkage group and position
params.sort<-params[order(as.numeric(params$chr), params$pos),]
params.sort$lg <- paste("LG",params.sort$chr, sep="")
params.sort$rs <- paste(params.sort$lg, params.sort$pos, sep="-")
params.sort$gene_name <- sapply(strsplit(params.sort$gene,split = ';'), `[`, 2)
params.sort$gene_name <- gsub('Name=','',params.sort$gene_name)
params.sort$eff <- abs(params.sort$beta*params.sort$gamma)

params.sort1 <- params.sort[abs(params.sort$beta) > 0,] # only variants with some effect
params.sort <- params.sort1
######################### plot of sum_gamma



# set up empty plot
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.ylim.png", sep=""), width=11.7,height=8.3,units="in",res=200)
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.png", sep=""), width=12,height=6)


# set up empty plot
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0,max(params.sort$gamma)),ylab=paste0("Top SNPs ", trait, " (PIP)"), xaxt="n")


# plot grey bands for chromosome/linkage groups

start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(params.sort[params.sort$chr==ch,])
  cat ("CH: ", ch, "\n")
  colour<-adjustcolor("lightgrey", alpha.f=0.1)
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=colour, border=colour)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  
  start<-start+size
}

# Add x axis labels
#axis(side=1,at=lab.pos,labels=chrs,tick=F)


# plot PIP for insignificant variants

#x<-match(params.sort$gamma[params.sort$gamma<threshold],params.sort$gamma)
top <- params.sort[quantile(params.sort$gamma,0.9999)[[1]] & params.sort$gamma > 0.01,]
x<- match(params.sort$gamma[!params.sort$id %in% top$id],params.sort$gamma)

# PIP
y<-params.sort$gamma[!params.sort$id %in% top$id]

# sparse effect size, used for dot size
z<-(params.sort[!params.sort$id %in% top$id,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))
new <- data.frame(x=x,y=y,z=z)
set.seed(123)  # Set a seed for reproducibility
new1 <- new[sample(nrow(new), 20), ]
symbols(new1$x,new1$y,circles=new1$z, bg="grey90",inches=1/26,fg=NULL,add=T)

# plot threshold line
#abline(h=threshold,lty=3,col="black")

# rank of high PIP variants across linkage groups with gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
z<-(params.sort[params.sort$id %in% top$id & !is.na(params.sort$gene_name),])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#4DBBD5",alpha.f = 0.8), inches=1/26,fg="black",add=T)

# rank of high PIP variants across linkage groups without gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
#z<-(params.sort[params.sort$gamma>=threshold,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#4DBBD5",alpha.f = 0.8), inches=1/26,fg=NA,add=T)

# add label high PIP variants
#text(x,y,labels=params.sort$gene[params.sort$gamma>=threshold], adj=c(0,0), cex=0.4)
top_20_rows <- head(top[!is.na(top$gene_name), ][order(-top$gamma[!is.na(top$gene_name)]), ], 20)


labels <- top_20_rows$gene_name
x1<-match(top_20_rows$gamma,params.sort$gamma)
y1<-top_20_rows$gamma

addTextLabels(x1,y1,tolower(labels), col.label="black", cex.label = 1.1,  col.background = NULL, lty=1)

#dev.off()


setwd('/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2023_bcftools_pnye3/gwas_jawshape/bslmm/')
# Reads input parameters
args<-commandArgs(trailingOnly=T)
#trait=args[1]
trait=c("PC1")


# library to speed up loading of big tables
library(data.table)

# Load parameters output

params<-fread(paste0("all.SNPs.filt.jaw.imp.bslmm.jaw.", trait, ".mean.params.anno.txt", sep=""),header=T,sep="\t", data.table=F)


params$pos <- params$start
params$chr <- as.numeric(gsub("chr","", params$seqnames))


# get list of linkage groups/chromosomes
chrs<-sort(as.numeric(unique(params$chr)))

# add linkage group column (chr)
chr<-gsub("lg|_.+","",params$chr)
params["chr"]<-chr


# sort by linkage group and position
params.sort<-params[order(as.numeric(params$chr), params$pos),]
params.sort$lg <- paste("LG",params.sort$chr, sep="")
params.sort$rs <- paste(params.sort$lg, params.sort$pos, sep="-")
params.sort$gene_name <- sapply(strsplit(params.sort$gene,split = ';'), `[`, 2)
params.sort$gene_name <- gsub('Name=','',params.sort$gene_name)
params.sort$eff <- abs(params.sort$beta*params.sort$gamma)

params.sort1 <- params.sort[abs(params.sort$beta) > 0,] # only variants with some effect
params.sort <- params.sort1
######################### plot of sum_gamma



# set up empty plot
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.ylim.png", sep=""), width=11.7,height=8.3,units="in",res=200)
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.png", sep=""), width=12,height=6)


# set up empty plot
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0,max(params.sort$gamma)),ylab=paste0("Top SNPs Craniofacial ", trait), xaxt="n")


# plot grey bands for chromosome/linkage groups

start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(params.sort[params.sort$chr==ch,])
  cat ("CH: ", ch, "\n")
  colour<-adjustcolor("lightgrey", alpha.f=0.1)
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=colour, border=colour)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  
  start<-start+size
}

# Add x axis labels
#axis(side=1,at=lab.pos,labels=chrs,tick=F)


# plot PIP for insignificant variants

#x<-match(params.sort$gamma[params.sort$gamma<threshold],params.sort$gamma)
top <- params.sort[quantile(params.sort$gamma,0.999)[[1]] & params.sort$gamma > 0.01,]
x<- match(params.sort$gamma[!params.sort$id %in% top$id],params.sort$gamma)

# PIP
y<-params.sort$gamma[!params.sort$id %in% top$id]

# sparse effect size, used for dot size
z<-(params.sort[!params.sort$id %in% top$id,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))
new <- data.frame(x=x,y=y,z=z)
set.seed(123)  # Set a seed for reproducibility
new1 <- new[sample(nrow(new), 20), ]
symbols(new1$x,new1$y,circles=new1$z, bg="grey90",inches=1/26,fg=NULL,add=T)


# plot threshold line
#abline(h=threshold,lty=3,col="black")

# rank of high PIP variants across linkage groups with gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
z<-(params.sort[params.sort$id %in% top$id & !is.na(params.sort$gene_name),])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#4DBBD5",alpha.f = 0.8), inches=1/26,fg="black",add=T)

# rank of high PIP variants across linkage groups without gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
#z<-(params.sort[params.sort$gamma>=threshold,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#4DBBD5",alpha.f = 0.8), inches=1/26,fg=NA,add=T)

# add label high PIP variants
#text(x,y,labels=params.sort$gene[params.sort$gamma>=threshold], adj=c(0,0), cex=0.4)
top_20_rows <- head(top[!is.na(top$gene_name), ][order(-top$gamma[!is.na(top$gene_name)]), ], 20)


labels <- top_20_rows$gene_name
x1<-match(top_20_rows$gamma,params.sort$gamma)
y1<-top_20_rows$gamma

addTextLabels(x1,y1,tolower(labels), col.label="black", cex.label = 1.1,  col.background = NULL, lty=1)

#dev.off()

#tophits <- params.sort[params.sort$gamma>=threshold,]
#write.table(tophits, paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, "99.9.PIPoutlier.SNP.txt", sep=""), quote=F, row.names=F, sep="\t")

dev.off()


########################################### SUPPMAT

## supp figures for paper


png(file="/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/manhattanplots/2025_figure3a_repeat.png", width=8,height=10,units="in",res=300)
##
par( mai=c(0.3, 0.6, 0.1, 0.1), mfrow=c(7,1))



setwd('/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_bcftools_pnye3_dentition_repeat/bslmm/') # to read in files

### multiplot 1a
# Reads input parameters
args<-commandArgs(trailingOnly=T)
#trait=args[1]
trait=c("CuspShape")


# library to speed up loading of big tables
library(data.table)

# Load parameters output

params<-fread(paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.", trait, ".mean.params.anno.txt", sep=""),header=T,sep="\t", data.table=F)


params$pos <- params$start
params$chr <- as.numeric(gsub("chr","", params$seqnames))


# get list of linkage groups/chromosomes
chrs<-sort(as.numeric(unique(params$chr)))

# add linkage group column (chr)
chr<-gsub("lg|_.+","",params$chr)
params["chr"]<-chr


# sort by linkage group and position
params.sort<-params[order(as.numeric(params$chr), params$pos),]
params.sort$lg <- paste("LG",params.sort$chr, sep="")
params.sort$rs <- paste(params.sort$lg, params.sort$pos, sep="-")
params.sort$gene_name <- sapply(strsplit(params.sort$gene,split = ';'), `[`, 2)
params.sort$gene_name <- gsub('Name=','',params.sort$gene_name)
params.sort$eff <- abs(params.sort$beta*params.sort$gamma)

params.sort1 <- params.sort[abs(params.sort$beta) > 0,] # only variants with some effect
params.sort <- params.sort1
######################### plot of sum_gamma



# set up empty plot
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.ylim.png", sep=""), width=11.7,height=8.3,units="in",res=200)
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.png", sep=""), width=12,height=6)


# set up empty plot
#plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0,max(params.sort$gamma)),ylab=paste0("Top SNPs ", trait, " (PIP)"), xaxt="n")
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0,1),ylab=paste0("Top SNPs ", trait, " (PIP)"), xaxt="n")


# plot grey bands for chromosome/linkage groups

start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(params.sort[params.sort$chr==ch,])
  cat ("CH: ", ch, "\n")
  colour<-adjustcolor("lightgrey", alpha.f=0.1)
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=colour, border=colour)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  
  start<-start+size
}

# Add x axis labels
#axis(side=1,at=lab.pos,labels=chrs,tick=F)


# highligh high PIP variants

#threshold <- 0.01
#threshold <- 0.01 
#abline(h=threshold,lty=3,col="black")

# plot PIP for insignificant variants

#x<-match(params.sort$gamma[params.sort$gamma<threshold],params.sort$gamma)
top <- params.sort[quantile(params.sort$gamma,0.9999)[[1]] & params.sort$gamma > 0.01,]
x<- match(params.sort$gamma[!params.sort$id %in% top$id],params.sort$gamma)

# PIP
y<-params.sort$gamma[!params.sort$id %in% top$id]

# sparse effect size, used for dot size
z<-(params.sort[!params.sort$id %in% top$id,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))
symbols(x,y,circles=z, bg="grey90",inches=1/26,fg=NULL,add=T)

# plot threshold line
#abline(h=threshold,lty=3,col="black")

# rank of high PIP variants across linkage groups with gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
z<-(params.sort[params.sort$id %in% top$id & !is.na(params.sort$gene_name),])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#4DBBD5",alpha.f = 0.8), inches=1/22,fg="black",add=T)

# rank of high PIP variants across linkage groups without gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
#z<-(params.sort[params.sort$gamma>=threshold,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#4DBBD5",alpha.f = 0.8), inches=1/22,fg=NA,add=T)

# add label high PIP variants
#text(x,y,labels=params.sort$gene[params.sort$gamma>=threshold], adj=c(0,0), cex=0.4)
top_20_rows <- head(top[!is.na(top$gene_name), ][order(-top$gamma[!is.na(top$gene_name)]), ], 20)


labels <- top_20_rows$gene_name
x1<-match(top_20_rows$gamma,params.sort$gamma)
y1<-top_20_rows$gamma

addTextLabels(x1,y1,tolower(labels), col.label="black", cex.label = 1.1,  col.background = NULL, lty=1)

#dev.off()


#tophits <- params.sort[params.sort$gamma>=threshold,]
#write.table(tophits, paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, "99.9.PIPoutlier.SNP.txt", sep=""), quote=F, row.names=F, sep="\t")


### multiplot 1b


# Reads input parameters
args<-commandArgs(trailingOnly=T)
#trait=args[1]
trait=c("InnertoothrowUJ")


# library to speed up loading of big tables
library(data.table)

# Load parameters output

params<-fread(paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.", trait, ".mean.params.anno.txt", sep=""),header=T,sep="\t", data.table=F)


params$pos <- params$start
params$chr <- as.numeric(gsub("chr","", params$seqnames))


# get list of linkage groups/chromosomes
chrs<-sort(as.numeric(unique(params$chr)))

# add linkage group column (chr)
chr<-gsub("lg|_.+","",params$chr)
params["chr"]<-chr


# sort by linkage group and position
params.sort<-params[order(as.numeric(params$chr), params$pos),]
params.sort$lg <- paste("LG",params.sort$chr, sep="")
params.sort$rs <- paste(params.sort$lg, params.sort$pos, sep="-")
params.sort$gene_name <- sapply(strsplit(params.sort$gene,split = ';'), `[`, 2)
params.sort$gene_name <- gsub('Name=','',params.sort$gene_name)
params.sort$eff <- abs(params.sort$beta*params.sort$gamma)

params.sort1 <- params.sort[abs(params.sort$beta) > 0,] # only variants with some effect
params.sort <- params.sort1
######################### plot of sum_gamma



# set up empty plot
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.ylim.png", sep=""), width=11.7,height=8.3,units="in",res=200)
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.png", sep=""), width=12,height=6)


# set up empty plot
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0,max(params.sort$gamma)),ylab=paste0("Top SNPs ", trait, " (PIP)"), xaxt="n")


# plot grey bands for chromosome/linkage groups

start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(params.sort[params.sort$chr==ch,])
  cat ("CH: ", ch, "\n")
  colour<-adjustcolor("lightgrey", alpha.f=0.1)
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=colour, border=colour)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  
  start<-start+size
}

# Add x axis labels
#axis(side=1,at=lab.pos,labels=chrs,tick=F)

# plot PIP for insignificant variants

#x<-match(params.sort$gamma[params.sort$gamma<threshold],params.sort$gamma)
top <- params.sort[quantile(params.sort$gamma,0.9999)[[1]] & params.sort$gamma > 0.01,]
x<- match(params.sort$gamma[!params.sort$id %in% top$id],params.sort$gamma)

# PIP
y<-params.sort$gamma[!params.sort$id %in% top$id]

# sparse effect size, used for dot size
z<-(params.sort[!params.sort$id %in% top$id,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))
symbols(x,y,circles=z, bg="grey90",inches=1/26,fg=NULL,add=T)

# plot threshold line
#abline(h=threshold,lty=3,col="black")

# rank of high PIP variants across linkage groups with gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
z<-(params.sort[params.sort$id %in% top$id & !is.na(params.sort$gene_name),])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#4DBBD5",alpha.f = 0.8), inches=1/26,fg="black",add=T)

# rank of high PIP variants across linkage groups without gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
#z<-(params.sort[params.sort$gamma>=threshold,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#4DBBD5",alpha.f = 0.8), inches=1/26,fg=NA,add=T)

# add label high PIP variants
#text(x,y,labels=params.sort$gene[params.sort$gamma>=threshold], adj=c(0,0), cex=0.4)
top_20_rows <- head(top[!is.na(top$gene_name), ][order(-top$gamma[!is.na(top$gene_name)]), ], 20)


labels <- top_20_rows$gene_name
x1<-match(top_20_rows$gamma,params.sort$gamma)
y1<-top_20_rows$gamma

addTextLabels(x1,y1,tolower(labels), col.label="black", cex.label = 1.1,  col.background = NULL, lty=1)

#dev.off()

#tophits <- params.sort[params.sort$gamma>=threshold,]
#write.table(tophits, paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, "99.9.PIPoutlier.SNP.txt", sep=""), quote=F, row.names=F, sep="\t")


###### multiplot begins ### multiplot 1c


# Reads input parameters
args<-commandArgs(trailingOnly=T)
#trait=args[1]
trait=c("LPJkeelDepth")


# library to speed up loading of big tables
library(data.table)

# Load parameters output

params<-fread(paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.", trait, ".mean.params.anno.txt", sep=""),header=T,sep="\t", data.table=F)


params$pos <- params$start
params$chr <- as.numeric(gsub("chr","", params$seqnames))


# get list of linkage groups/chromosomes
chrs<-sort(as.numeric(unique(params$chr)))

# add linkage group column (chr)
chr<-gsub("lg|_.+","",params$chr)
params["chr"]<-chr


# sort by linkage group and position
params.sort<-params[order(as.numeric(params$chr), params$pos),]
params.sort$lg <- paste("LG",params.sort$chr, sep="")
params.sort$rs <- paste(params.sort$lg, params.sort$pos, sep="-")
params.sort$gene_name <- sapply(strsplit(params.sort$gene,split = ';'), `[`, 2)
params.sort$gene_name <- gsub('Name=','',params.sort$gene_name)
params.sort$eff <- abs(params.sort$beta*params.sort$gamma)

params.sort1 <- params.sort[abs(params.sort$beta) > 0,] # only variants with some effect
params.sort <- params.sort1
######################### plot of sum_gamma



# set up empty plot
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.ylim.png", sep=""), width=11.7,height=8.3,units="in",res=200)
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.png", sep=""), width=12,height=6)


# set up empty plot
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0,max(params.sort$gamma)),ylab=paste0("Top SNPs ", trait, " (PIP)"), xaxt="n")


# plot grey bands for chromosome/linkage groups

start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(params.sort[params.sort$chr==ch,])
  cat ("CH: ", ch, "\n")
  colour<-adjustcolor("lightgrey", alpha.f=0.1)
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=colour, border=colour)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  
  start<-start+size
}

# Add x axis labels
#axis(side=1,at=lab.pos,labels=chrs,tick=F)


# plot PIP for insignificant variants

#x<-match(params.sort$gamma[params.sort$gamma<threshold],params.sort$gamma)
top <- params.sort[quantile(params.sort$gamma,0.9999)[[1]] & params.sort$gamma > 0.01,]
x<- match(params.sort$gamma[!params.sort$id %in% top$id],params.sort$gamma)

# PIP
y<-params.sort$gamma[!params.sort$id %in% top$id]

# sparse effect size, used for dot size
z<-(params.sort[!params.sort$id %in% top$id,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))
symbols(x,y,circles=z, bg="grey90",inches=1/26,fg=NULL,add=T)

# plot threshold line
#abline(h=threshold,lty=3,col="black")

# rank of high PIP variants across linkage groups with gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
z<-(params.sort[params.sort$id %in% top$id & !is.na(params.sort$gene_name),])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#4DBBD5",alpha.f = 0.8), inches=1/26,fg="black",add=T)

# rank of high PIP variants across linkage groups without gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
#z<-(params.sort[params.sort$gamma>=threshold,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#4DBBD5",alpha.f = 0.8), inches=1/26,fg=NA,add=T)

# add label high PIP variants
#text(x,y,labels=params.sort$gene[params.sort$gamma>=threshold], adj=c(0,0), cex=0.4)
top_20_rows <- head(top[!is.na(top$gene_name), ][order(-top$gamma[!is.na(top$gene_name)]), ], 20)


labels <- top_20_rows$gene_name
x1<-match(top_20_rows$gamma,params.sort$gamma)
y1<-top_20_rows$gamma

addTextLabels(x1,y1,tolower(labels), col.label="black", cex.label = 1.1,  col.background = NULL, lty=1)

#dev.off()

#tophits <- params.sort[params.sort$gamma>=threshold,]
#write.table(tophits, paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, "99.9.PIPoutlier.SNP.txt", sep=""), quote=F, row.names=F, sep="\t")


### multiplot 1d YF


# Reads input parameters
args<-commandArgs(trailingOnly=T)
#trait=args[1]
trait=c("LPJtoothMolarisation")


# library to speed up loading of big tables
library(data.table)

# Load parameters output

params<-fread(paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.", trait, ".mean.params.anno.txt", sep=""),header=T,sep="\t", data.table=F)


params$pos <- params$start
params$chr <- as.numeric(gsub("chr","", params$seqnames))


# get list of linkage groups/chromosomes
chrs<-sort(as.numeric(unique(params$chr)))

# add linkage group column (chr)
chr<-gsub("lg|_.+","",params$chr)
params["chr"]<-chr


# sort by linkage group and position
params.sort<-params[order(as.numeric(params$chr), params$pos),]
params.sort$lg <- paste("LG",params.sort$chr, sep="")
params.sort$rs <- paste(params.sort$lg, params.sort$pos, sep="-")
params.sort$gene_name <- sapply(strsplit(params.sort$gene,split = ';'), `[`, 2)
params.sort$gene_name <- gsub('Name=','',params.sort$gene_name)
params.sort$eff <- abs(params.sort$beta*params.sort$gamma)

params.sort1 <- params.sort[abs(params.sort$beta) > 0,] # only variants with some effect
params.sort <- params.sort1
######################### plot of sum_gamma



# set up empty plot
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.ylim.png", sep=""), width=11.7,height=8.3,units="in",res=200)
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.png", sep=""), width=12,height=6)


# set up empty plot
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0,max(params.sort$gamma)),ylab=paste0("Top SNPs ", trait, " (PIP)"), xaxt="n")


# plot grey bands for chromosome/linkage groups

start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(params.sort[params.sort$chr==ch,])
  cat ("CH: ", ch, "\n")
  colour<-adjustcolor("lightgrey", alpha.f=0.1)
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=colour, border=colour)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  
  start<-start+size
}

# Add x axis labels
#axis(side=1,at=lab.pos,labels=chrs,tick=F)

# plot PIP for insignificant variants

#x<-match(params.sort$gamma[params.sort$gamma<threshold],params.sort$gamma)
top <- params.sort[quantile(params.sort$gamma,0.9999)[[1]] & params.sort$gamma > 0.01,]
x<- match(params.sort$gamma[!params.sort$id %in% top$id],params.sort$gamma)

# PIP
y<-params.sort$gamma[!params.sort$id %in% top$id]

# sparse effect size, used for dot size
z<-(params.sort[!params.sort$id %in% top$id,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))
symbols(x,y,circles=z, bg="grey90",inches=1/26,fg=NULL,add=T)

# plot threshold line
#abline(h=threshold,lty=3,col="black")

# rank of high PIP variants across linkage groups with gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
z<-(params.sort[params.sort$id %in% top$id & !is.na(params.sort$gene_name),])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#4DBBD5",alpha.f = 0.8), inches=1/26,fg="black",add=T)

# rank of high PIP variants across linkage groups without gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
#z<-(params.sort[params.sort$gamma>=threshold,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#4DBBD5",alpha.f = 0.8), inches=1/26,fg=NA,add=T)

# add label high PIP variants
#text(x,y,labels=params.sort$gene[params.sort$gamma>=threshold], adj=c(0,0), cex=0.4)
top_20_rows <- head(top[!is.na(top$gene_name), ][order(-top$gamma[!is.na(top$gene_name)]), ], 20)


labels <- top_20_rows$gene_name
x1<-match(top_20_rows$gamma,params.sort$gamma)
y1<-top_20_rows$gamma

addTextLabels(x1,y1,tolower(labels), col.label="black", cex.label = 1.1,  col.background = NULL, lty=1)

#tophits <- params.sort[params.sort$gamma>=threshold,]
#write.table(tophits, paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, "99.9.PIPoutlier.SNP.txt", sep=""), quote=F, row.names=F, sep="\t")



###### multiplot begins 1e

setwd('/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2023_bcftools_pnye3/gwas_jawshape/bslmm/')
# Reads input parameters
args<-commandArgs(trailingOnly=T)
#trait=args[1]
trait=c("PC1")


# library to speed up loading of big tables
library(data.table)

# Load parameters output

params<-fread(paste0("all.SNPs.filt.jaw.imp.bslmm.jaw.", trait, ".mean.params.anno.txt", sep=""),header=T,sep="\t", data.table=F)


params$pos <- params$start
params$chr <- as.numeric(gsub("chr","", params$seqnames))


# get list of linkage groups/chromosomes
chrs<-sort(as.numeric(unique(params$chr)))

# add linkage group column (chr)
chr<-gsub("lg|_.+","",params$chr)
params["chr"]<-chr


# sort by linkage group and position
params.sort<-params[order(as.numeric(params$chr), params$pos),]
params.sort$lg <- paste("LG",params.sort$chr, sep="")
params.sort$rs <- paste(params.sort$lg, params.sort$pos, sep="-")
params.sort$gene_name <- sapply(strsplit(params.sort$gene,split = ';'), `[`, 2)
params.sort$gene_name <- gsub('Name=','',params.sort$gene_name)
params.sort$eff <- abs(params.sort$beta*params.sort$gamma)

params.sort1 <- params.sort[abs(params.sort$beta) > 0,] # only variants with some effect
params.sort <- params.sort1
######################### plot of sum_gamma



# set up empty plot
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.ylim.png", sep=""), width=11.7,height=8.3,units="in",res=200)
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.png", sep=""), width=12,height=6)


# set up empty plot
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0,max(params.sort$gamma)),ylab=paste0("Top SNPs Craniofacial ", trait), xaxt="n")


# plot grey bands for chromosome/linkage groups

start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(params.sort[params.sort$chr==ch,])
  cat ("CH: ", ch, "\n")
  colour<-adjustcolor("lightgrey", alpha.f=0.1)
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=colour, border=colour)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  
  start<-start+size
}

# Add x axis labels
#axis(side=1,at=lab.pos,labels=chrs,tick=F)


# plot PIP for insignificant variants

#x<-match(params.sort$gamma[params.sort$gamma<threshold],params.sort$gamma)
top <- params.sort[quantile(params.sort$gamma,0.999)[[1]] & params.sort$gamma > 0.01,]
x<- match(params.sort$gamma[!params.sort$id %in% top$id],params.sort$gamma)

# PIP
y<-params.sort$gamma[!params.sort$id %in% top$id]

# sparse effect size, used for dot size
z<-(params.sort[!params.sort$id %in% top$id,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))
symbols(x,y,circles=z, bg="grey90",inches=1/26,fg=NULL,add=T)

# plot threshold line
#abline(h=threshold,lty=3,col="black")

# rank of high PIP variants across linkage groups with gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
z<-(params.sort[params.sort$id %in% top$id & !is.na(params.sort$gene_name),])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#4DBBD5",alpha.f = 0.8), inches=1/26,fg="black",add=T)

# rank of high PIP variants across linkage groups without gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
#z<-(params.sort[params.sort$gamma>=threshold,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#4DBBD5",alpha.f = 0.8), inches=1/26,fg=NA,add=T)

# add label high PIP variants
#text(x,y,labels=params.sort$gene[params.sort$gamma>=threshold], adj=c(0,0), cex=0.4)
top_20_rows <- head(top[!is.na(top$gene_name), ][order(-top$gamma[!is.na(top$gene_name)]), ], 20)


labels <- top_20_rows$gene_name
x1<-match(top_20_rows$gamma,params.sort$gamma)
y1<-top_20_rows$gamma

addTextLabels(x1,y1,tolower(labels), col.label="black", cex.label = 1.1,  col.background = NULL, lty=1)

#dev.off()

#tophits <- params.sort[params.sort$gamma>=threshold,]
#write.table(tophits, paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, "99.9.PIPoutlier.SNP.txt", sep=""), quote=F, row.names=F, sep="\t")

#dev.off()

###### multiplot begins 1f

# Reads input parameters
args<-commandArgs(trailingOnly=T)
#trait=args[1]
trait=c("PC2")


# library to speed up loading of big tables
library(data.table)

# Load parameters output

params<-fread(paste0("all.SNPs.filt.jaw.imp.bslmm.jaw.", trait, ".mean.params.anno.txt", sep=""),header=T,sep="\t", data.table=F)


params$pos <- params$start
params$chr <- as.numeric(gsub("chr","", params$seqnames))


# get list of linkage groups/chromosomes
chrs<-sort(as.numeric(unique(params$chr)))

# add linkage group column (chr)
chr<-gsub("lg|_.+","",params$chr)
params["chr"]<-chr


# sort by linkage group and position
params.sort<-params[order(as.numeric(params$chr), params$pos),]
params.sort$lg <- paste("LG",params.sort$chr, sep="")
params.sort$rs <- paste(params.sort$lg, params.sort$pos, sep="-")
params.sort$gene_name <- sapply(strsplit(params.sort$gene,split = ';'), `[`, 2)
params.sort$gene_name <- gsub('Name=','',params.sort$gene_name)
params.sort$eff <- abs(params.sort$beta*params.sort$gamma)

params.sort1 <- params.sort[abs(params.sort$beta) > 0,] # only variants with some effect
params.sort <- params.sort1
######################### plot of sum_gamma



# set up empty plot
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.ylim.png", sep=""), width=11.7,height=8.3,units="in",res=200)
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.png", sep=""), width=12,height=6)


# set up empty plot
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0,max(params.sort$gamma)),ylab=paste0("Top SNPs Craniofacial ", trait), xaxt="n")


# plot grey bands for chromosome/linkage groups

start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(params.sort[params.sort$chr==ch,])
  cat ("CH: ", ch, "\n")
  colour<-adjustcolor("lightgrey", alpha.f=0.1)
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=colour, border=colour)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  
  start<-start+size
}

# Add x axis labels
#axis(side=1,at=lab.pos,labels=chrs,tick=F)

# plot PIP for insignificant variants

#x<-match(params.sort$gamma[params.sort$gamma<threshold],params.sort$gamma)
top <- params.sort[quantile(params.sort$gamma,0.999)[[1]] & params.sort$gamma > 0.01,]
x<- match(params.sort$gamma[!params.sort$id %in% top$id],params.sort$gamma)

# PIP
y<-params.sort$gamma[!params.sort$id %in% top$id]

# sparse effect size, used for dot size
z<-(params.sort[!params.sort$id %in% top$id,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))
symbols(x,y,circles=z, bg="grey90",inches=1/26,fg=NULL,add=T)

# plot threshold line
#abline(h=threshold,lty=3,col="black")

# rank of high PIP variants across linkage groups with gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
z<-(params.sort[params.sort$id %in% top$id & !is.na(params.sort$gene_name),])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#4DBBD5",alpha.f = 0.8), inches=1/26,fg="black",add=T)

# rank of high PIP variants across linkage groups without gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
#z<-(params.sort[params.sort$gamma>=threshold,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#4DBBD5",alpha.f = 0.8), inches=1/26,fg=NA,add=T)

# add label high PIP variants
#text(x,y,labels=params.sort$gene[params.sort$gamma>=threshold], adj=c(0,0), cex=0.4)
top_20_rows <- head(top[!is.na(top$gene_name), ][order(-top$gamma[!is.na(top$gene_name)]), ], 20)


labels <- top_20_rows$gene_name
x1<-match(top_20_rows$gamma,params.sort$gamma)
y1<-top_20_rows$gamma

addTextLabels(x1,y1,tolower(labels), col.label="black", cex.label = 1.1,  col.background = NULL, lty=1)


#tophits <- params.sort[params.sort$gamma>=threshold,]
#write.table(tophits, paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, "99.9.PIPoutlier.SNP.txt", sep=""), quote=F, row.names=F, sep="\t")



###### multiplot begins 1g
setwd('/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_repeatcolour/bslmm/')
# Reads input parameters
args<-commandArgs(trailingOnly=T)
#trait=args[1]
trait=c("VerticalBars")


# library to speed up loading of big tables
library(data.table)

# Load parameters output

params<-fread(paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".mean.params.anno.txt", sep=""),header=T,sep="\t", data.table=F)


params$pos <- params$start
params$chr <- as.numeric(gsub("chr","", params$seqnames))


# get list of linkage groups/chromosomes
chrs<-sort(as.numeric(unique(params$chr)))

# add linkage group column (chr)
chr<-gsub("lg|_.+","",params$chr)
params["chr"]<-chr


# sort by linkage group and position
params.sort<-params[order(as.numeric(params$chr), params$pos),]
params.sort$lg <- paste("LG",params.sort$chr, sep="")
params.sort$rs <- paste(params.sort$lg, params.sort$pos, sep="-")
params.sort$gene_name <- sapply(strsplit(params.sort$gene,split = ';'), `[`, 2)
params.sort$gene_name <- gsub('Name=','',params.sort$gene_name)
params.sort$eff <- abs(params.sort$beta*params.sort$gamma)

params.sort1 <- params.sort[abs(params.sort$beta) > 0,] # only variants with some effect
params.sort <- params.sort1
######################### plot of sum_gamma



# set up empty plot
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.ylim.png", sep=""), width=11.7,height=8.3,units="in",res=200)
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.png", sep=""), width=12,height=6)


# set up empty plot
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0,max(params.sort$gamma)),ylab=paste0("Top SNPs ", trait, " (PIP)"), xaxt="n")


# plot grey bands for chromosome/linkage groups

start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(params.sort[params.sort$chr==ch,])
  cat ("CH: ", ch, "\n")
  colour<-adjustcolor("lightgrey", alpha.f=0.1)
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=colour, border=colour)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  
  start<-start+size
}

# Add x axis labels
axis(side=1,at=lab.pos,labels=chrs,tick=F)


# plot PIP for insignificant variants

#x<-match(params.sort$gamma[params.sort$gamma<threshold],params.sort$gamma)
top <- params.sort[quantile(params.sort$gamma,0.9999)[[1]] & params.sort$gamma > 0.01,]
x<- match(params.sort$gamma[!params.sort$id %in% top$id],params.sort$gamma)

# PIP
y<-params.sort$gamma[!params.sort$id %in% top$id]

# sparse effect size, used for dot size
z<-(params.sort[!params.sort$id %in% top$id,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))
symbols(x,y,circles=z, bg="grey90",inches=1/26,fg=NULL,add=T)

# plot threshold line
#abline(h=threshold,lty=3,col="black")

# rank of high PIP variants across linkage groups with gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
z<-(params.sort[params.sort$id %in% top$id & !is.na(params.sort$gene_name),])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#009E73",alpha.f = 0.8), inches=1/26,fg="black",add=T)

# rank of high PIP variants across linkage groups without gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
#z<-(params.sort[params.sort$gamma>=threshold,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#009E73",alpha.f = 0.8), inches=1/26,fg=NA,add=T)

# add label high PIP variants
#text(x,y,labels=params.sort$gene[params.sort$gamma>=threshold], adj=c(0,0), cex=0.4)
top_20_rows <- head(top[!is.na(top$gene_name), ][order(-top$gamma[!is.na(top$gene_name)]), ], 20)


labels <- top_20_rows$gene_name
x1<-match(top_20_rows$gamma,params.sort$gamma)
y1<-top_20_rows$gamma

addTextLabels(x1,y1,tolower(labels), col.label="black", cex.label = 1.1,  col.background = NULL, lty=1)

#tophits <- params.sort[params.sort$gamma>=threshold,]
#write.table(tophits, paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, "99.9.PIPoutlier.SNP.txt", sep=""), quote=F, row.names=F, sep="\t")

dev.off()


############################## multiplot begins 2b 


png(file="/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/manhattanplots/2025_figure3b_repeat.png", width=8,height=10,units="in",res=300)
##
par( mai=c(0.3, 0.6, 0.1, 0.1), mfrow=c(7,1))
setwd('/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_repeatcolour/bslmm')

options(echo=TRUE)
library(zoo) 
library(basicPlotteR)



#### 2a
# Reads input parameters
args<-commandArgs(trailingOnly=T)
#trait=args[1]
trait=c("MidLateralBand")

setwd('/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_repeatcolour/bslmm/')

# library to speed up loading of big tables
library(data.table)

# Load parameters output

params<-fread(paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".mean.params.anno.txt", sep=""),header=T,sep="\t", data.table=F)


params$pos <- params$start
params$chr <- as.numeric(gsub("chr","", params$seqnames))


# get list of linkage groups/chromosomes
chrs<-sort(as.numeric(unique(params$chr)))

# add linkage group column (chr)
chr<-gsub("lg|_.+","",params$chr)
params["chr"]<-chr


# sort by linkage group and position
params.sort<-params[order(as.numeric(params$chr), params$pos),]
params.sort$lg <- paste("LG",params.sort$chr, sep="")
params.sort$rs <- paste(params.sort$lg, params.sort$pos, sep="-")
params.sort$gene_name <- sapply(strsplit(params.sort$gene,split = ';'), `[`, 2)
params.sort$gene_name <- gsub('Name=','',params.sort$gene_name)
params.sort$eff <- abs(params.sort$beta*params.sort$gamma)


params.sort1 <- params.sort[abs(params.sort$beta) > 0,] # only variants with some effect
params.sort <- params.sort1
######################### plot of sum_gamma



# set up empty plot
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.ylim.png", sep=""), width=11.7,height=8.3,units="in",res=200)
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.png", sep=""), width=12,height=6)


# set up empty plot
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0,max(params.sort$gamma)),ylab=paste0("Top SNPs ", trait, " (PIP)"),xlab="linkage group", xaxt="n")


# plot grey bands for chromosome/linkage groups

start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(params.sort[params.sort$chr==ch,])
  cat ("CH: ", ch, "\n")
  colour<-adjustcolor("lightgrey", alpha.f=0.1)
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=colour, border=colour)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  
  start<-start+size
}

# Add x axis labels
axis(side=1,at=lab.pos,labels=chrs,tick=F)


# plot PIP for insignificant variants

#x<-match(params.sort$gamma[params.sort$gamma<threshold],params.sort$gamma)
top <- params.sort[quantile(params.sort$gamma,0.9999)[[1]] & params.sort$gamma > 0.01,]
x<- match(params.sort$gamma[!params.sort$id %in% top$id],params.sort$gamma)

# PIP
y<-params.sort$gamma[!params.sort$id %in% top$id]

# sparse effect size, used for dot size
z<-(params.sort[!params.sort$id %in% top$id,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))
symbols(x,y,circles=z, bg="grey90",inches=1/26,fg=NULL,add=T)

# plot threshold line
#abline(h=threshold,lty=3,col="black")

# rank of high PIP variants across linkage groups with gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
z<-(params.sort[params.sort$id %in% top$id & !is.na(params.sort$gene_name),])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#009E73",alpha.f = 0.8), inches=1/26,fg="black",add=T)

# rank of high PIP variants across linkage groups without gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
#z<-(params.sort[params.sort$gamma>=threshold,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#009E73",alpha.f = 0.8), inches=1/26,fg=NA,add=T)

# add label high PIP variants
#text(x,y,labels=params.sort$gene[params.sort$gamma>=threshold], adj=c(0,0), cex=0.4)
top_20_rows <- head(top[!is.na(top$gene_name), ][order(-top$gamma[!is.na(top$gene_name)]), ], 20)


labels <- top_20_rows$gene_name
x1<-match(top_20_rows$gamma,params.sort$gamma)
y1<-top_20_rows$gamma

addTextLabels(x1,y1,tolower(labels), col.label="black", cex.label = 1.1,  col.background = NULL, lty=1)

#dev.off()

#tophits <- params.sort[params.sort$gamma>=threshold,]
#write.table(tophits, paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, "99.9.PIPoutlier.SNP.txt", sep=""), quote=F, row.names=F, sep="\t")


#### 2b


# Reads input parameters
args<-commandArgs(trailingOnly=T)
#trait=args[1]
trait=c("DorsalLateralBand")

# library to speed up loading of big tables
library(data.table)

# Load parameters output

params<-fread(paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".mean.params.anno.txt", sep=""),header=T,sep="\t", data.table=F)


params$pos <- params$start
params$chr <- as.numeric(gsub("chr","", params$seqnames))


# get list of linkage groups/chromosomes
chrs<-sort(as.numeric(unique(params$chr)))

# add linkage group column (chr)
chr<-gsub("lg|_.+","",params$chr)
params["chr"]<-chr


# sort by linkage group and position
params.sort<-params[order(as.numeric(params$chr), params$pos),]
params.sort$lg <- paste("LG",params.sort$chr, sep="")
params.sort$rs <- paste(params.sort$lg, params.sort$pos, sep="-")
params.sort$gene_name <- sapply(strsplit(params.sort$gene,split = ';'), `[`, 2)
params.sort$gene_name <- gsub('Name=','',params.sort$gene_name)
params.sort$eff <- abs(params.sort$beta*params.sort$gamma)

params.sort1 <- params.sort[abs(params.sort$beta) > 0,] # only variants with some effect
params.sort <- params.sort1
######################### plot of sum_gamma



# set up empty plot
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.ylim.png", sep=""), width=11.7,height=8.3,units="in",res=200)
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.png", sep=""), width=12,height=6)


# set up empty plot
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0,max(params.sort$gamma)),ylab=paste0("Top SNPs ", trait, " (PIP)"),xlab="linkage group", xaxt="n")


# plot grey bands for chromosome/linkage groups

start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(params.sort[params.sort$chr==ch,])
  cat ("CH: ", ch, "\n")
  colour<-adjustcolor("lightgrey", alpha.f=0.1)
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=colour, border=colour)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  
  start<-start+size
}

# Add x axis labels
axis(side=1,at=lab.pos,labels=chrs,tick=F)


# plot PIP for insignificant variants

#x<-match(params.sort$gamma[params.sort$gamma<threshold],params.sort$gamma)
top <- params.sort[quantile(params.sort$gamma,0.9999)[[1]] & params.sort$gamma > 0.01,]
x<- match(params.sort$gamma[!params.sort$id %in% top$id],params.sort$gamma)

# PIP
y<-params.sort$gamma[!params.sort$id %in% top$id]

# sparse effect size, used for dot size
z<-(params.sort[!params.sort$id %in% top$id,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))
symbols(x,y,circles=z, bg="grey90",inches=1/26,fg=NULL,add=T)

# plot threshold line
#abline(h=threshold,lty=3,col="black")

# rank of high PIP variants across linkage groups with gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
z<-(params.sort[params.sort$id %in% top$id & !is.na(params.sort$gene_name),])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#009E73",alpha.f = 0.8), inches=1/26,fg="black",add=T)

# rank of high PIP variants across linkage groups without gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
#z<-(params.sort[params.sort$gamma>=threshold,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#009E73",alpha.f = 0.8), inches=1/26,fg=NA,add=T)

# add label high PIP variants
#text(x,y,labels=params.sort$gene[params.sort$gamma>=threshold], adj=c(0,0), cex=0.4)
top_20_rows <- head(top[!is.na(top$gene_name), ][order(-top$gamma[!is.na(top$gene_name)]), ], 20)


labels <- top_20_rows$gene_name
x1<-match(top_20_rows$gamma,params.sort$gamma)
y1<-top_20_rows$gamma

addTextLabels(x1,y1,tolower(labels), col.label="black", cex.label = 1.1,  col.background = NULL, lty=1)


#dev.off()

#tophits <- params.sort[params.sort$gamma>=threshold,]
#write.table(tophits, paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, "99.9.PIPoutlier.SNP.txt", sep=""), quote=F, row.names=F, sep="\t")

#### 2c


# Reads input parameters
args<-commandArgs(trailingOnly=T)
#trait=args[1]
trait=c("Flameback")

# library to speed up loading of big tables
library(data.table)

# Load parameters output

params<-fread(paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".mean.params.anno.txt", sep=""),header=T,sep="\t", data.table=F)


params$pos <- params$start
params$chr <- as.numeric(gsub("chr","", params$seqnames))


# get list of linkage groups/chromosomes
chrs<-sort(as.numeric(unique(params$chr)))

# add linkage group column (chr)
chr<-gsub("lg|_.+","",params$chr)
params["chr"]<-chr


# sort by linkage group and position
params.sort<-params[order(as.numeric(params$chr), params$pos),]
params.sort$lg <- paste("LG",params.sort$chr, sep="")
params.sort$rs <- paste(params.sort$lg, params.sort$pos, sep="-")
params.sort$gene_name <- sapply(strsplit(params.sort$gene,split = ';'), `[`, 2)
params.sort$gene_name <- gsub('Name=','',params.sort$gene_name)
params.sort$eff <- abs(params.sort$beta*params.sort$gamma)

params.sort1 <- params.sort[abs(params.sort$beta) > 0,] # only variants with some effect
params.sort <- params.sort1
######################### plot of sum_gamma



# set up empty plot
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.ylim.png", sep=""), width=11.7,height=8.3,units="in",res=200)
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.png", sep=""), width=12,height=6)


# set up empty plot
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0,max(params.sort$gamma)),ylab=paste0("Top SNPs ", trait, " (PIP)"),xlab="linkage group", xaxt="n")


# plot grey bands for chromosome/linkage groups

start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(params.sort[params.sort$chr==ch,])
  cat ("CH: ", ch, "\n")
  colour<-adjustcolor("lightgrey", alpha.f=0.1)
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=colour, border=colour)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  
  start<-start+size
}

# Add x axis labels
axis(side=1,at=lab.pos,labels=chrs,tick=F)


# plot PIP for insignificant variants

#x<-match(params.sort$gamma[params.sort$gamma<threshold],params.sort$gamma)
top <- params.sort[quantile(params.sort$gamma,0.9999)[[1]] & params.sort$gamma > 0.01,]
x<- match(params.sort$gamma[!params.sort$id %in% top$id],params.sort$gamma)

# PIP
y<-params.sort$gamma[!params.sort$id %in% top$id]

# sparse effect size, used for dot size
z<-(params.sort[!params.sort$id %in% top$id,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))
symbols(x,y,circles=z, bg="grey90",inches=1/26,fg=NULL,add=T)

# plot threshold line
#abline(h=threshold,lty=3,col="black")

# rank of high PIP variants across linkage groups with gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
z<-(params.sort[params.sort$id %in% top$id & !is.na(params.sort$gene_name),])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#E64B35",alpha.f = 0.8), inches=1/26,fg="black",add=T)

# rank of high PIP variants across linkage groups without gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
#z<-(params.sort[params.sort$gamma>=threshold,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#E64B35",alpha.f = 0.8), inches=1/26,fg=NA,add=T)

# add label high PIP variants
#text(x,y,labels=params.sort$gene[params.sort$gamma>=threshold], adj=c(0,0), cex=0.4)
top_20_rows <- head(top[!is.na(top$gene_name), ][order(-top$gamma[!is.na(top$gene_name)]), ], 20)


labels <- top_20_rows$gene_name
x1<-match(top_20_rows$gamma,params.sort$gamma)
y1<-top_20_rows$gamma

addTextLabels(x1,y1,tolower(labels), col.label="black", cex.label = 1.1,  col.background = NULL, lty=1)


#dev.off()

#tophits <- params.sort[params.sort$gamma>=threshold,]
#write.table(tophits, paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, "99.9.PIPoutlier.SNP.txt", sep=""), quote=F, row.names=F, sep="\t")

#### 2d


# Reads input parameters
args<-commandArgs(trailingOnly=T)
#trait=args[1]
trait=c("YellowFlank")

# library to speed up loading of big tables
library(data.table)

# Load parameters output

params<-fread(paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".mean.params.anno.txt", sep=""),header=T,sep="\t", data.table=F)


params$pos <- params$start
params$chr <- as.numeric(gsub("chr","", params$seqnames))


# get list of linkage groups/chromosomes
chrs<-sort(as.numeric(unique(params$chr)))

# add linkage group column (chr)
chr<-gsub("lg|_.+","",params$chr)
params["chr"]<-chr


# sort by linkage group and position
params.sort<-params[order(as.numeric(params$chr), params$pos),]
params.sort$lg <- paste("LG",params.sort$chr, sep="")
params.sort$rs <- paste(params.sort$lg, params.sort$pos, sep="-")
params.sort$gene_name <- sapply(strsplit(params.sort$gene,split = ';'), `[`, 2)
params.sort$gene_name <- gsub('Name=','',params.sort$gene_name)
params.sort$eff <- abs(params.sort$beta*params.sort$gamma)

params.sort1 <- params.sort[abs(params.sort$beta) > 0,] # only variants with some effect
params.sort <- params.sort1
######################### plot of sum_gamma



# set up empty plot
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.ylim.png", sep=""), width=11.7,height=8.3,units="in",res=200)
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.png", sep=""), width=12,height=6)


# set up empty plot
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0,max(params.sort$gamma)),ylab=paste0("Top SNPs ", trait, " (PIP)"),xlab="linkage group", xaxt="n")


# plot grey bands for chromosome/linkage groups

start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(params.sort[params.sort$chr==ch,])
  cat ("CH: ", ch, "\n")
  colour<-adjustcolor("lightgrey", alpha.f=0.1)
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=colour, border=colour)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  
  start<-start+size
}

# Add x axis labels
axis(side=1,at=lab.pos,labels=chrs,tick=F)


# plot PIP for insignificant variants

#x<-match(params.sort$gamma[params.sort$gamma<threshold],params.sort$gamma)
top <- params.sort[quantile(params.sort$gamma,0.9999)[[1]] & params.sort$gamma > 0.01,]
x<- match(params.sort$gamma[!params.sort$id %in% top$id],params.sort$gamma)

# PIP
y<-params.sort$gamma[!params.sort$id %in% top$id]

# sparse effect size, used for dot size
z<-(params.sort[!params.sort$id %in% top$id,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))
symbols(x,y,circles=z, bg="grey90",inches=1/26,fg=NULL,add=T)

# plot threshold line
#abline(h=threshold,lty=3,col="black")

# rank of high PIP variants across linkage groups with gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
z<-(params.sort[params.sort$id %in% top$id & !is.na(params.sort$gene_name),])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#E64B35",alpha.f = 0.8), inches=1/26,fg="black",add=T)

# rank of high PIP variants across linkage groups without gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
#z<-(params.sort[params.sort$gamma>=threshold,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#E64B35",alpha.f = 0.8), inches=1/26,fg=NA,add=T)

# add label high PIP variants
#text(x,y,labels=params.sort$gene[params.sort$gamma>=threshold], adj=c(0,0), cex=0.4)
top_20_rows <- head(top[!is.na(top$gene_name), ][order(-top$gamma[!is.na(top$gene_name)]), ], 20)


labels <- top_20_rows$gene_name
x1<-match(top_20_rows$gamma,params.sort$gamma)
y1<-top_20_rows$gamma

addTextLabels(x1,y1,tolower(labels), col.label="black", cex.label = 1.1,  col.background = NULL, lty=1)

#dev.off()

#tophits <- params.sort[params.sort$gamma>=threshold,]
#write.table(tophits, paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, "99.9.PIPoutlier.SNP.txt", sep=""), quote=F, row.names=F, sep="\t")

#### 2e


# Reads input parameters
args<-commandArgs(trailingOnly=T)
#trait=args[1]
trait=c("RedChest")

# library to speed up loading of big tables
library(data.table)

# Load parameters output

params<-fread(paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".mean.params.anno.txt", sep=""),header=T,sep="\t", data.table=F)


params$pos <- params$start
params$chr <- as.numeric(gsub("chr","", params$seqnames))


# get list of linkage groups/chromosomes
chrs<-sort(as.numeric(unique(params$chr)))

# add linkage group column (chr)
chr<-gsub("lg|_.+","",params$chr)
params["chr"]<-chr


# sort by linkage group and position
params.sort<-params[order(as.numeric(params$chr), params$pos),]
params.sort$lg <- paste("LG",params.sort$chr, sep="")
params.sort$rs <- paste(params.sort$lg, params.sort$pos, sep="-")
params.sort$gene_name <- sapply(strsplit(params.sort$gene,split = ';'), `[`, 2)
params.sort$gene_name <- gsub('Name=','',params.sort$gene_name)
params.sort$eff <- abs(params.sort$beta*params.sort$gamma)

params.sort1 <- params.sort[abs(params.sort$beta) > 0,] # only variants with some effect
params.sort <- params.sort1
######################### plot of sum_gamma



# set up empty plot
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.ylim.png", sep=""), width=11.7,height=8.3,units="in",res=200)
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.png", sep=""), width=12,height=6)


# set up empty plot
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0,max(params.sort$gamma)),ylab=paste0("Top SNPs ", trait, " (PIP)"),xlab="linkage group", xaxt="n")


# plot grey bands for chromosome/linkage groups

start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(params.sort[params.sort$chr==ch,])
  cat ("CH: ", ch, "\n")
  colour<-adjustcolor("lightgrey", alpha.f=0.1)
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=colour, border=colour)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  
  start<-start+size
}

# Add x axis labels
axis(side=1,at=lab.pos,labels=chrs,tick=F)


# plot PIP for insignificant variants

#x<-match(params.sort$gamma[params.sort$gamma<threshold],params.sort$gamma)
top <- params.sort[quantile(params.sort$gamma,0.9999)[[1]] & params.sort$gamma > 0.01,]
x<- match(params.sort$gamma[!params.sort$id %in% top$id],params.sort$gamma)

# PIP
y<-params.sort$gamma[!params.sort$id %in% top$id]

# sparse effect size, used for dot size
z<-(params.sort[!params.sort$id %in% top$id,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))
symbols(x,y,circles=z, bg="grey90",inches=1/26,fg=NULL,add=T)

# plot threshold line
#abline(h=threshold,lty=3,col="black")

# rank of high PIP variants across linkage groups with gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
z<-(params.sort[params.sort$id %in% top$id & !is.na(params.sort$gene_name),])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#E64B35",alpha.f = 0.8), inches=1/26,fg="black",add=T)

# rank of high PIP variants across linkage groups without gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
#z<-(params.sort[params.sort$gamma>=threshold,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#E64B35",alpha.f = 0.8), inches=1/26,fg=NA,add=T)

# add label high PIP variants
#text(x,y,labels=params.sort$gene[params.sort$gamma>=threshold], adj=c(0,0), cex=0.4)
top_20_rows <- head(top[!is.na(top$gene_name), ][order(-top$gamma[!is.na(top$gene_name)]), ], 20)


labels <- top_20_rows$gene_name
x1<-match(top_20_rows$gamma,params.sort$gamma)
y1<-top_20_rows$gamma

addTextLabels(x1,y1,tolower(labels), col.label="black", cex.label = 1.1,  col.background = NULL, lty=1)

#dev.off()

#tophits <- params.sort[params.sort$gamma>=threshold,]
#write.table(tophits, paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, "99.9.PIPoutlier.SNP.txt", sep=""), quote=F, row.names=F, sep="\t")




#### 2f


# Reads input parameters
args<-commandArgs(trailingOnly=T)
#trait=args[1]
trait=c("MelanicBody")

# library to speed up loading of big tables
library(data.table)

# Load parameters output

params<-fread(paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".mean.params.anno.txt", sep=""),header=T,sep="\t", data.table=F)


params$pos <- params$start
params$chr <- as.numeric(gsub("chr","", params$seqnames))


# get list of linkage groups/chromosomes
chrs<-sort(as.numeric(unique(params$chr)))

# add linkage group column (chr)
chr<-gsub("lg|_.+","",params$chr)
params["chr"]<-chr


# sort by linkage group and position
params.sort<-params[order(as.numeric(params$chr), params$pos),]
params.sort$lg <- paste("LG",params.sort$chr, sep="")
params.sort$rs <- paste(params.sort$lg, params.sort$pos, sep="-")
params.sort$gene_name <- sapply(strsplit(params.sort$gene,split = ';'), `[`, 2)
params.sort$gene_name <- gsub('Name=','',params.sort$gene_name)
params.sort$eff <- abs(params.sort$beta*params.sort$gamma)

params.sort1 <- params.sort[abs(params.sort$beta) > 0,] # only variants with some effect
params.sort <- params.sort1
######################### plot of sum_gamma



# set up empty plot
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.ylim.png", sep=""), width=11.7,height=8.3,units="in",res=200)
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.png", sep=""), width=12,height=6)


# set up empty plot
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0,max(params.sort$gamma)),ylab=paste0("Top SNPs ", trait, " (PIP)"),xlab="linkage group", xaxt="n")


# plot grey bands for chromosome/linkage groups

start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(params.sort[params.sort$chr==ch,])
  cat ("CH: ", ch, "\n")
  colour<-adjustcolor("lightgrey", alpha.f=0.1)
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=colour, border=colour)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  
  start<-start+size
}

# Add x axis labels
axis(side=1,at=lab.pos,labels=chrs,tick=F)


# plot PIP for insignificant variants

#x<-match(params.sort$gamma[params.sort$gamma<threshold],params.sort$gamma)
top <- params.sort[quantile(params.sort$gamma,0.9999)[[1]] & params.sort$gamma > 0.01,]
x<- match(params.sort$gamma[!params.sort$id %in% top$id],params.sort$gamma)

# PIP
y<-params.sort$gamma[!params.sort$id %in% top$id]

# sparse effect size, used for dot size
z<-(params.sort[!params.sort$id %in% top$id,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))
symbols(x,y,circles=z, bg="grey90",inches=1/26,fg=NULL,add=T)

# plot threshold line
#abline(h=threshold,lty=3,col="black")

# rank of high PIP variants across linkage groups with gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
z<-(params.sort[params.sort$id %in% top$id & !is.na(params.sort$gene_name),])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#E64B35",alpha.f = 0.8), inches=1/26,fg="black",add=T)

# rank of high PIP variants across linkage groups without gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
#z<-(params.sort[params.sort$gamma>=threshold,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#E64B35",alpha.f = 0.8), inches=1/26,fg=NA,add=T)

# add label high PIP variants
#text(x,y,labels=params.sort$gene[params.sort$gamma>=threshold], adj=c(0,0), cex=0.4)
top_20_rows <- head(top[!is.na(top$gene_name), ][order(-top$gamma[!is.na(top$gene_name)]), ], 20)


labels <- top_20_rows$gene_name
x1<-match(top_20_rows$gamma,params.sort$gamma)
y1<-top_20_rows$gamma

addTextLabels(x1,y1,tolower(labels), col.label="black", cex.label = 1.1,  col.background = NULL, lty=1)

#dev.off()

#tophits <- params.sort[params.sort$gamma>=threshold,]
#write.table(tophits, paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, "99.9.PIPoutlier.SNP.txt", sep=""), quote=F, row.names=F, sep="\t")




#### 2g


# Reads input parameters
args<-commandArgs(trailingOnly=T)
#trait=args[1]
trait=c("BlueBody")

# library to speed up loading of big tables
library(data.table)

# Load parameters output

params<-fread(paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".mean.params.anno.txt", sep=""),header=T,sep="\t", data.table=F)


params$pos <- params$start
params$chr <- as.numeric(gsub("chr","", params$seqnames))


# get list of linkage groups/chromosomes
chrs<-sort(as.numeric(unique(params$chr)))

# add linkage group column (chr)
chr<-gsub("lg|_.+","",params$chr)
params["chr"]<-chr


# sort by linkage group and position
params.sort<-params[order(as.numeric(params$chr), params$pos),]
params.sort$lg <- paste("LG",params.sort$chr, sep="")
params.sort$rs <- paste(params.sort$lg, params.sort$pos, sep="-")
params.sort$gene_name <- sapply(strsplit(params.sort$gene,split = ';'), `[`, 2)
params.sort$gene_name <- gsub('Name=','',params.sort$gene_name)
params.sort$eff <- abs(params.sort$beta*params.sort$gamma)

params.sort1 <- params.sort[abs(params.sort$beta) > 0,] # only variants with some effect
params.sort <- params.sort1
######################### plot of sum_gamma



# set up empty plot
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.ylim.png", sep=""), width=11.7,height=8.3,units="in",res=200)
#png(file=paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, ".pip_plot.SNP.png", sep=""), width=12,height=6)


# set up empty plot
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0,max(params.sort$gamma)),ylab=paste0("Top SNPs ", trait, " (PIP)"),xlab="linkage group", xaxt="n")


# plot grey bands for chromosome/linkage groups

start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(params.sort[params.sort$chr==ch,])
  cat ("CH: ", ch, "\n")
  colour<-adjustcolor("lightgrey", alpha.f=0.1)
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=colour, border=colour)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  
  start<-start+size
}

# Add x axis labels
axis(side=1,at=lab.pos,labels=chrs,tick=F)


# plot PIP for insignificant variants

#x<-match(params.sort$gamma[params.sort$gamma<threshold],params.sort$gamma)
top <- params.sort[quantile(params.sort$gamma,0.9999)[[1]] & params.sort$gamma > 0.01,]
x<- match(params.sort$gamma[!params.sort$id %in% top$id],params.sort$gamma)

# PIP
y<-params.sort$gamma[!params.sort$id %in% top$id]

# sparse effect size, used for dot size
z<-(params.sort[!params.sort$id %in% top$id,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))
symbols(x,y,circles=z, bg="grey90",inches=1/26,fg=NULL,add=T)

# plot threshold line
#abline(h=threshold,lty=3,col="black")

# rank of high PIP variants across linkage groups with gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & !is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
z<-(params.sort[params.sort$id %in% top$id & !is.na(params.sort$gene_name),])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#E64B35",alpha.f = 0.8), inches=1/26,fg="black",add=T)

# rank of high PIP variants across linkage groups without gene names
x<-match(params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)],params.sort$gamma)

y<-params.sort$gamma[params.sort$id %in% top$id & is.na(params.sort$gene_name)]
# sparse effect size, used for dot size
#z<-(params.sort[params.sort$gamma>=threshold,])$eff
#z<-1/abs(log(z))
z <- rep(0.01, length(x))  
symbols(x,y, circles=z, bg=adjustcolor("#E64B35",alpha.f = 0.8), inches=1/26,fg=NA,add=T)

# add label high PIP variants
#text(x,y,labels=params.sort$gene[params.sort$gamma>=threshold], adj=c(0,0), cex=0.4)
top_20_rows <- head(top[!is.na(top$gene_name), ][order(-top$gamma[!is.na(top$gene_name)]), ], 20)


labels <- top_20_rows$gene_name
x1<-match(top_20_rows$gamma,params.sort$gamma)
y1<-top_20_rows$gamma

addTextLabels(x1,y1,tolower(labels), col.label="black", cex.label = 1.1,  col.background = NULL, lty=1)

#dev.off()

#tophits <- params.sort[params.sort$gamma>=threshold,]
#write.table(tophits, paste0("all.SNPs.filt.colour.imp.bslmm.colour.", trait, "99.9.PIPoutlier.SNP.txt", sep=""), quote=F, row.names=F, sep="\t")


dev.off()



