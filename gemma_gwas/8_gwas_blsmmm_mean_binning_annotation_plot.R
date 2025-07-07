#gwas_gemma_20kb_binning_annotation.R
#pooja.singh09@gmail.com
#this plots sum of gamma/pip in 20kb bins and add's cutoffs of 99, 95, and 99.9 percentiles




setwd('/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_dentition/bslmm')

options(echo=TRUE)
library(zoo) 
library(basicPlotteR)

# Reads input parameters
args<-commandArgs(trailingOnly=T)
trait=args[1]
#trait=c("YellowFlank")


# library to speed up loading of big tables
library(data.table)

# Load parameters output

params<-fread(paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.", trait, ".mean.params.20kb.anno.txt", sep=""),header=T,sep="\t", data.table=F)


params$pos <- (params$start + params$end)/2
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
params.sort$eff_20kb <- abs(params.sort$sum_gamma_20kb*params.sort$sum_beta_20kb) # variants with sparse effect
params.sort$gene_name <- sapply(strsplit(params.sort$gene,split = ';'), `[`, 2)
params.sort$gene_name <- gsub('Name=','',params.sort$gene_name)

pip_cutoff <- params.sort[params.sort$sum_gamma_20kb > 0.01,]
write.table(pip_cutoff, paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.", trait, ".PIPcutoff0.01.txt", sep=""), quote=F, row.names=F, sep="\t")

pip_cutoff <- params.sort[params.sort$sum_gamma_20kb > 0.05,]
write.table(pip_cutoff, paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.", trait, ".PIPcutoff0.05.txt", sep=""), quote=F, row.names=F, sep="\t")

################################################# plot of sum_gamma



# set up empty plot
png(file=paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.", trait, ".20kbpip_plot.png", sep=""), width=11.7,height=8.3,units="in",res=200)
#svg(file=paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.", trait, ".20kbpip_plot.svg", sep=""), width=11.7,height=8.3)


# set up empty plot
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0,max(params.sort$sum_gamma_20kb)),ylab=paste0("20Kb window PIP sum", trait),xlab="linkage group", xaxt="n", main=paste0(trait, "top GWAS hits BSLMM PIP 20KB win"))


# plot grey bands for chromosome/linkage groups

start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(params.sort[params.sort$chr==ch,])
  cat ("CH: ", ch, "\n")
  dentition<-"light grey"
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=dentition, border=dentition)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  
  start<-start+size
}

# Add x axis labels
axis(side=1,at=lab.pos,labels=chrs,tick=F)


# highligh high PIP variants

threshold <- quantile(params.sort$sum_gamma_20kb,0.999)[[1]]
abline(h=threshold,lty=3,col="black")

# plot PIP for insignificant variants

# rank of high PIP variants across linkage groups
x<-match(params.sort$sum_gamma_20kb[params.sort$sum_gamma_20kb<threshold],params.sort$sum_gamma_20kb)
# PIP
y<-params.sort$sum_gamma_20kb[params.sort$sum_gamma_20kb<threshold]
# sparse effect size, used for dot size
z<-(params.sort[params.sort$sum_gamma_20kb<threshold,])$eff_20kb
z<-1/abs(log(z))

symbols(x,y,circles=z, bg="grey",inches=1/8,fg=NULL,add=T)

# plot threshold line
abline(h=threshold,lty=3,col="black")
# rank of high PIP variants across linkage groups
x<-match(params.sort$sum_gamma_20kb[params.sort$sum_gamma_20kb>=threshold],params.sort$sum_gamma_20kb)
# PIP
y<-params.sort$sum_gamma_20kb[params.sort$sum_gamma_20kb>=threshold]
# sparse effect size, used for dot size
z<-(params.sort[params.sort$sum_gamma_20kb>=threshold,])$eff_20kb
z<-1/abs(log(z))
  
symbols(x,y,circles=z, bg="orange",inches=1/8,fg=NULL,add=T)

# plot threshold line


# add label high PIP variants
#text(x,y,labels=params.sort$gene[params.sort$sum_gamma_20kb>=threshold], adj=c(0,0), cex=0.4)
labels=params.sort$gene_name[params.sort$sum_gamma_20kb>=threshold]
addTextLabels(x,y,labels, col.label="black", cex.label = 0.8,  col.background = NULL, lty=1)

dev.off()

tophits <- params.sort[params.sort$sum_gamma_20kb>=threshold,]
write.table(tophits, paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.", trait, "99.9.PIPoutlier.txt", sep=""), quote=F, row.names=F, sep="\t")

