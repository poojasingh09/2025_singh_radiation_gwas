##pooja.singh09@gmail.com
# this code summarises gemma blsmm hyperparameters output
#adapted from http://romainvilloutreix.alwaysdata.net/romainvilloutreix/wp-content/uploads/2017/01/gwas_gemma-2017-01-17.pdf
# usage Rscript --vanilla gwas_gemma_param_summary_plot.R


options(echo=TRUE)

# Reads input parameters
args<-commandArgs(trailingOnly=T)
trait=args[1]
#source=args[2]

#trait=c("MelanicBody")


#source is either markgwas or dentitiongwas
#trait is one of the traits in the sources


if(length(args)<1){
  stop(paste0("Aborted, not enough arguments given.\nUsage: ",usage,detail))
}

setwd("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_bcftools_pnye3_dentition_repeat/bslmm")
print(trait)

# library to speed up loading of big tables
library(data.table)

# Load parameters output

params<-fread(paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.", trait, ".mean.params.txt", sep=""),header=T,sep="\t", data.table=F)



# Get variants with sparse effect size on phenotypes 

# add sparse effect size (= beta * gamma) to data frame
params["eff"]<-abs(params$beta_mean*params$gamma_mean)

# get variants with effect size > 0
params.effects<-params[params$eff>0,]

# show number of variants with measurable effect
nrow(params.effects)

# sort by descending effect size
params.effects.sort<-params.effects[order(-params.effects$eff),]

# show top 10 variants with highest effect
head(params.effects.sort, 10)

# variants with the highest sparse effects

# top 1% variants (above 99% quantile)
top1<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.99),]
# top 0.1% variants (above 99.9% quantile)
top01<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.999),]
# top 0.01% variants (above 99.99% quantile)
top001<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.9999),]


# write tables
write.table(top1, file=paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.", trait, ".top1eff.dsv", sep=""), quote=F, row.names=F, sep="\t")
write.table(top01, file=paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.", trait, ".top0.1eff.dsv", sep=""), quote=F, row.names=F, sep="\t")
write.table(top001, file=paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.", trait, ".top0.01eff.dsv", sep=""), quote=F, row.names=F, sep="\t")



# Get variants with high Posterior Inclusion Probability (PIP) == gamma

# PIP is the frequency a variant is estimated to have a sparse effect in the MCMC

# sort variants by descending PIP
params.pipsort<-params[order(-params$gamma_mean),]

# Show top 10 variants with highest PIP
head(params.pipsort,10)

# sets of variants above a certain threshold, the high the pip, the better! Parameters with gamma values 0.01, 0.1, 0.25, and 0.5 have posterior inclusion probabilities of 1%, 10%, 25% and 50%
# the parameters that have a pip (or gamma ) that are greater than 0.01, .1 and 0.5 so 1%, 10% and 50% chance of being included in the model, the higher the PIP the higher the prop that SNP explains var. in pheno

# variants with effect in 1% MCMC samples or more
pip01<-params.pipsort[params.pipsort$gamma_mean>=0.01,]

# variants with effect in 1% MCMC samples or more
pip05<-params.pipsort[params.pipsort$gamma_mean>=0.05,]

# variants with effect in 10% MCMC samples or more
pip10<-params.pipsort[params.pipsort$gamma_mean>=0.10,]
# variants with effect in 25% MCMC samples or more
pip25<-params.pipsort[params.pipsort$gamma_mean>=0.25,]
# variants with effect in 50% MCMC samples or more
pip50<-params.pipsort[params.pipsort$gamma_mean>=0.50,]

# write tables
write.table(pip01, file=paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.", trait, ".pip01.dsv", sep=""), quote=F, row.names=F, sep="\t")
write.table(pip05, file=paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.", trait, ".pip05.dsv", sep=""), quote=F, row.names=F, sep="\t")
write.table(pip10, file=paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.", trait, ".pip10.dsv", sep=""), quote=F, row.names=F, sep="\t")
write.table(pip25, file=paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.", trait, ".pip25.dsv", sep=""), quote=F, row.names=F, sep="\t")
write.table(pip50, file=paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.", trait, ".pip50.dsv", sep=""), quote=F, row.names=F, sep="\t")



# plot variants PIPs across linkage groups/chromosomes

# Prepare data

# add linkage group column (chr)
chr<-gsub("lg|_.+","",params$chr)
params["chr"]<-chr


# sort by linkage group and position
params.sort<-params[order(as.numeric(params$chr), params$ps),]
params.sort$lg <- paste("LG",params.sort$chr, sep="")
params.sort$rs <- paste(params.sort$lg, params.sort$ps, sep="-")

####### annotate SNPs

library(plyranges)
library(dplyr)
wk_dir <- "."

sizes <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2023_bcftools_pnye3/pnye_v3.fa.fai.2col", header=FALSE, stringsAsFactors = F)
#genes <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2023_bcftools_pnye3/Bouillabase_Tilapia_NCBI_nyerereitrack.pundcross.gapsEstimated.genes2pnye2.bed", stringsAsFactors = FALSE)
genes <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/pnye_v3_assembly/punda.v3.annotation_v2.gene.bed.txt", stringsAsFactors = FALSE)


#colnames(sizes) <- c("seqnames", "end")
#sizes$start <- 0
#bins <- sizes %>% as_granges() %>% tile_ranges(width = 50000L)


### overlap with genes ### top EFF 1%

colnames(genes) <- c("seqnames", "start","end","gene")
head(genes)
genes <- genes %>% as_granges()

  snps <- params.sort
  colnames(snps) <- c("chr", "pos","alpha","beta","gamma", "eff", "lg", "id")
  snps$start <- snps$pos
  snps$end <- snps$pos
  snps$seqnames <- paste0("chr", snps$chr)
  snps <- snps %>% as_granges()
  snps1 <- join_overlap_inner(snps, genes)
  snps1 <- as(snps1, "data.frame")
  snps1 <- snps1[order(snps1$gamma, decreasing = TRUE),]
  write.table(snps1, paste(wk_dir,"all.SNPs.filt.dentition.imp.bslmm.dentition.",trait,".mean.params.anno.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = '\t')


## clean genes names
library(stringr)
params.sort <- snps1
test <- str_split(params.sort$gene, "_\\(LOC") ### edit this based on your annotation
params.sort$gene_name <- unlist(lapply(test, "[[", 1))


####################



# get list of linkage groups/chromosomes
chrs<-sort(as.numeric(unique(params$chr)))

# Plot to a png file because the number of dots is very high
# drawing this kind of plot over the network is very slow
# also opening vectorial files with many objects is slow


## sort again by chr, pos because earlier you sorted by gamma for the text file

params.sort<-params.sort[order(as.numeric(params.sort$chr), params.sort$pos),]

#set up png

png(file=paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.", trait, ".pip_plot.png", sep=""), width=11.7,height=8.3,units="in",res=200)

# set up empty plot
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0,max(params.sort$gamma)),ylab=paste0("PIP ", trait),xlab="linkage group", xaxt="n", main=paste0(trait, " top GWAS hits BSLMM"))



# plot grey bands for chromosome/linkage groups
chrs <- 1:22
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


#set threshold
threshold <- quantile(params.sort$gamma, prob=c(0.9999))

# plot threshold line
abline(h=threshold,lty=3,col="black")


## # plot PIP for below threshold variations

# rank of variants across linkage groups
x<-match(params.sort$gamma[params.sort$gamma<threshold],params.sort$gamma)
# PIP
y<-params.sort$gamma[params.sort$gamma<threshold]
# sparse effect size, used for dot size
z<-params.sort$eff[params.sort$gamma<threshold]
z<-1/abs(log(z))
symbols(x,y,circles=z, bg="grey",inches=1/10, fg=NULL,add=T)


# rank of high PIP variants across linkage groups
x<-match(params.sort$gamma[params.sort$gamma>=threshold],params.sort$gamma)
# PIP
y<-params.sort$gamma[params.sort$gamma>=threshold]
# sparse effect size, used for dot size
z<-params.sort$eff[params.sort$gamma>=threshold]
z<-1/abs(log(z))
  
symbols(x,y,circles=z, bg="green",inches=1/10,fg=NULL,add=T)


# add label high PIP variants
text(x,y,labels=params.sort$gene_name[params.sort$gamma>=threshold], adj=c(0,0), cex=0.4)


# close device
dev.off()
