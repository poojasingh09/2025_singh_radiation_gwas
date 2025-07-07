#pooja.singh09@gmail.com
#march2023


library(qqman)

# set input file read

args<-commandArgs(trailingOnly=T)
trait=args[1]

# read file
a1 <- read.table(paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.",trait,".mean.params.anno.txt"), header=T)

## select only needed cols

a2 <- a1[,c(1,6,9,10,11)] # cols chr, pos, gamma, eff, id, gene
a2 <- a2[order(a2$chr, a2$pos),]
a2$id <- paste0(a2$chr, "_", a2$pos)


cutoff1 <- quantile(a2$gamma, probs=c(0.99))
print(cutoff1)
cutoff2 <- quantile(a2$gamma, probs=c(0.95))
print(cutoff2)
lim <- max(a2$gamma)

## plot manhattan with snp id of gemma gamma from 10 bslmm runs

png(file=paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.",trait,".mean.params.manhattan.snpid.png"), width=14,height=8.3,units="in",res=200)
manhattan(a2, chr = "chr", bp = "pos", p = "gamma", snp = "id", ylim=c(0,lim), ylab = paste0(trait, " PIP BLSMM"), col = c("gray10", "gray60"), suggestiveline = cutoff2, genomewideline = cutoff1, logp = FALSE, main=paste0(trait, " manhattan plot"), annotatePval = 0.01, annotateTop = TRUE)
dev.off()

## plot manhattan without snp id
png(file=paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.",trait,".mean.params.manhattan.png"), width=14,height=8.3,units="in",res=200)
manhattan(a2, chr = "chr", bp = "pos", p = "gamma", snp = "id", ylim=c(0,lim), ylab = paste0(trait, " PIP BLSMM"), col = c("gray10", "gray60"), suggestiveline = cutoff2, genomewideline =cutoff1, logp = FALSE, main=paste0(trait, " manhattan plot"))
dev.off()

## plot manhattan per chr
chrs <- 1:22

pdf(file=paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.",trait,".mean.params.manhattan.chrs.pdf"), onefile=TRUE, width=14,height=8.3)

for (ch in chrs){
	manhattan(subset(a2, chr==ch), chr = "chr", bp = "pos", p = "gamma", snp = "id", ylim=c(0,lim), ylab = paste0(trait, " PIP BLSMM"), col = c("gray10", "gray60"), suggestiveline = cutoff2, genomewideline = cutoff1, logp = FALSE, main=paste0(trait, " manhattan plot"), annotatePval = 0.01, annotateTop = FALSE)	      	    }
dev.off()

## plot qqplot

png(file=paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.",trait,".mean.params.qqplot.png"), width=11.7,height=8.3,units="in",res=200)
qq(a2$gamma, main = paste0(trait, " Q-Q plot of GWAS p-values"))
dev.off()
