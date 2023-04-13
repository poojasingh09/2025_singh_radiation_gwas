# pooja.singh09@gmail.com
# making manhattan plots of the most significant hits from the gwas mlmm or lmm
# march 2023

library(stringr)

args<-commandArgs(trailingOnly=T)
trait=args[1]

a1 <- read.table(paste0("all.SNPs.filt.jaw.imp.lmm.jaw.", trait, ".assoc.top1percent.anno.txt"), header=T)
a1$lg <- as.numeric(sapply(strsplit(a1$id,"_"), `[`, 1))
a <- a1[order(a1$lg, a1$pos),]
a$logs <- -log10(a$p_score)
test <- str_split(a$gene, "_\\(LOC") ### edit this based on your annotation
a$gene_name <- unlist(lapply(test, "[[", 1))
b <- a[!is.na(a$gene_name),] ## removes SNPs with no gene annotation
a <- b

ylim_l <- min(a$logs)
ylim_u <- max(a$logs)

png(paste0("all.SNPs.filt.jaw.imp.lmm.jaw.", trait, ".assoc.top1percent.anno.manhattan.png"), width=11.7,height=8.3,units="in",res=200)

# set up empty plot

plot(-1,-1,xlim=c(0,nrow(a)),ylim=c(ylim_l,ylim_u),ylab="-log[p.value]",xlab="chromosome", xaxt="n", main=paste0(trait," top GWAS hits LVRS LMM \n symbol size = AF"))


# plot grey bands for chromosome/linkage groups

chrs <- 1:22
start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(a[a$lg == ch,])
  cat ("CH: ", ch, "\n")
  jaw<-"light grey"
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(ylim_l,ylim_u,ylim_u,ylim_l,ylim_l), col=jaw, border=jaw)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  
  start<-start+size
}

# Add variants outside linkage groups
#chrs<-c(chrs,"NA")
#size<-nrow(a[a$chr=="NA",])
#lab.pos<-c(lab.pos, start+size/2)


# Add x axis labels
axis(side=1,at=lab.pos,labels=chrs,tick=F)

# plot PIP for all variants

# rank of variants across linkage groups
x<-seq(1,length(a$logs),1)
rownames(a) <- x 

# use pvalue for y axis
y<- a$logs
# use allele frequency for DOT size
z<-a$af
# log-transform to enhance visibility
z[z==0]<-0.00000000001
z<-1/abs(log(z))
# plot
symbols(x,y,circles=z, bg="darkgrey",inches=1/10, fg=NULL,add=T)


# highligh high variants (top 1 percentile of -log pvalue set by threshold)

threshold <- quantile(a$logs, prob=c(0.99))

# plot threshold line
abline(h=threshold,lty=3,col="dark grey")

# add label high pvalue variants

high <- a[a$logs>=threshold,]
high$z <- high$af
high$z[high$z==0]<-0.00000000001
high$zz<-1/abs(log(high$z))
symbols(rownames(high),high$logs,circles=high$z, bg="palevioletred1",col="black", inches=1/10,fg=NULL,add=T)

text(rownames(high),high$logs,labels=high$gene_name, adj=c(0,0), cex=0.4)


# close device
dev.off()
