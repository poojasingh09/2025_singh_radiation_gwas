##pooja.singh09@gmail.com
# this code summarises gemma blsmm hyperparameters output
#adapted from http://romainvilloutreix.alwaysdata.net/romainvilloutreix/wp-content/uploads/2017/01/gwas_gemma-2017-01-17.pdf
# usage Rscript --vanilla gwas_gemma_hyperparam_summary.R

options(echo=TRUE)

library(ggplot2)
library(gridExtra)
library(dplyr)


# Reads input parameters
args<-commandArgs(trailingOnly=T)
trait=args[1]

print(trait)

#source is either markgwas or dentitiongwas
#trait is one of the traits in the sources

if(length(args)<1){
  stop(paste0("Aborted, not enough arguments given.\nUsage: ",usage,detail))
}

setwd("./")

# Load hyperparameter files

dataFiles <- lapply(Sys.glob(paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.", trait, ".*.hyp.txt", sep="")), read.table, header=T)


# Combine all hyperparameter files to calculate overall distributions

hyp.params <- bind_rows(dataFiles, .id = "column_label")
write.table(hyp.params, file=paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.", trait, ".hyp.combined.txt", sep=""), sep="\t", quote=F)


########################


# Get mean, median, and 95% ETPI of hyperparameters
# rho-> approximation to proportion of genetic variance explained by variants with major effect (PGE)
# rho=0 -> pure LMM, highly polygenic basis
# rho=1 -> pure BVSR, few major effect loci

rho<-c("rho",mean(hyp.params$rho),quantile(hyp.params$rho, probs=c(0.5,0.025,0.975)))

# pve -> proportion of phenotypic variance explained by the genotypes
pve<-c("PVE", mean(hyp.params$pve),quantile(hyp.params$pve, probs=c(0.5,0.025,0.975)))

# pge -> proportion of genetic variance explained by major effect loci
pge<-c("PGE",mean(hyp.params$pge),quantile(hyp.params$pge, probs=c(0.5,0.025,0.975)))

# pi -> proportion of variants with non-zero effects
pi<-c("pi",mean(hyp.params$pi),quantile(hyp.params$pi, probs=c(0.5,0.025,0.975)))

# n.gamma -> number of variants with major effect
n.gamma<-c("n.gamma",mean(hyp.params$n_gamma),quantile(hyp.params$n_gamma, probs=c(0.5,0.025,0.975)))


# make table of hyperparameters

hyp.params.table<-as.data.frame(rbind(rho,pve,pge,pi,n.gamma),row.names=F)
colnames(hyp.params.table)<-c("hyperparam", "mean","median","2.5%", "97.5%")

# show table
hyp.params.table
hyp.params.table$phenotpe <- trait
# write table to file
write.table(hyp.params.table, file=paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.", trait, ".hyp.combined.overallsummary.txt", sep=""), sep="\t", quote=F)


# plot traces and distributions of hyperparameters

pdf(file=paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.", trait, ".hyp.combined-traces-hyperparam.pdf", sep=""), width=8.3,height=11.7)
layout(matrix(c(1,1,2,3,4,4,5,6), 4, 2, byrow = TRUE))

# PVE
# ------------------------------------------------------------------------------
plot(hyp.params$pve, type="l", ylab="PVE", main="PVE - trace")
hist(hyp.params$pve, main="PVE - posterior distribution", xlab="PVE")
plot(density(hyp.params$pve), main="PVE - posterior distribution", xlab="PVE")
# ------------------------------------------------------------------------------

# PGE
# ------------------------------------------------------------------------------
plot(hyp.params$pge, type="l", ylab="PGE", main="PGE - trace")
hist(hyp.params$pge, main="PGE - posterior distribution", xlab="PGE")
plot(density(hyp.params$pge), main="PGE - posterior distribution", xlab="PGE")
# ------------------------------------------------------------------------------

# pi
# ------------------------------------------------------------------------------
plot(hyp.params$pi, type="l", ylab="pi", main="pi")
hist(hyp.params$pi, main="pi", xlab="pi")
plot(density(hyp.params$pi), main="pi", xlab="pi")
# ------------------------------------------------------------------------------

# No gamma
# ------------------------------------------------------------------------------
plot(hyp.params$n_gamma, type="l", ylab="n_gamma", main="n_gamma - trace")
hist(hyp.params$n_gamma, main="n_gamma - posterior distribution", xlab="n_gamma")
plot(density(hyp.params$n_gamma), main="n_gamma - posterior distribution", xlab="n_gamma")
# ------------------------------------------------------------------------------
dev.off()



# plot distributions of hyperparameters for all MCMC runs

pdf(file=paste0("all.SNPs.filt.dentition.imp.bslmm.dentition.", trait, ".hyp.overall-hyperparam.pdf", sep=""), width=15,height=5)
#layout(matrix(c(1, 2, 3), 1,3, byrow=TRUE), respect = TRUE)

# PVE
# ------------------------------------------------------------------------------

p <- ggplot(data=hyp.params, mapping = aes(x = pve)) + theme_classic()
p <- p + geom_density(color="blue")
# ------------------------------------------------------------------------------



# PGE
# ------------------------------------------------------------------------------
p1 <- ggplot(data=hyp.params, mapping = aes(x = pge)) + theme_classic()
p1 <- p1 + geom_density(color="blue")
# ------------------------------------------------------------------------------



# No gamma
# ------------------------------------------------------------------------------
p2 <- ggplot(data=hyp.params, mapping = aes(x = n_gamma)) + theme_classic()
p2 <- p2 + geom_density(color="blue")
# ------------------------------------------------------------------------------


grid.arrange(p, p1, p2, nrow = 1)
dev.off()

