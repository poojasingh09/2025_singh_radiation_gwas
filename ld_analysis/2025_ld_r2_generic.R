


## pooja.singh09@gmail.com
## feb 2025
## calculate interchronmosome LD
## https://eacooper400.github.io/gen8900/exercises/ld-1.html


setwd("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/ld")
library(dplyr)

####################@@@@@@@@@@@FUNCTIONS

### Calculate Allele Frequencies from Genotype Counts
allele.freq <- function(genotypeCounts) {
  n = sum(genotypeCounts) - genotypeCounts["NN"]
  p = ((2*genotypeCounts["AA"]) + genotypeCounts["Aa"])/(2*n)
  q = 1-p
  freqs = c(p,q)
  names(freqs) = c("p", "q")
  return(freqs)
}
### Calculate the LD parameter R-squared from 2 rows of a VCF file
calc_r2 <- function(row1, row2) {
  g1 = get.field(row1[10:length(row1)], row1[9], "GT")
  g2 = get.field(row2[10:length(row2)], row2[9], "GT")
  pA = unname(allele.freq(count.genotypes(g1))["p"])
  pB = unname(allele.freq(count.genotypes(g2))["p"])
  h = get.haplotypes(g1, g2)
  pAB = (length(h[h=="00"]))/(length(h))
  D = pAB - (pA * pB)
  rsq = (D**2)/(pA*(1-pA)*pB*(1-pB))
  return(rsq)
}
### Count the AA, Aa, and aa genotypes in a sample
count.genotypes <- function(genotypes) {
  genotypes = gsub("(\\||/)", "", genotypes) 
  gen.patterns = c("00", "01", "10", "11", "..") 
  my.counts=table(factor(genotypes, levels=gen.patterns)) 
  final.counts = c(my.counts[1], (my.counts[2]+my.counts[3]), my.counts[4:5]) 
  names(final.counts) = c("AA", "Aa", "aa", "NN") 
  return(final.counts)
}
### Count the number of derived alleles for a VCF row
derivedCount <- function(row) {
  row=as.vector(row, mode="character")
  x=count.genotypes(get.field(row[10:length(row)], row[9], "GT"))
  dc=(2*x["aa"])+x["Aa"]
  return(unname(dc))
}
### Calculate Nucleotide Divergence (Dxy)
dxy <- function(vcf1, vcf2, perBP=TRUE) {
  g1=t(apply(vcf1, 1, function(x) get.field(x[10:length(x)], x[9], "GT")))
  g2=t(apply(vcf2, 1, function(x) get.field(x[10:length(x)], x[9], "GT")))
  af1=apply(g1, 1, function(x) allele.freq(count.genotypes(x))["p"])
  af2=apply(g2, 1, function(x) allele.freq(count.genotypes(x))["p"])
  ## Let x be the allele frequency (p) in pop1 * (1-p) in Pop2
  x = af1 * (1-af2)
  ## Let y be the allele frequency (p) in pop2 * (1-p) in Pop1
  y = af2 * (1-af1)
  dxy=sum((x+y))
  if (perBP) {
    c = unique(vcf1$CHROM)
    s = sapply(c, function(x,y) min(y[which(y$CHROM==x),2]), y=vcf1)
    e = sapply(c, function(x,y) max(y[which(y$CHROM==x),2]), y=vcf1)
    bp=sum(e-s)
    return(dxy/bp)
  } else { 
    return(dxy) 
  }
}
### Calculate Expected Het. and return, p, n, and H
expected.het <- function(genotypes) {
  obs.counts = count.genotypes(genotypes)
  n = sum(obs.counts) - obs.counts["NN"]
  freqs = allele.freq(obs.counts)
  Hexp = 2 * freqs[1] * freqs[2]
  res = c(freqs["p"], n, Hexp)
  res=as.numeric(unname(res))
  return(res)
}
### Determine which SNPs are Polymorphic vs Fixed in 2 species
fixed.poly <- function(vcf1, vcf2) {
  g1=t(apply(vcf1, 1, function(x) get.field(x[10:length(x)], x[9], "GT")))
  g2=t(apply(vcf2, 1, function(x) get.field(x[10:length(x)], x[9], "GT")))
  af1=apply(g1, 1, function(x) allele.freq(count.genotypes(x))["p"])
  af2=apply(g2, 1, function(x) allele.freq(count.genotypes(x))["p"])
  res = rep("Polymorphic", nrow(g1))
  res[which(abs(af1-af2)==1)] = "Fixed"
  return(res)
}
### Get a Specified Field From a VCF Sample/Genotype String
get.field <- function(samples, format, fieldName) {
  x=strsplit(samples, split=":")
  fields=unlist(strsplit(format, split=":")) 
  i=which(fields==fieldName)
  if (!(fieldName %in% fields)) stop('fieldName not found in format fields') 
  return(sapply(x, `[[`, i)) 
}
### Get the HAPLOTYPES for a pair of genotype strings
get.haplotypes <- function(genotypes1, genotypes2) {
  a1 = gsub("\\|", "", genotypes1) 
  a2 = gsub("\\|", "", genotypes2)
  a1=unlist(strsplit(paste0(a1, collapse=""), split="")) 
  a2=unlist(strsplit(paste0(a2, collapse=""), split=""))
  haps = paste0(a1,a2)
  return(haps)
}
### Calculate minor allele frequency
maf <- function(vcf.row) {
  temp=as.vector(vcf.row, mode="character")
  af=allele.freq(count.genotypes(get.field(temp[10:length(temp)], temp[9], "GT")))
  maf=min(unname(af))
  return(maf)
}
### Calculate Nucleotide Diversity (Pi)
pi.diversity <- function(vcf, perBP=TRUE) {
  J=apply(vcf, 1, derivedCount)
  N=apply(vcf, 1, function(x) sum(count.genotypes(get.field(x[10:length(x)], x[9], "GT"))))
  C=2*N
  pi = sum((2*J*(C-J))/(C*(C-1)))
  if (perBP) {
    c = unique(vcf$CHROM)
    s = sapply(c, function(x,y) min(y[which(y$CHROM==x),2]), y=vcf)
    e = sapply(c, function(x,y) max(y[which(y$CHROM==x),2]), y=vcf)
    bp=sum(e-s)
    return(pi/bp)
  } else { return(pi) }
}	      
### Read in a VCF file as a table
read.vcf <- function(file, special.char="##", ...) {
  my.search.term=paste0(special.char, ".*")
  all.lines=readLines(file)
  clean.lines=gsub(my.search.term, "",  all.lines)
  clean.lines=gsub("#CHROM", "CHROM", clean.lines)
  read.table(..., text=paste(clean.lines, collapse="\n"))
}
### Calculate the variance for Tajima's D
variance.d <- function(n,S) {
  a1=sum(1/(seq(from=1, to=(n-1), by=1)))
  a2=sum(1/((seq(from=1, to=(n-1), by=1))**2))
  b1=(n+1)/(3*(n-1))
  b2=(2*((n**2)+n+3))/((9*n)*(n-1))
  c1=b1 - (1/a1)
  c2=b2-((n+2)/(a1*n)) + (a2/(a1**2))
  e1=c1/a1
  e2=c2/((a1**2)+a2)
  var=(e1*S) + (e2*S*(S-1))
  return(var)
}
### Calculate Waterson's theta
waterson.theta <- function(data, perBP=TRUE) {
  num.bp=nrow(data)
  check=apply(data, 1, FUN=maf)
  filter=data[which(check>0),]
  Sn=nrow(filter)
  n.ind=ncol(filter)-9
  i=seq(1, ((2*n.ind)-1))
  theta=Sn/sum(1/i)
  if (perBP) { return(theta/num.bp) }
  else { return(theta) }
}


#########@@@@@@@@@@@@@

##setwd

setwd('/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/ld')

## read vcf for top SNPs across all comparisons

my.data <- read.vcf("top_rwas_all_99.9_and_pip01.vcf", header=TRUE, stringsAsFactors=FALSE, check.names = FALSE)
core <- my.data[1:9]

### read in species pairs
samples <- read.table("lilo", header=F)

### species pair in consideration
my.data1 <- my.data[,names(my.data) %in% samples$V1]
my.data1 <- cbind(core, my.data1)

### read in top hits species pair in consideration
tophits <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/fst/allcomparisons_top_gwas_pip01_SNP_0.95fst_tophits_SNPoverlaponly.txt", header=T)
comparison <- c("stitch")
tophits_species <- tophits[tophits$comparison %in% comparison, ]
my.data2 <- merge(my.data1, tophits_species, by.x=c("CHROM", "POS"), by.y=c("seqnames", "start"))
dim(tophits_species)
dim(my.data2)
column_name <- "comparison"
my.data2 <- my.data2[, 1:(which(names(my.data2) == column_name) - 1)]
dim(my.data2)



### Extract the first 2 rows
row1=as.vector(my.data2[1,], mode="character")
row2=as.vector(my.data2[2,], mode="character")

### Calculate the distance !adjust for multichromosomes
#dist = as.numeric(row2[2]) - as.numeric(row1[2])

### Get pA and pB
g1 = get.field(row1[10:length(row1)], row1[9], "GT") # Using the fxn isn't necessary since I don't have other fields, but I'm going to have it in the code for the future when I might have other stuff there and still want to calc LD
g2 = get.field(row2[10:length(row2)], row2[9], "GT")

pA = unname(allele.freq(count.genotypes(g1))["p"])
pB = unname(allele.freq(count.genotypes(g2))["p"])

### Calculate pAB:
### First, get all of the HAPLOTYPES
get.haplotypes <- function(genotypes1, genotypes2) {
  a1 = gsub("\\|", "", genotypes1) # Get rid of the "|" symbol between alleles
  a2 = gsub("\\|", "", genotypes2)
  a1=unlist(strsplit(paste0(a1, collapse=""), split="")) # Collapse everything into 1 string, then break it apart again but now in a way such that each allele is a separate item in the list
  a2=unlist(strsplit(paste0(a2, collapse=""), split=""))
  haps = paste0(a1,a2)
  return(haps)
}

h=get.haplotypes(g1,g2)

### Now, just count the fraction of times "00" is the hap:
pAB = (length(h[h=="00"]))/length(h)

### Calculate D and r-squared
D = pAB - (pA * pB)
rsq = (D**2)/(pA*(1-pA)*pB*(1-pB))

### Check the other sites:
my.results = data.frame(matrix(0, nrow=nrow(my.data2)*nrow(my.data2), ncol=8))
colnames(my.results)=c("SNP1", "SNP2", "SNP1chr", "SNP1pos", "SNP2chr", "SNP2pos","D", "Rsq")
counter=1

for (i in 1:((nrow(my.data2))-1)) {
  row1=as.vector(my.data2[i,], mode="character")
  row1id <- paste(row1[1], row1[2], sep="_")
  row1chr <- row1[1]
  row1pos <- row1[2]
  g1 = get.field(row1[10:length(row1)], row1[9], "GT")
  pA = unname(allele.freq(count.genotypes(g1))["p"])
  
  for (j in (i+1):(nrow(my.data2))) {
    row2=as.vector(my.data2[j,], mode="character")
    row2id <- paste(row2[1], row2[2], sep="_")
    row2chr <- row2[1]
    row2pos <- row2[2]
    g2 = get.field(row2[10:length(row2)], row2[9], "GT")
    pB = unname(allele.freq(count.genotypes(g2))["p"])
    #dist = as.numeric(row2[2]) - as.numeric(row1[2])
    h=get.haplotypes(g1,g2)
    pAB = (length(h[h=="00"]))/length(h)
    D = pAB - (pA * pB)
    rsq = (D**2)/(pA*(1-pA)*pB*(1-pB))
    my.results[counter,] = c(row1id, row2id,row1chr, row1pos, row2chr, row2pos, D, rsq)
    counter = counter+1
  }
}

head(my.results)

dim(my.results)
my.results.clean <- my.results[my.results$SNP1 != 0 & my.results$Rsq != "NaN",]
dim(my.results.clean)
my.results.clean2 <- my.results.clean[my.results.clean$SNP1 != my.results.clean$SNP2,]
dim(my.results.clean2)
boxplot(as.numeric(my.results.clean2$Rsq), outline=F)
my.results.clean2$comparison <- paste0(comparison, "_top")

my.results.clean2.top <- my.results.clean2

head(my.results.clean2.top)

my.results.clean2.top_unique <- my.results.clean2.top[!duplicated(my.results.clean2.top), ]


###########################


## read vcf for Background SNPs across all comparisons

my.data <- read.vcf("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_dentition/ld/all.SNPs.filt2.colour.imp.random1m.recode.vcf", header=TRUE, stringsAsFactors=FALSE, check.names = FALSE)
my.data <- sample_n(my.data, 100)
core <- my.data[1:9]

### subset species pair in consideration
my.data1 <- my.data[,names(my.data) %in% samples$V1]
my.data1 <- cbind(core, my.data1)


### Extract the first 2 rows
row1=as.vector(my.data1[1,], mode="character")
row2=as.vector(my.data1[2,], mode="character")

### Calculate the distance !adjust for multichromosomes
#dist = as.numeric(row2[2]) - as.numeric(row1[2])

### Get pA and pB
g1 = get.field(row1[10:length(row1)], row1[9], "GT") # Using the fxn isn't necessary since I don't have other fields, but I'm going to have it in the code for the future when I might have other stuff there and still want to calc LD
g2 = get.field(row2[10:length(row2)], row2[9], "GT")

pA = unname(allele.freq(count.genotypes(g1))["p"])
pB = unname(allele.freq(count.genotypes(g2))["p"])

### Calculate pAB:
### First, get all of the HAPLOTYPES
get.haplotypes <- function(genotypes1, genotypes2) {
  a1 = gsub("\\|", "", genotypes1) # Get rid of the "|" symbol between alleles
  a2 = gsub("\\|", "", genotypes2)
  a1=unlist(strsplit(paste0(a1, collapse=""), split="")) # Collapse everything into 1 string, then break it apart again but now in a way such that each allele is a separate item in the list
  a2=unlist(strsplit(paste0(a2, collapse=""), split=""))
  haps = paste0(a1,a2)
  return(haps)
}

h=get.haplotypes(g1,g2)

### Now, just count the fraction of times "00" is the hap:
pAB = (length(h[h=="00"]))/length(h)

### Calculate D and r-squared
D = pAB - (pA * pB)
rsq = (D**2)/(pA*(1-pA)*pB*(1-pB))

### Check the other sites:
my.results = data.frame(matrix(0, nrow=nrow(my.data1)*nrow(my.data1), ncol=8))
colnames(my.results)=c("SNP1", "SNP2", "SNP1chr", "SNP1pos", "SNP2chr", "SNP2pos","D", "Rsq")
counter=1

for (i in 1:((nrow(my.data1))-1)) {
  row1=as.vector(my.data1[i,], mode="character")
  row1id <- paste(row1[1], row1[2], sep="_")
  row1chr <- row1[1]
  row1pos <- row1[2]
  g1 = get.field(row1[10:length(row1)], row1[9], "GT")
  pA = unname(allele.freq(count.genotypes(g1))["p"])
  
  for (j in (i+1):(nrow(my.data1))) {
    row2=as.vector(my.data1[j,], mode="character")
    row2id <- paste(row2[1], row2[2], sep="_")
    row2chr <- row2[1]
    row2pos <- row2[2]
    g2 = get.field(row2[10:length(row2)], row2[9], "GT")
    pB = unname(allele.freq(count.genotypes(g2))["p"])
    #dist = as.numeric(row2[2]) - as.numeric(row1[2])
    h=get.haplotypes(g1,g2)
    pAB = (length(h[h=="00"]))/length(h)
    D = pAB - (pA * pB)
    rsq = (D**2)/(pA*(1-pA)*pB*(1-pB))
    my.results[counter,] = c(row1id, row2id,row1chr, row1pos, row2chr, row2pos, D, rsq)
    counter = counter+1
    print(counter)
  }
}

head(my.results)

dim(my.results)
my.results.clean <- my.results[my.results$SNP1 != 0 & my.results$Rsq != "NaN",]
dim(my.results.clean)
my.results.clean2 <- my.results.clean[my.results.clean$SNP1 != my.results.clean$SNP2,]
dim(my.results.clean2)
boxplot(as.numeric(my.results.clean2$Rsq), outline=F)
my.results.clean2$comparison <- paste0(comparison, "_bg")
my.results.clean2.bg <- my.results.clean2
head(my.results.clean2.bg)

my.results.clean2.bg_unique <- my.results.clean2.bg[!duplicated(my.results.clean2.bg), ]

## save results
write.table(rbind(my.results.clean2.bg_unique, my.results.clean2.top_unique), paste0(comparison, "_LD_r2_results.txt"), sep="\t", quote=F, row.names=F)

## both plots
png(paste0(comparison, "_LD_boxplot.png"))
boxplot(as.numeric(my.results.clean2.bg_unique$Rsq), as.numeric(my.results.clean2.top_unique$Rsq), names=c("genome-wide","top"), ylab="LD (R^2)", xlab=comparison, outline=F, horizontal=F, col=c("grey","pink"), notch=F, ylim=c(0,1))
dev.off()


stats <- wilcox.test(as.numeric(my.results.clean2.bg_unique$Rsq), as.numeric(my.results.clean2.top_unique$Rsq))
write.csv(stats$p.value, paste0(comparison, "_LD_r2_results_wilcox.txt"), quote=F, row.names=F)
