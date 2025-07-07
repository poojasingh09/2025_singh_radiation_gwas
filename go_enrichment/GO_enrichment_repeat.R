## pooja.singh09@gmail.com
## feb 2025
## GO enrichment using Pnye v3 genome for Singh et al 

## set wd

setwd("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/go_enrichment/")

# load libraries

library(rtracklayer)
library(tidyverse)
library(tidyr)
library(dplyr)
library(topGO)
library(Rgraphviz)

# 
# ## build gene2go mappings
# 
# my_tags <- c("Name", "Parent", "Ontology_term")
# gff <- readGFF("punda.v3.annotation_v2.gff.GO", tags=my_tags)
# head(gff)
# gene2go1 <- unique(data.frame(unlist(paste(gff$Parent, gff$Ontology_term))))
# names(gene2go1) <- "mapping"
# gene2go <- unique(data.frame(gff[,c(10,11)]))
# 
# head(gene2go)
# head(gene2go1)
# 
# #write.table(gene2go, file="pnye3_gene2go.txt", quote=F, row.names = F)
# #write.table(gene2go1, file="pnye2_gene2go_alt.txt", quote=F, row.names = F, sep="\t")
# 
# ### read in gene2go (this was a bit of a headache to wrangle)
# 
# gene2go1 <- read.table("pnye3_gene2go_alt.txt", sep="\t")
# gene2go <- chop(gene2go1, V2)
# names(gene2go) <- c("gene", "go_id")
# fwrite(list(gene2go$V2), file = "b.csv")
# fwrite(list(gene2go$V1), file = "a.csv")
# #paste a.csv b.csv | sed 's/|G/, G/g'> pnye3_gene2go.txt  <- in shell

geneID2GO <- readMappings(file= "pnye3_gene2go.txt")
str(head(geneID2GO))

background_gwas <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_repeatcolour/bslmm/background_genes", header=F)

## update gene2go
new <- geneID2GO[names(geneID2GO) %in% background_gwas$V1]
geneID2GO <- new

##gwas background

geneNames <- names(geneID2GO)
head(geneNames)


files <- list.files(full.names = T, path="/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/snp_effect_overlaps", pattern = "\\_pip01.txt$")
files
traits <- c("bb", "cs", "dlb", "fb", "ir", "kd", "mb", "mlb","pc1", "pc2", "rc","tm","vb","yf")

counter=0
for (file in files) {
  counter=counter+1
  traitname <- traits[counter]
  print(traitname)
  myInterestingGenes <- read.table(files[counter], header=T)
  dim(myInterestingGenes)
  head(myInterestingGenes)
  myInterestingGenes3 <- as.character(sapply(strsplit(myInterestingGenes$gene, ";"), `[[`, 1))
  myInterestingGenes <- gsub("ID=", "", myInterestingGenes3)

  head(myInterestingGenes)
  length(myInterestingGenes)
  myInterestingGenes <- myInterestingGenes[!is.na(myInterestingGenes)]
  geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList) <- geneNames
  str(geneList)


  GOdata_BP <- new("topGOdata", ontology = "BP", allGenes = geneList , annot = annFUN.gene2GO, gene2GO = geneID2GO) 
  GOdata_BP
  numGenes(GOdata_BP)
  sg <- sigGenes(GOdata_BP)
  str(sg)
  
  myGOdata <- GOdata_BP
  
  # run the Fisher's exact tests
  resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")
  resultElim <- runTest(myGOdata, algorithm="elim", statistic="fisher")
  resultTopgo <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
  resultParentchild <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
  
  # see how many results we get where weight01 gives a P-value <= 0.05:
  #mysummary <- summary(attributes(resultTopgo)$score <= 0.01)
  #numsignif <- as.integer(mysummary[[3]])
  #dim(numsignif)
  
  # print out the top 'numsignif' results:
  allRes <- GenTable(myGOdata, classicFisher = resultClassic, elimFisher = resultElim, topgoFisher = resultTopgo, parentchildFisher = resultParentchild, orderBy = "topgoFisher", ranksOf = "classicFisher", numChar=1000)
  allRes$topgoFisher.padjust <- p.adjust(allRes$topgoFisher, method="fdr")
  
  write.table(allRes[allRes$topgoFisher.padjust < 0.05,], file = paste(traitname, "topgo_fisher_weight_tests_BP_q0.05.txt", sep="_"), quote = FALSE, row.names = FALSE, sep = "\t")
  write.table(allRes, file = paste(traitname, "topgo_fisher_weight_tests_BP_full.txt", sep="_"), quote = FALSE, row.names = FALSE, sep = "\t")
  
  
  
  
  ### genes associated with go enriched terms
  
  
  ## read in annotation
  
  annos <- read.table("punda.v3.annotation_v2.gene.bed.edited.txt")
  anno <- annos[,c(5,4)]
  head(anno)
  colnames(anno) <- c("gene_id","gene_info")
  anno$gene_id <- gsub("ID=","" , anno$gene_id ,ignore.case = TRUE)
  
  sig_genes <- data.frame(myInterestingGenes)
  big <- merge(sig_genes, anno, by.x="myInterestingGenes", by.y="gene_id")
  
  gene2go_alt <- read.table("pnye3_gene2go_alt.txt")
  go_of_gwas_hits <- merge(gene2go_alt, big, by.x="V1", by.y="myInterestingGenes")
  write.table(go_of_gwas_hits, file = paste(traitname, "gene_name_GOid.txt", sep="_"), quote = FALSE, row.names = FALSE, sep = "\t")
  
  my_list <- list()
  
  myterms = (allRes[allRes$topgoFisher.padjust < 0.05,])$GO.ID
  
  if (length(myterms) > 0){
  mygenes <- genesInTerm(myGOdata, myterms)
  for (j in 1:length(myterms))
  {
    myterm <- myterms[j]
    mygenesforterm <- mygenes[myterm][[1]]
    myfactor <- mygenesforterm %in% myInterestingGenes # find the genes that are in the list of genes of interest
    mygenesforterm2 <- mygenesforterm[myfactor == TRUE] 
    mygenesforterm3 <- data.frame(mygenesforterm2)
    mygenesforterm3$GO_id <- myterms[j]
    colnames(mygenesforterm3) <- c("gene_id", "GO_id")
    full <- merge(mygenesforterm3, anno, by="gene_id")
    my_list[[j]] <- full
  }
  
  # Write out single data frames
  lst2 <- lapply(my_list,function(x) cbind(rowname=rownames(x),x))
  df1 <- Reduce(function(x,y) merge(x,y,all=T),lst2)
  out <- df1[,c(2:4)]
  out <- out[order(out$GO_id, out$gene_id),]
  
  details = (allRes[allRes$topgoFisher.padjust < 0.05,])[,c(1,2)]
  details
  
  out1 <- merge(out, details, by.x="GO_id", by.y="GO.ID")
  
  write.table(out1, file = paste(traitname, "topgo_fisher_weight_tests_BP_q0.05sig_genes.txt", sep="_"), quote = FALSE, row.names = FALSE, sep = "\t")
  }
 
  myterms = allRes$GO.ID
  if (length(myterms) > 0){
    mygenes <- genesInTerm(myGOdata, myterms)
    for (j in 1:length(myterms))
    {
      myterm <- myterms[j]
      mygenesforterm <- mygenes[myterm][[1]]
      myfactor <- mygenesforterm %in% myInterestingGenes # find the genes that are in the list of genes of interest
      mygenesforterm2 <- mygenesforterm[myfactor == TRUE] 
      mygenesforterm3 <- data.frame(mygenesforterm2)
      mygenesforterm3$GO_id <- myterms[j]
      colnames(mygenesforterm3) <- c("gene_id", "GO_id")
      full <- merge(mygenesforterm3, anno, by="gene_id")
      my_list[[j]] <- full
    }
    
    # Write out single data frames
    lst2 <- lapply(my_list,function(x) cbind(rowname=rownames(x),x))
    df1 <- Reduce(function(x,y) merge(x,y,all=T),lst2)
    out <- df1[,c(2:4)]
    out <- out[order(out$GO_id, out$gene_id),]
    
    details = allRes[,c(1,2)]
    details
    
    out1 <- merge(out, details, by.x="GO_id", by.y="GO.ID")
    
    write.table(out1, file = paste(traitname, "topgo_fisher_weight_tests_BP_fullsig_genes.txt", sep="_"), quote = FALSE, row.names = FALSE, sep = "\t")
  }
}

