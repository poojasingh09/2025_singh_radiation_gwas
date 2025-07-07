

### ancestor matching
setwd("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/ancestor_matching/")

## but FIRST let's get all the good stuff we need

################ summarise all top hits - like you know

# all samples and run some heatmap and PCAs
library(ggplot2)
library(ggpubr)
theme_set(theme_pubr())
library("gridExtra")
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)
library(cowplot)
library(stringr)

setwd("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/snp_effect_overlaps/")

cs <- read.table("cs_top_99.9_and_pip01.txt", header=T)
ir <- read.table("ir_top_99.9_and_pip01.txt", header=T)
kd <- read.table("kd_top_99.9_and_pip01.txt", header=T)
tm <- read.table("tm_top_99.9_and_pip01.txt", header=T)
mb <- read.table("mb_top_99.9_and_pip01.txt", header=T)
fb <- read.table("fb_top_99.9_and_pip01.txt", header=T)
yf <- read.table("yf_top_99.9_and_pip01.txt", header=T)
rc <- read.table("rc_top_99.9_and_pip01.txt", header=T)
bb <- read.table("bb_top_99.9_and_pip01.txt", header=T)
mlb <- read.table("mlb_top_99.9_and_pip01.txt", header=T)
dlb <- read.table("dlb_top_99.9_and_pip01.txt", header=T)
vb <- read.table("vb_top_99.9_and_pip01.txt", header=T)
pc1 <- read.table("pc1_top_99.9_and_pip01.txt", header=T)
pc2 <- read.table("pc2_top_99.9_and_pip01.txt", header=T)


## wd
setwd("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/ancestor_matching/")


## prep files

mb$trait <- "mb"
fb$trait <- "fb"
yf$trait <- "yf"
rc$trait <- "rc"
bb$trait <- "bb"
mlb$trait <- "mlb"
dlb$trait <- "dlb"
vb$trait <- "vb"
pc1$trait <- "pc1"
pc2$trait <- "pc2"

cs$trait <- "cs"
kd$trait <- "kd"
tm$trait <- "tm"
ir$trait <- "ir"

all <- rbind(cs, ir, kd, tm, pc1, pc2,
                              vb, mlb, dlb,
                              fb, yf, rc, mb, bb)
dim(all)
write.table(all, "top_rwas_all_99.9_and_pip01_v2.txt", row.names=F, quote=F, sep="\t")


################ LETS FIRST LOOK AT allele frequencies of top hits across all samples
## be careful here because the frq file from all.filt2 vcf does not contain the ref and alt allele, you have to use the raw file for that
af <- read.table("../snp_effect_overlaps/top_rwas_all_99.9_and_pip01.frq.txt")
af1 <- af[,c(1,2,2,3,4,5,6)]
af <- af1
colnames(af) <- c("chr","start","end","n_alleles","n_chr","af1","af2")
af[c('ref', 'ref_af')] <- str_split_fixed(af$af1, ':', 2)
af[c('alt', 'alt_af')] <- str_split_fixed(af$af2, ':', 2)
af$ref_af <- as.numeric(af$ref_af)
af$alt_af <- as.numeric(af$alt_af)
af$id <- paste(af$chr, af$start, sep="-")
af$id <- gsub("chr", "",af$id)
head(af)

all_af <- merge(af, all, by="id")
all_af1 <- unique(all_af)
dim(all_af)
dim(all_af1)

## this comand selects the higher AF so the major allele
all_af1$major_af <-ifelse(all_af1$ref_af > all_af1$alt_af, all_af1$ref_af,ifelse(all_af1$ref_af < all_af1$alt_af, all_af1$alt_af, "No"))
all_af1$major_af <- as.numeric(all_af1$major_af)

## this comand selects the lower AF so the rarer allele: THIS IS THE GWAS EFFECT ALLELE FROM GEMMA SEE https://www.xzlab.org/software/GEMMAmanual.pdf
all_af1$minor_af <-ifelse(all_af1$ref_af < all_af1$alt_af, all_af1$ref_af,ifelse(all_af1$ref_af > all_af1$alt_af, all_af1$alt_af, "No"))
all_af1$minor_af <- as.numeric(all_af1$minor_af)

write.table(all_af1, "top_rwas_all_99.9_and_pip01.allthegoodstuff.txt", quote=F, row.names=F, sep="\t")
dim(all_af1)


#### read effect allele in 
af <- read.table("top_rwas_all_99.9_and_pip01.allthegoodstuff.txt", header=T, check.names=F)
#colnames(af) <- c("chr","start","end","n_alleles","n_chr","af1","af2")

#
# af[c('ref', 'ref_af')] <- str_split_fixed(af$af1, ':', 2)
# af[c('alt', 'alt_af')] <- str_split_fixed(af$af2, ':', 2)
# af$ref_af <- as.numeric(af$ref_af)
# af$alt_af <- as.numeric(af$alt_af)
# af$id <- paste(af$chr, af$start, sep="-")
# af$id <- gsub("chr", "",af$id)
# head(af)

#### read individuals genotypes in ancestors
all <- read.table("../snp_effect_overlaps/top_rwas_all_99.9_and_pip01.anc.GT.txt", header=T, check.names = FALSE) # this includes genotypes of some samples that didnt go into the RWAS but it ensured that this file had all the ancestors
head(all)
all$id <- paste(all$CHROM, all$POS, sep="-")
all$id <- gsub("chr", "",all$id)
head(all)
all1 <- all[c("id", '71001', '71014', '81032', '81035', '81342', '81343', '81344', '81345')] # '81035' remove, shit quality genome, tpharygalis

all_af <- merge(af, all1, by="id")
all_af1 <- unique(all_af)
dim(all_af)
dim(all_af1)
all_af1$special <- paste0(all_af1$chr,"_", all_af1$start.x)
all_af2 <- all_af1[!duplicated(all_af1$special),]
rownames(all_af2) <- all_af2$special 
dim(all_af2)
dim(all_af1)

## this comand selects the higher AF so the major allele
all_af1$major_af <-ifelse(all_af1$ref_af > all_af1$alt_af, all_af1$ref_af,ifelse(all_af1$ref_af < all_af1$alt_af, all_af1$alt_af, "No"))
all_af1$major_af <- as.numeric(all_af1$major_af)
all_af1$major_a <-ifelse(all_af1$ref_af > all_af1$alt_af, all_af1$ref,ifelse(all_af1$ref_af < all_af1$alt_af, all_af1$alt, "No"))

## this comand selects the lower AF so the rarer allele: THIS IS THE GWAS EFFECT ALLELE FROM GEMMA SEE https://www.xzlab.org/software/GEMMAmanual.pdf
all_af1$minor_af <-ifelse(all_af1$ref_af < all_af1$alt_af, all_af1$ref_af,ifelse(all_af1$ref_af > all_af1$alt_af, all_af1$alt_af, "No"))
all_af1$minor_af <- as.numeric(all_af1$minor_af)
all_af1$minor_a <-ifelse(all_af1$ref_af < all_af1$alt_af, all_af1$ref,ifelse(all_af1$ref_af > all_af1$alt_af, all_af1$alt, "No"))
head(all_af1)

######## ancestor matching

all_af1_a <- all_af1[all_af1$minor_a == all_af1$alt,] # is the minor allele alt
all_af1_r <- all_af1[all_af1$minor_a == all_af1$ref,] # is the minor allele ref
dim(all_af1_a)
dim(all_af1_r)
  


all_af1_a$gwas_allele_TP <- ifelse(all_af1_a$`71001` == "1" | all_af1_a$`71001` == "2" | 
                                     all_af1_a$`71014` == "1" | all_af1_a$`71014` == "2" | 
                                     all_af1_a$`81032` == "1" | all_af1_a$`81032` == "2",
                                   "yes" , "no")

all_af1_a$gwas_allele_AS <- ifelse(all_af1_a$`81342` == "1" | all_af1_a$`81342` == "2" | 
                                     all_af1_a$`81343` == "1" | all_af1_a$`81343` == "2" | 
                                     all_af1_a$`81344` == "1" | all_af1_a$`81344` == "2" | 
                                     all_af1_a$`81345` == "1" | all_af1_a$`81345` == "2" ,
                                   "yes" , "no")

all_af1_r$gwas_allele_TP <- ifelse(all_af1_r$`71001` == "1" | all_af1_r$`71001` == "0" | 
                                     all_af1_r$`71014` == "1" | all_af1_r$`71014` == "0" | 
                                     all_af1_r$`81032` == "1" | all_af1_r$`81032` == "0",
                                   "yes" , "no")

all_af1_r$gwas_allele_AS <- ifelse(all_af1_r$`81342` == "1" | all_af1_r$`81342` == "0" | 
                                     all_af1_r$`81343` == "1" | all_af1_r$`81343` == "0" | 
                                     all_af1_r$`81344` == "1" | all_af1_r$`81344` == "0" | 
                                     all_af1_r$`81345` == "1" | all_af1_r$`81345` == "0" ,
                                   "yes" , "no")


all_af2 <- rbind(all_af1_a, all_af1_r)
  
all_af2$gwas_allele_TP[is.na(all_af2$gwas_allele_TP)] <- "no"
all_af2$gwas_allele_TP[is.na(all_af2$gwas_allele_AS)] <- "no"
head(all_af2)
dim(all_af2)

write.table(all_af2, "ancestor_matching_results_99.9_pip01.txt", quote=F)

traits <- c("cs", "ir", "kd", "tm", "pc1", "pc2",
             "vb", "mlb", "dlb",
             "fb", "yf", "rc", "mb", "bb")

Res<-list()

for (trait in traits){
  all_af3 <- all_af2[all_af2$trait == trait,]
  res <- data.frame(matrix(ncol = 4, nrow = 5))
  colnames(res) <- c("trait", "source", "value", "percent")
  res$trait <- trait
  
  res$source <- c("ancestor_Tpharyngalis", "ancestor_Astappersi", "both", "neither", "total")
  res[1, 3] <-nrow(all_af3[all_af3$gwas_allele_TP == "yes" & all_af3$gwas_allele_AS == "no",])
  res[2, 3] <-nrow(all_af3[all_af3$gwas_allele_TP == "no" & all_af3$gwas_allele_AS == "yes",])
  res[3, 3] <- nrow(all_af3[all_af3$gwas_allele_TP == "yes" & all_af3$gwas_allele_AS == "yes",])
  res[4, 3] <-nrow(all_af3[all_af3$gwas_allele_TP == "no" & all_af3$gwas_allele_AS == "no",])
  res[5, 3] <- sum(res[c(1:4),3])
  res$percent <- (res$value/res[5,3])*100
  print(res)
  Res[[trait]]<- res

}


Res1 <- do.call(rbind, Res)

##plot it

library(ggplot2)

Res1$trait <- factor(Res1$trait, levels = unique(Res1$trait)) ## this prevents reodering of bars
Res1a <- Res1[Res1$source != "total",]




  
## stappersi is congolese and pharyngalis is nilotic
  ## set expectations


  counter=0
  res <- data.frame()
  iterations <- c(1:100)
  for (iter in iterations) {
    counter=counter+1
    data1 <- read.table("../../2024_bcftools_pnye3_dentition/ancestor_matching/all.SNPs.raw.ancestorsonly.random100k.vcf.GT.012", header=T, check.names = F)
    data2 <- data1[,c(3:10)]
    data2$id <- paste(data1$CHROM, data1$POS, sep="_")
    data <- data2[c("id", '71001', '71014', '81032', '81342', '81343', '81344', '81345')] # '81035' remove, shit quality genome, tpharygalis
    
    data1 <- data[sample(nrow(data), 200), ]

    data1$gwas_allele_TP <- ifelse(data1$`71001` == "1" | data1$`71001` == "2" | 
                                     data1$`71014` == "1" | data1$`71014` == "2" | 
                                     data1$`81032` == "1" | data1$`81032` == "2",
                                   "yes" , "no")
    
    data1$gwas_allele_AS <- ifelse(data1$`81342` == "1" | data1$`81342` == "2" | 
                                     data1$`81343` == "1" | data1$`81343` == "2" | 
                                     data1$`81344` == "1" | data1$`81344` == "2" | 
                                     data1$`81345` == "1" | data1$`81345` == "2" ,
                                   "yes" , "no")
    
    data1$gwas_allele_TP[is.na(data1$gwas_allele_TP)] <- "no"
    data1$gwas_allele_TP[is.na(data1$gwas_allele_AS)] <- "no"

    
    res[iter, 1] <-nrow(data1[data1$gwas_allele_TP == "yes" & data1$gwas_allele_AS == "no",])
    res[iter, 2] <-nrow(data1[data1$gwas_allele_TP == "no" & data1$gwas_allele_AS == "yes",])
    res[iter, 3] <- nrow(data1[data1$gwas_allele_TP == "yes" & data1$gwas_allele_AS == "yes",])
    res[iter, 4] <-nrow(data1[data1$gwas_allele_TP == "no" & data1$gwas_allele_AS == "no",])
    print(iter)
  }
  res
  #rownames(res) <- c("resample1", "resample2", "resample3", "resample4", "resample5", "resample6", "resample7", "resample8", "resample9", "resample10")
  colnames(res) <- c("ancestor_Tpharyngalis", "ancestor_Astappersi", "both", "neither")

  
  write.table(res, "ancestor_matching_resampling_table.txt", quote=F, sep="\t")
  res$t <- (res$ancestor_Tpharyngalis/200)*100
  res$a <- (res$ancestor_Astappersi/200)*100
  res$b <- (res$both/200)*100
  res$n <- (res$neither/200)*100
  
  
  quantile(res$a, 0.99)
  quantile(res$t, 0.99)
  quantile(res$b, 0.99)
  quantile(res$n, 0.99)
  
  a <- Res1[Res1$source == "ancestor_Astappersi",]
  dim(a[a$percent > quantile(res$a, 0.99),])
  
  t <- Res1[Res1$source == "ancestor_Tpharyngalis",]
  dim(t[t$percent > quantile(res$t, 0.99),])
  
  b <- Res1[Res1$source == "both",]
  dim(b[b$percent > quantile(res$b, 0.99),])
  
  n <- Res1[Res1$source == "neither",]
  dim(n[n$percent > quantile(res$n, 0.99),])
  
  a <- Res1[Res1$source == "ancestor_Astappersi",]
  dim(a[a$percent > quantile(res$a, 0.95),])
  
  t <- Res1[Res1$source == "ancestor_Tpharyngalis",]
  dim(t[t$percent > quantile(res$t, 0.95),])
  
  b <- Res1[Res1$source == "both",]
  dim(b[b$percent > quantile(res$b, 0.95),])
  
  n <- Res1[Res1$source == "neither",]
  dim(n[n$percent > quantile(res$n, 0.95),])
  
  

### effect size of admixture derived SNPs
  library(ggplot2)
  library(reshape2)
  library(forcats)  # For factor ordering
  library(dplyr)
  
a <- read.table("ancestor_matching_results_99.9_pip01.txt", header=T, check.names=F)
  
b <- a[!is.na(a$gwas_allele_TP) & !is.na(a$gwas_allele_AS) & a$gwas_allele_TP != a$gwas_allele_AS, ]
head(b)
dim(b)
c <- a[is.na(a$gwas_allele_TP) | is.na(a$gwas_allele_AS) | a$gwas_allele_TP == a$gwas_allele_AS, ]

# Define colors for traits
# Define colors for traits
set_colors <- c("CS" = "#4DBBD5", "IR" = "#4DBBD5", 
                "KD" = "#4DBBD5", "TM" = "#4DBBD5", 
                "PC1" = "#4DBBD5", "PC2" = "#4DBBD5", 
                "VB" = "#009E73", "MLB" = "#009E73", 
                "DLB" = "#009E73", "FB" = "#E64B35", 
                "YF" = "#E64B35", "RC" = "#E64B35", 
                "MB" = "#E64B35", "BF" = "#E64B35")

# Define trait order
trait_order <- c("CS", "IR", "KD", "TM", "PC1", "PC2",
                 "VB", "MLB", "DLB", "FB", "YF", "RC", "MB", "BF")

# Combine b and c datasets into one (with a 'dataset' column to distinguish them)
b$c_data <- "admixture"  # Update to 'admixture' for b
c$c_data <- "non-admixture"  # Update to 'non-admixture' for c

# Combine the data
combined_data <- rbind(
  melt(b, id.vars = "trait", measure.vars = "beta")[, c(1, 3)], 
  melt(c, id.vars = "trait", measure.vars = "beta")[, c(1, 3)]
)
combined_data$c_data <- c(rep("admixture", nrow(b)), rep("non-admixture", nrow(c)))

# Rename trait categories for combined dataset
combined_data$traits <- recode(combined_data$trait,  
                               "cs" = "CS",
                               "ir" = "IR",
                               "kd" = "KD",
                               "tm" = "TM",
                               "pc1" = "PC1",
                               "pc2" = "PC2",
                               "vb" = "VB",
                               "mlb" = "MLB",
                               "dlb" = "DLB",
                               "fb" = "FB",
                               "yf" = "YF",
                               "rc" = "RC",
                               "mb" = "MB",
                               "bb" = "BF")

# Ensure traits are factors with correct order
combined_data$traits <- factor(combined_data$traits, levels = trait_order)

# Create a new column combining the trait and admixture status
combined_data$trait_label <- interaction(combined_data$traits, combined_data$c_data, sep = " ")

# Ensure the trait_label respects the order of traits
combined_data$trait_label <- factor(combined_data$trait_label, levels = paste(rep(trait_order, each = 2), 
                                                                              c("admixture", "non-admixture"), sep = " "))

# Boxplot with color patterns and separate boxes for admixture and non-admixture
p_combined <- ggplot(combined_data, aes(x = trait_label, y = abs(value), fill = traits)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # Dodge to separate boxes for 'admixture' and 'non-admixture'
  scale_fill_manual(values = set_colors) +  # Use the set_colors for filling the traits
  theme_classic() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1),  # Rotate for readability
        axis.text.y = element_text(size = 18), 
        axis.title.y = element_text(size = 18),
        legend.position = "none") +  # Remove legend
  ylab(expression("SNP effect size |" ~ beta ~ "|")) +
  xlab("")  # No x-axis label as trait is already included in labels


svg("snp_effectsize_admixed_nonadmixed.svg")
# Display the plot
print(p_combined)
dev.off()

## check sig 
# Create an empty dataframe to store Wilcoxon test results
wilcox_results <- data.frame(trait = character(), p_value = numeric(), significance = character(), stringsAsFactors = FALSE)

# Loop through each trait and perform the Wilcoxon test
for (trait in trait_order) {
  # Subset data for admixture and non-admixture groups
  trait_data_admixture <- combined_data[combined_data$traits == trait & combined_data$c_data == "admixture", ]
  trait_data_non_admixture <- combined_data[combined_data$traits == trait & combined_data$c_data == "non-admixture", ]
  
  # Perform the Wilcoxon rank-sum test
  wilcox_test <- wilcox.test(trait_data_admixture$value, trait_data_non_admixture$value)
  
  # Store the results
  wilcox_results <- rbind(wilcox_results, 
                          data.frame(trait = trait, 
                                     p_value = wilcox_test$p.value, 
                                     significance = ifelse(wilcox_test$p.value < 0.05, "*", "ns")))
}

# Print the Wilcoxon test results
print(wilcox_results)

